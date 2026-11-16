# -*- coding: utf-8 -*-
# Place reads onto the consensus map, producing ra_*.tsv and fr_*.tsv.
#   Step: placement.

# -*- coding: utf-8 -*-
# Place reads onto the consensus map and detect flank-bearing reads directly.
#
#  The design:
#
#  1) Two 50 kb flanks are included in the reference.
#     The reference is [haplotype 1][haplotype 2][proximal 50k][distal 50k].
#     - A read carrying flank aligns to the flank reference, so the classification comes straight from the alignment.
#     - Because the flank portion has somewhere to align, the aligned fraction can be read at face value
#       (without the flank reference it is clipped, and the aligned fraction comes out low, 31-89%).
#
#  2) Read classes, decided from the alignment alone:
#       SPAN: touches both the proximal and the distal flank, and the array as well.
#       SINGLE_P: the proximal flank only.
#       SINGLE_D: the distal flank only.
#       ARRAY_ONLY: no flank.
#
#  3) There are two outputs:
#       ra_*.tsv, the read count and methylation per unit, for the map;
#       fr_*.tsv, the assignment and bases per read, for checking the flank reads.
#
#  Both haplotypes must be included together; with only one, the other haplotype's reads align to it as well.
#  The orientation is taken from the revcomp flag of CONSENSUS_SUM_G6.tsv and never recomputed.
#  minimap2 -y carries MM/ML across; without it every beta comes out as "-".
#  MAPQ is reported alongside: the two haplotypes are about 99.6% identical, so where they cannot be told apart MAPQ is 0.
import io
import os
import re
import sys
import gzip
import bisect
import subprocess
import collections
import pysam

W = os.path.expanduser("~/ralign")
os.makedirs(W, exist_ok=True)
D = os.path.expanduser("~/cmp4d")
R = os.path.expanduser("~/s5s_tuc/ref")
NL = "\n"
TH = os.environ.get("TH", "8")
TAG = os.environ.get("TAG", "x")
CMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
PROX = "PROX50k"
DIST = "DIST50k"
# Deciding whether a read carries flank: length alone is not enough.
#   Measured: requiring only 300 bp of alignment called 4,041 of COLO829-BL's 4,229 reads (96%)
#   as carrying flank. These are false positives aligning to Alu and L1 elements inside the 50 kb flank.
#   So the distance to the junction is examined as well.
#     FLANK_proximal_50000: the junction is at its end (CHM13 227,694,885-227,744,884).
#     FLANK_distal_50000: the junction is at its start (CHM13 228,025,900-228,075,899).
#   A genuine flank read continues into the array, so it must align at the junction end.
#   The distance threshold was set by measurement, over 13 pilot individuals against the canonical table.
#     The distribution of junction distances is clearly bimodal:
#       proximal: 327 reads at 0 with a tail to 44, then a gap, then 168 (9 reads) and 23,641-33,533.
#       distal: 307 reads at 0 with a tail to 44, then a gap, then 5,652-7,503.
#     All 18 reads the canonical table calls SPAN have a junction distance of 17 or less on both sides.
#     The reads at 168 align to a repeat block of about 5,922 bp inside the proximal flank; eight of them
#     match in length and distance as well (5,916-5,924 bp, junction distance exactly 168), which makes them plainly false positives.
#     At a threshold of 2000 the SPAN count is overcalled at 27; at 50 it is 19,
#       and the canonical SPAN reads are recovered 18 of 18 either way.
MINFL = 1000         # the minimum alignment length against the flank reference
JXWIN = 50           # the alignment must reach within this distance of the junction, the boundary measured in the bimodal distribution


def rc(s):
    return s.translate(CMP)[::-1]


def rdfa(p):
    return "".join(x.strip() for x in
                   io.open(p, encoding="utf8").read().split(NL)[1:])


FL_P = rdfa(os.path.join(R, "FLANK_proximal_50000.fa"))
FL_D = rdfa(os.path.join(R, "FLANK_distal_50000.fa"))


def find_reads(s):
    """Read source: the ONT fastq first, and if absent the HiFi FASTA of TUC.

    The two are not the same data and are recorded separately:
      ONT   `*.hits.fastq.gz` / `*.array.fastq.gz`
            - long (median tens to hundreds of kb) -> flank-to-flank spanning is possible
            - MM/ML in the header, so methylation can be reported
      HIFI  `~/s5s_tuc/cand/<sample>.fa`  (output of the TUC pipeline)
            - PacBio HiFi, median 18.5 kb (about 8 units) -> spanning is almost impossible in principle
            - FASTA, so no MM/ML -> beta is "-" throughout
            - higher base accuracy, which suits sequence comparison better
    Only 7 individuals have no data at all in either source;
      the others have at least one of the two.
    """
    for p in ("~/s5s_ont/hits/%s.hits.fastq.gz" % s,
              "~/s5s_new/hits/%s.hits.fastq.gz" % s,
              "~/kgp_ont/out/%s.array.fastq.gz" % s,
              "~/giab_hits/%s.hits.fastq.gz" % s,
              "~/cancer5s/hits/%s.hits.fastq.gz" % s):
        q = os.path.expanduser(p)
        if os.path.exists(q) and os.path.getsize(q) > 1000:
            return q, "ONT"
    q = os.path.expanduser("~/s5s_tuc/cand/%s.fa" % s)
    if os.path.exists(q) and os.path.getsize(q) > 1000:
        return q, "HIFI"
    return None, None


def sniff_platform(path, guess):
    """The platform is decided from the read name, not from the directory.

    A directory such as `~/cancer5s/hits/` can hold PacBio HiFi reads (m84039_...).
    Aligning those with map-ont lowers the alignment identity:
    measured 69% for COLO829, against 97-99% in the normal individuals.
    PacBio names carry `m<digits>_<date>_...` or `/ccs`.
    """
    try:
        op = gzip.open if path.endswith(".gz") else io.open
        with op(path, "rt") as fh:
            for L in fh:
                L = L.strip()
                if not L or L[0] not in "@>":
                    continue
                nm = L[1:].split()[0]
                if "|" in nm:
                    nm = nm.split("|", 1)[1]
                if re.match(r"^m\d+_\d+", nm) or "/ccs" in nm:
                    return "HIFI"
                return "ONT" if guess == "ONT" else guess
    except Exception:
        pass
    return guess


SEQF = {}
for L in io.open(os.path.expanduser("~/hap_seq_node.tsv"), encoding="utf8"):
    f = L.rstrip(NL).split("\t")
    if len(f) >= 3 and f[0] != "sample":
        SEQF[(f[0], f[1])] = f[2]
RCF = {}
UNI = {}
for i, L in enumerate(io.open(os.path.expanduser("~/CONSENSUS_SUM_G6.tsv"), encoding="utf8")):
    if i == 0:
        continue
    f = L.rstrip(NL).split("\t")
    if len(f) >= 6:
        RCF[(f[0], f[1])] = (f[5] == "Y")
        UNI[(f[0], f[1])] = int(f[4])

UB = collections.defaultdict(list)
for i, L in enumerate(io.open(os.path.expanduser("~/CONSENSUS_MAP_G6.tsv"), encoding="utf8")):
    if i == 0:
        continue
    f = L.rstrip(NL).split("\t")
    if len(f) < 8 or f[4] == "PARTIAL_PROX":
        continue
    UB[(f[0], f[1])].append((int(f[3]), int(f[5]) - 1, int(f[6])))

OUT = io.open(os.path.join(W, "ra_%s.tsv" % TAG), "a", encoding="utf8", newline=NL)
FR = io.open(os.path.join(W, "fr_%s.tsv" % TAG), "a", encoding="utf8", newline=NL)
smps = [x.strip() for x in io.open(sys.argv[1], encoding="utf8") if x.strip()]
DONE = os.path.join(W, "done")
os.makedirs(DONE, exist_ok=True)

for smp in smps:
    mk = os.path.join(DONE, smp + ".done")
    if os.path.exists(mk):
        continue
    q, rsrc = find_reads(smp)
    if q:
        rsrc = sniff_platform(q, rsrc)      # keyed on the read name rather than the directory
    slots = sorted(k[1] for k in SEQF if k[0] == smp)
    if not q or not slots:
        sys.stderr.write("  건너뜀 %s (read %s · 슬롯 %d)\n" % (smp, bool(q), len(slots)))
        continue
    ref = os.path.join(W, "%s.ref.fa" % smp)
    names = {}
    reflen = {}
    with io.open(ref, "w", newline=NL) as g:
        for sl in slots:
            p = os.path.join(D, SEQF[(smp, sl)])
            if not os.path.exists(p):
                continue
            sq = rdfa(p)
            if RCF.get((smp, sl)):
                sq = rc(sq)
            nm = "%s_h%s" % (smp, sl)
            names[nm] = sl
            reflen[nm] = len(sq)
            g.write(">" + nm + NL + sq + NL)
        # the flank reference, so that a read's flank content is decided directly by alignment
        g.write(">" + PROX + NL + FL_P + NL)
        g.write(">" + DIST + NL + FL_D + NL)
    if not names:
        continue
    bam = os.path.join(W, "%s.bam" % smp)
    # the preset matched to the platform; using map-ont on HiFi throws away its accuracy
    preset = "map-ont" if rsrc == "ONT" else "map-hifi"
    p1 = subprocess.Popen(["minimap2", "-ax", preset, "-y", "-t", TH,
                           "--secondary=no", ref, q],
                          stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    p2 = subprocess.Popen(["samtools", "sort", "-@", "2", "-o", bam],
                          stdin=p1.stdout, stderr=subprocess.DEVNULL)
    p1.stdout.close()
    p2.wait()
    subprocess.run(["samtools", "index", bam], check=True, capture_output=True)
    bf = pysam.AlignmentFile(bam, "rb")

    # ---- 1) per-read classification: which references it touched ----
    HIT = collections.defaultdict(dict)      # read id -> reference name -> total aligned length
    JX = collections.defaultdict(dict)       # read id -> reference name -> minimum distance to the junction
    BEST = {}                                # read id -> its primary alignment on a haplotype
    NSUP = collections.Counter()
    LP, LD = len(FL_P), len(FL_D)
    for a in bf.fetch(until_eof=True):
        if a.is_unmapped or a.is_secondary:
            continue
        rid = a.query_name
        rn = a.reference_name
        HIT[rid][rn] = HIT[rid].get(rn, 0) + a.query_alignment_length
        if rn == PROX:                       # the junction is at the end
            d = LP - a.reference_end
        elif rn == DIST:                     # the junction is at the start
            d = a.reference_start
        else:
            d = None
        if d is not None:
            JX[rid][rn] = min(JX[rid].get(rn, 10 ** 9), max(0, d))
        if a.is_supplementary:
            NSUP[rid] += 1
        elif rn in names:
            BEST[rid] = a
    KLS = {}
    for rid, h in HIT.items():
        hp = (h.get(PROX, 0) >= MINFL and JX[rid].get(PROX, 10 ** 9) <= JXWIN)
        hd = (h.get(DIST, 0) >= MINFL and JX[rid].get(DIST, 10 ** 9) <= JXWIN)
        KLS[rid] = ("SPAN" if (hp and hd) else "SINGLE_P" if hp
                    else "SINGLE_D" if hd else "ARRAY_ONLY")
    for rid, a in BEST.items():
        cl = KLS.get(rid, "ARRAY_ONLY")
        nmt = a.get_tag("NM") if a.has_tag("NM") else -1
        al = a.query_alignment_length
        idt = 100.0 * (1.0 - nmt / float(al)) if (al and nmt >= 0) else -1
        rl = a.infer_read_length() or len(a.query_sequence or "")
        # the aligned fraction is taken including supplementary alignments, because the flank portion aligns separately
        tot_al = sum(HIT[rid].values())
        sl = names[a.reference_name]
        rcv = 100.0 * (a.reference_end - a.reference_start) / reflen[a.reference_name]
        jp = JX[rid].get(PROX, -1)
        jd = JX[rid].get(DIST, -1)
        FR.write("%s\t%s\t%s\t%s\t%s\t%d\t%d\t%.3f\t%.1f\t%.1f\t%.1f\t%d\t%d"
                 "\t%d\t%d\t%d\t%d%s" %
                 (smp, rid[:12], rsrc, cl, sl, UNI.get((smp, sl), 0),
                  a.mapping_quality, idt,
                  100.0 * al / rl if rl else 0,             # the primary aligned fraction
                  100.0 * tot_al / rl if rl else 0,         # the aligned fraction including supplementary alignments
                  rcv, NSUP[rid], rl,
                  HIT[rid].get(PROX, 0), jp,                # proximal flank length and junction distance
                  HIT[rid].get(DIST, 0), jd, NL))           # distal flank length and junction distance

    # ---- 2) reads and methylation per unit ----
    for nm, sl in names.items():
        ub = UB.get((smp, sl), [])
        if not ub:
            continue
        span = collections.Counter()
        over = collections.Counter()
        spq = collections.Counter()
        spf = collections.Counter()          # reads carrying flank only
        msum = collections.defaultdict(lambda: [0.0, 0])
        asum = collections.defaultdict(lambda: [0.0, 0])
        qsum = collections.defaultdict(lambda: [0.0, 0])
        ntot = 0
        nq0 = 0
        for a in bf.fetch(nm):
            if a.is_unmapped or a.is_secondary or a.is_supplementary:
                continue
            ntot += 1
            hiq = (a.mapping_quality >= 20)
            hasfl = KLS.get(a.query_name, "ARRAY_ONLY") != "ARRAY_ONLY"
            if a.mapping_quality == 0:
                nq0 += 1
            lo, hi = a.reference_start, a.reference_end
            try:
                mb = a.modified_bases
            except Exception:
                mb = None
            rmod = []
            if mb:
                q2r = dict(a.get_aligned_pairs(matches_only=True))
                for key, lst in mb.items():
                    if len(key) < 3 or key[2] not in ("m", 27551):
                        continue
                    for qp, ql in lst:
                        rp = q2r.get(qp)
                        if rp is not None:
                            rmod.append((rp, 1.0 if ql >= 128 else 0.0))
                rmod.sort()
            # Scanning the whole modification list for every unit is O(units x calls) and took 5.3 minutes per individual
            #   (154 units x 10,000 calls x 3,000 reads = 4.6 billion operations).
            #   Both lists are in coordinate order, so a binary search cuts out just the interval needed.
            RP = [x[0] for x in rmod]
            for pos, s0, e0 in ub:
                if hi <= s0 or lo >= e0:
                    continue
                full = (lo <= s0 and hi >= e0)
                over[pos] += 1
                if full:
                    span[pos] += 1
                    if hiq:
                        spq[pos] += 1
                    if hasfl:
                        spf[pos] += 1
                if not RP:
                    continue
                i0 = bisect.bisect_left(RP, s0)
                i1 = bisect.bisect_left(RP, e0)
                if i1 <= i0:
                    continue
                sv = 0.0
                for j in range(i0, i1):
                    sv += rmod[j][1]
                nv = i1 - i0
                asum[pos][0] += sv
                asum[pos][1] += nv
                if full:
                    msum[pos][0] += sv
                    msum[pos][1] += nv
                    if hiq:
                        qsum[pos][0] += sv
                        qsum[pos][1] += nv
        for pos, s0, e0 in ub:
            bb = msum[pos]
            aa = asum[pos]
            qq = qsum[pos]
            OUT.write("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%d\t%d\t%d%s" %
                      (smp, sl, rsrc, pos, span[pos], over[pos] - span[pos], ntot,
                       ("%.4f" % (bb[0] / bb[1])) if bb[1] else "-",
                       ("%.4f" % (aa[0] / aa[1])) if aa[1] else "-",
                       bb[1], spq[pos],
                       ("%.4f" % (qq[0] / qq[1])) if qq[1] else "-",
                       qq[1], nq0, spf[pos], NL))
    bf.close()
    for x in (bam, bam + ".bai", ref):
        try:
            os.remove(x)
        except OSError:
            pass
    io.open(mk, "w").write("ok")
    OUT.flush()
    FR.flush()
    sys.stderr.write("  %s [%s] 완료 (read %d · flank 보유 %d · SPAN %d)\n" %
                     (smp, rsrc, len(HIT),
                      sum(1 for v in KLS.values() if v != "ARRAY_ONLY"),
                      sum(1 for v in KLS.values() if v == "SPAN")))
    sys.stderr.flush()
OUT.close()
FR.close()
print("done")
