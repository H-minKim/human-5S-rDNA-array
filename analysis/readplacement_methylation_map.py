# -*- coding: utf-8 -*-
# A two-layer methylation map: the aggregate, and the per-read values.
#   Step: placement and methylation.
#   Method document : analysis/METHOD_CONSENSUS_MAP_G8_REBUILD.md

# -*- coding: utf-8 -*-
# Build the methylation map by aligning array reads to the haplotype references, in two layers.
#
#  The reference is that individual's array sequence per haplotype:
#     tier 1 is the array interval of an array-spanning read;
#     tier 2 is a hifiasm reassembly of the array, verified at 100.0000% identity against T2T Q100 in HG002.
#
#  Both haplotypes are put in one reference and aligned together,
#    so that each read goes to its own haplotype: the primary alignment is the assignment.
#    Aligning to one at a time forces the other haplotype's reads onto it as well.
#
#  No MAPQ filter. In a tandem array, -q 20 erases the array window entirely
#    (measured: 63% of reads were misclassified as NO_ARRAY). Primary alignments only, with -F 0x900.
#
#  MM/ML is read through pysam's modified_bases, which corrects the coordinates itself
#    for a reverse alignment. Parsing the tag by hand makes the strand easy to get wrong.
#
#  The two layers:
#    1) the aggregate layer methmap_agg, with the mean beta over reads and its standard deviation per haplotype and unit;
#    2) the per-read layer methmap_read, with the beta per unit for every single read, which shows the methylation as it came off the read.
import io
import os
import sys
import gzip
import subprocess
import collections
import math

W = os.path.expanduser("~/methmap")
os.makedirs(W, exist_ok=True)
GENE = os.path.expanduser("~/s5s_tuc/ref/5S_gene_query.fa")
TH = os.environ.get("TH", "6")
TAG = os.environ.get("TAG", "x")
NL = "\n"
import pysam


def read_fa(p):
    nm = None
    buf = []
    out = {}
    for L in io.open(p, encoding="utf8"):
        if L[0] == ">":
            if nm:
                out[nm] = "".join(buf)
            nm = L[1:].strip()
            buf = []
        else:
            buf.append(L.strip())
    if nm:
        out[nm] = "".join(buf)
    return out


def unit_bounds(seq, tag):
    """cut the reference array into units -> [(start,end), ...]"""
    fa = os.path.join(W, "b_%s.fa" % tag)
    io.open(fa, "w", newline=NL).write(">x" + NL + seq + NL)
    db = os.path.join(W, "bd_%s" % tag)
    subprocess.run(["makeblastdb", "-in", fa, "-dbtype", "nucl", "-out", db],
                   check=True, capture_output=True)
    r = subprocess.run(["blastn", "-task", "blastn", "-query", GENE, "-db", db,
                        "-outfmt", "6 sstart send", "-dust", "no", "-word_size", "9",
                        "-evalue", "1e-5", "-max_target_seqs", "100000000",
                        "-num_threads", TH], capture_output=True, text=True, check=True)
    st = set()
    for L in r.stdout.splitlines():
        f = L.split("\t")
        if len(f) >= 2:
            st.add(min(int(f[0]), int(f[1])) - 1)
    K = sorted(st)
    for e in (".ndb", ".nhr", ".nin", ".not", ".nsq", ".ntf", ".nto"):
        try:
            os.remove(db + e)
        except OSError:
            pass
    os.remove(fa)
    return [(K[i], K[i + 1]) for i in range(len(K) - 1)
            if 2100 <= K[i + 1] - K[i] <= 2400]


def find_reads(smp):
    for p in (os.path.expanduser("~/s5s_ont/hits/%s.hits.fastq.gz" % smp),
              os.path.expanduser("~/kgp_ont/out/%s.array.fastq.gz" % smp),
              os.path.expanduser("~/giab_hits/%s.hits.fastq.gz" % smp),
              os.path.expanduser("~/cancer5s/hits/%s.hits.fastq.gz" % smp)):
        if os.path.exists(p):
            return p
    return None


AGG = io.open(os.path.join(W, "methmap_agg_%s.tsv" % TAG), "a",
              encoding="utf8", newline=NL)
RDL = io.open(os.path.join(W, "methmap_read_%s.tsv" % TAG), "a",
              encoding="utf8", newline=NL)

for spec in sys.argv[1:]:
    # spec = "sample:ref1.fa[:ref2.fa]"; the references are the array sequences, one per haplotype
    parts = spec.split(":")
    smp = parts[0]
    refs = [p for p in parts[1:] if os.path.exists(p)]
    if not refs:
        sys.stderr.write("  참조 없음 %s\n" % smp)
        continue
    q = find_reads(smp)
    if not q:
        sys.stderr.write("  read 없음 %s\n" % smp)
        continue

    # the two haplotypes are merged into one reference
    ref = os.path.join(W, "ref_%s.fa" % smp)
    UB = {}
    with io.open(ref, "w", newline=NL) as g:
        for hi, rp in enumerate(refs):
            d = read_fa(rp)
            nm0, sq = list(d.items())[0]
            name = "%s_hap%d" % (smp, hi + 1)
            g.write(u">%s\n%s\n" % (name, sq))
            UB[name] = (unit_bounds(sq, "%s_%d" % (smp, hi + 1)), len(sq), nm0)
    for name, (ub, L, nm0) in UB.items():
        sys.stderr.write("  참조 %s : %d bp · 유닛 %d  [%s]\n" % (name, L, len(ub), nm0[:50]))

    bam = os.path.join(W, "%s.bam" % smp)
    # -y moves MM/ML from the fastq comment into SAM tags, confirmed on minimap2 2.26
    p1 = subprocess.Popen(["minimap2", "-y", "-ax", "map-ont", "-t", TH,
                           "--secondary=no", ref, q], stdout=subprocess.PIPE,
                          stderr=subprocess.DEVNULL)
    p2 = subprocess.Popen(["samtools", "sort", "-@", "2", "-o", bam],
                          stdin=p1.stdout, stderr=subprocess.DEVNULL)
    p1.stdout.close()
    p2.wait()
    subprocess.run(["samtools", "index", bam], check=True, capture_output=True)

    bf = pysam.AlignmentFile(bam, "rb")
    nr = collections.Counter()
    per = collections.defaultdict(list)   # (hap, unit) -> [beta...]
    for a in bf.fetch(until_eof=False):
        if a.is_unmapped or a.is_secondary or a.is_supplementary:
            continue
        hap = a.reference_name
        nr[hap] += 1
        mb = None
        try:
            mb = a.modified_bases
        except Exception:
            mb = None
        if not mb:
            continue
        # 5mC only (C+m); the value is a probability on 0-255
        calls = []
        for k, v in mb.items():
            if len(k) >= 3 and (k[2] == "m" or k[2] == 77):
                calls.extend(v)
        if not calls:
            continue
        # read coordinate -> reference coordinate
        r2t = dict(a.get_aligned_pairs(matches_only=True))
        ub = UB[hap][0]
        bypos = collections.defaultdict(list)
        for rp, ql in calls:
            tp = r2t.get(rp)
            if tp is None:
                continue
            bypos[tp].append(ql / 255.0)
        if not bypos:
            continue
        # collected per unit
        import bisect
        starts = [x[0] for x in ub]
        acc = collections.defaultdict(list)
        for tp, vs in bypos.items():
            i = bisect.bisect_right(starts, tp) - 1
            if i < 0 or i >= len(ub) or tp >= ub[i][1]:
                continue
            acc[i].extend(vs)
        for ui, vs in acc.items():
            b = sum(vs) / len(vs)
            per[(hap, ui)].append(b)
            RDL.write("%s\t%s\t%s\t%d\t%d\t%d\t%.4f\t%d%s"
                      % (smp, hap, a.query_name, ui + 1, ub[ui][0], ub[ui][1],
                         b, len(vs), NL))
    bf.close()
    RDL.flush()
    for hap in sorted(UB):
        ub, L, nm0 = UB[hap]
        for ui in range(len(ub)):
            vs = per.get((hap, ui), [])
            if not vs:
                AGG.write("%s\t%s\t%d\t%d\t%d\t0\t-\t-%s"
                          % (smp, hap, ui + 1, ub[ui][0], ub[ui][1], NL))
                continue
            m = sum(vs) / len(vs)
            sd = math.sqrt(sum((x - m) ** 2 for x in vs) / len(vs)) if len(vs) > 1 else 0.0
            AGG.write("%s\t%s\t%d\t%d\t%d\t%d\t%.4f\t%.4f%s"
                      % (smp, hap, ui + 1, ub[ui][0], ub[ui][1], len(vs), m, sd, NL))
        sys.stderr.write("    %s : 붙은 read %d · 메틸 있는 유닛 %d/%d\n"
                         % (hap, nr[hap],
                            sum(1 for ui in range(len(ub)) if per.get((hap, ui))),
                            len(ub)))
    AGG.flush()
    sys.stderr.flush()
    os.remove(ref)
AGG.close()
RDL.close()
print("done")
