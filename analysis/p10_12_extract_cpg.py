# -*- coding: utf-8 -*-
# P10 batch BAM -> per-CpG table (5S and D4Z4 4q/10q only): read x unit x CHM13 unit coordinate c x methylation x ML.
# The unit frame, the landmarks and the polarity rule are identical to p10_08_extract_frame.py; only the output differs.
# Purpose: compare the transition at the read level (the 5S distal end against the D4Z4 proximal end), the sparsity of CpG sites, and the ML confidence.
import sys, io, os, glob, gzip, hashlib, bisect, collections, statistics
import pysam

def _need(_v):
    """No defaults are set: the script stops if a variable is missing.
    The roots of the data this script reads are listed in the
    "Running this code" section of README.md."""
    import os as _o, sys as _s
    _x = _o.environ.get(_v)
    if not _x:
        _s.exit(
            "environment variable %s is not set; it is required by %s"
            % (_v, _o.path.basename(__file__))
            + "\n(a module stops the script that imports it)."
            + "\nThe deposited copy has no default; see README.md"
            + " section \"Running this code\"."
        )
    return _x


P = os.path.join(_need("RDNA_WORK_ROOT"), "P10")
OUTD = os.path.join(_need("RDNA_WORK_ROOT"), "out", "P10")
K, BIN = 31, 40
FAM_UNIT = {"5S": "5S", "D4Z4_4q": "D4Z4_4q", "D4Z4_10q": "D4Z4_10q", "DXZ4": "DXZ4",
            "RNU2": "RNU2", "SST1": "SST1_a"}
# The reference unit length is the measured period of the CHM13 unit sequence (the most frequent offset of 31-mer self-matches).
#   Using the published lengths (RNU2 6100, D4Z4 3300, DXZ4 3000, SST1 2523) truncates the end of the longer units
#   (for RNU2 this cut 71 bp, which removed the upstream promoter region of the next U2 copy).
ULEN = {"5S": 2240, "D4Z4_4q": 3298, "D4Z4_10q": 3323, "DXZ4": 2972, "RNU2": 6171, "SST1": 2527}
want = sys.argv[1] if len(sys.argv) > 1 else "ALL"
tag = sys.argv[2] if len(sys.argv) > 2 else "ALL"


def rc(s):
    return s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


UN, name = {}, None
for L in io.open(P + "/data/P10_UNITS_CHM13.fa", encoding="utf8"):
    if L[0] == ">":
        name = L[1:].split("_unit2")[0]
        UN[name] = []
    else:
        UN[name].append(L.strip().upper())
UN = {k: "".join(v) for k, v in UN.items()}

KM = {}                                        # fam -> {kmer: (strand, i)}
for fam, un in FAM_UNIT.items():
    L = ULEN[fam]
    u = UN[un][:L]
    circ = u + u[:K - 1]                       # the unit is circular
    cnt = collections.Counter(circ[i:i + K] for i in range(L))
    d = {}
    for i in range(L):
        k = circ[i:i + K]
        if cnt[k] != 1 or "N" in k:
            continue
        r = rc(k)
        if r in cnt:
            continue
        d[k] = ("+", i)
        d[r] = ("-", i)
    KM[fam] = d
    print("앵커 %-9s 유닛 %d bp · 유일 31-mer %d" % (fam, L, len(d) // 2), flush=True)

REG = {}
for i, L in enumerate(io.open(P + "/out/P10_SEGMENT_REGISTRY.tsv", encoding="utf8")):
    if i == 0:
        continue
    f = L.rstrip("\n").split("\t")
    REG[(f[0], f[1], f[2])] = (f[15], f[16])        # final_label, status
POL = {}
for i, L in enumerate(io.open(P + "/out/P10_PLACED_CLEAN.tsv", encoding="utf8")):
    if i == 0:
        continue
    f = L.rstrip("\n").split("\t")
    POL[(f[0], f[3])] = f[6]


def landmarks(seq, fam):
    """Strand is decided per landmark (one piece can hold a sub-array in the opposite orientation, as in SST1).
    Returns : list of unit intervals [(a, b, strand)] (piece coordinates, a<=rp<b), and the fwd/rev hit counts"""
    d, L = KM[fam], ULEN[fam]
    est = {"+": [], "-": []}
    nf = nr = 0
    for j in range(0, len(seq) - K + 1):
        h = d.get(seq[j:j + K])
        if h is None:
            continue
        if h[0] == "+":
            nf += 1
            est["+"].append(j - h[1])
        else:
            nr += 1
            est["-"].append(j + h[1] + K - 1)
    units = []
    for st in "+-":
        v = sorted(est[st])
        lm, cur = [], []
        for x in v:
            if cur and x - cur[-1] > 40:
                if len(cur) >= 20:
                    lm.append(int(statistics.median(cur)))
                cur = []
            cur.append(x)
        if len(cur) >= 20:
            lm.append(int(statistics.median(cur)))
        for k in range(len(lm) - 1):
            g = lm[k + 1] - lm[k]
            if not (0.8 * L <= g <= 1.2 * L):
                continue
            if st == "+":
                units.append((lm[k], lm[k + 1], "+"))       # c = rp - a
            else:
                units.append((lm[k] + 1, lm[k + 1] + 1, "-"))   # c = (b-1) - rp
    units.sort()
    # overlapping units (boundaries in opposite orientation) are discarded
    clean = []
    for u in units:
        if clean and u[0] < clean[-1][1]:
            clean.pop()
            continue
        clean.append(u)
    return clean, nf, nr


def calls(a):
    mm = a.get_tag("MM") if a.has_tag("MM") else (a.get_tag("Mm") if a.has_tag("Mm") else None)
    if mm is None:
        return {}
    ml = list(a.get_tag("ML") if a.has_tag("ML") else a.get_tag("Ml"))
    seq = a.query_sequence
    if not seq:
        return {}
    orig = rc(seq) if a.is_reverse else seq
    o, base, idx = {}, 0, None
    for blk in [b for b in mm.split(";") if b]:
        f = blk.split(",")
        dd = [int(x) for x in f[1:] if x != ""]
        if f[0].startswith("C+m"):
            if idx is None:
                idx = [i for i, c in enumerate(orig) if c == "C"]
            j = -1
            for t, x in enumerate(dd):
                j += x + 1
                if j >= len(idx):
                    break
                op = idx[j]
                sp = (len(seq) - 1 - op) if a.is_reverse else op
                o[sp] = ml[base + t] if base + t < len(ml) else 0
        base += len(dd)
    return o


import gzip as _gz
FO = _gz.open("%s/P10_CPG_%s.tsv.gz" % (OUTD, tag), "wt", encoding="utf8", newline="\n")
FO.write("sample\thap\tfamily\tseg_md5\tunit_strand\tpolarity_used\tn_units\tread\tread_rev\t"
         "unit_idx\tunit_from_distal\tc\tmeth\tml\n")
FS = io.open("%s/P10_CPG_SEGS_%s.tsv" % (OUTD, tag), "w", encoding="utf8", newline="\n")
FS.write("sample\thap\tfamily\tseg_md5\tseg_len\tseg_strand\tfwd_hits\trev_hits\tpolarity_flank\t"
         "polarity_used\tn_units\tspacing_med\tspacing_iqr\n")
seen = set()
folders = sorted(glob.glob(P + "/data/place/*"))
for fd in folders:
    smp = os.path.basename(fd).split("_")[0]
    if want != "ALL" and smp not in want.split(","):
        continue
    ref, bam = fd + "/ref.fa", fd + "/placed.bam"
    if not (os.path.exists(ref) and os.path.exists(bam)):
        continue
    SEG = {}
    nm, buf = None, []
    for L in io.open(ref):
        if L[0] == ">":
            if nm:
                SEG[nm] = "".join(buf)
            nm, buf = L[1:].strip(), []
        else:
            buf.append(L.strip())
    if nm:
        SEG[nm] = "".join(buf)
    INFO = {}
    for nm, s in SEG.items():
        hap = nm.split("_")[2]
        md = hashlib.md5(s.encode()).hexdigest()[:12]
        lab, status = REG.get((smp, hap, md), ("?", "?"))
        if status != "ok" or lab not in ULEN or lab not in ("5S", "D4Z4_4q", "D4Z4_10q"):
            continue
        units, nf, nr = landmarks(s.upper(), lab)
        if len(units) < 2:
            continue
        sp = sorted(b - a for a, b, _s in units)
        # A block's orientation is set when at least 95% of the hits agree; D4Z4_4q has 0.5% of hits in the other orientation and was otherwise classified as mixed.
        tot = float(nf + nr)
        if tot and max(nf, nr) / tot >= 0.95:
            st = "+" if nf >= nr else "-"
        else:
            st = "mix"
        pol = POL.get((smp, md), "?")
        # Polarity: the flank call wins. If the flank is '?' and the block is in one orientation, fill it in from the unit strand.
        #        (In the pilot, the unit strand equalled the flank polarity in 24 of 24 single-orientation blocks.)
        if pol in ("+", "-"):
            pol_used = pol if (st == "mix" or st == pol) else "CONFLICT"
        else:
            pol_used = st if st in ("+", "-") else "NA"
        FS.write("%s\t%s\t%s\t%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%d-%d\n"
                 % (smp, hap, lab, md, len(s), st, nf, nr, pol, pol_used, len(units),
                    statistics.median(sp), sp[len(sp) // 4], sp[3 * len(sp) // 4]))
        INFO[nm] = (hap, lab, md, st, units, pol, pol_used)
    if not INFO:
        continue
    bf = pysam.AlignmentFile(bam)
    for a in bf.fetch():
        if a.is_secondary or a.is_supplementary or a.is_unmapped:
            continue
        inf = INFO.get(a.reference_name)
        if inf is None:
            continue
        hap, lab, md, st, units, pol, pol_used = inf
        key = (smp, md, a.query_name)
        if key in seen:
            continue
        seen.add(key)
        cl = calls(a)
        if not cl:
            continue
        L = ULEN[lab]
        starts = [u[0] for u in units]
        pairs = dict((q, r) for q, r in a.get_aligned_pairs(matches_only=True))
        per = collections.defaultdict(list)
        for q, v in cl.items():
            rp = pairs.get(q)
            if rp is None:
                continue
            u = bisect.bisect_right(starts, rp) - 1
            if u < 0:
                continue
            ua, ub, us = units[u]
            if rp >= ub:
                continue
            c = (rp - ua) if us == "+" else (ub - 1 - rp)
            if c < 0 or c >= L:
                continue
            per[u].append((c, 1 if v >= 128 else 0, v))
        nun = len(units)
        rv = 1 if a.is_reverse else 0
        for u in sorted(per):
            if len(per[u]) < 5:
                continue
            if pol_used == "+":
                ufd = nun - u
            elif pol_used == "-":
                ufd = u + 1
            else:
                ufd = "NA"
            for c, m, ml in sorted(per[u]):
                FO.write("%s\t%s\t%s\t%s\t%s\t%s\t%d\t%s\t%d\t%d\t%s\t%d\t%d\t%d\n"
                         % (smp, hap, lab, md, units[u][2], pol_used, nun, a.query_name[:16], rv,
                            u, ufd, c, m, ml))
    print("%s 끝" % fd.split("/")[-1], flush=True)
FO.close()
FS.close()
print("saved P10_CPG_%s.tsv.gz" % tag)
