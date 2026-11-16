# -*- coding: utf-8 -*-
# Per-unit m6A and CpG of the placed Fiber-seq reads (v12)
#   The skeleton is mkbam_ao_v11.py, for building the BAM and the reference, plus meth_units_ao_v11.py, for the coordinate rules, carried across unchanged.
#   Settings, through environment variables:
#     D     the working folder
#     CFG   AOJV_FS or AOJV_PB: only the SAMs of this cfg are read
#     LABEL UCSC_ONT_FIBERSEQ or PACBIO_FIBERSEQ, in place of the per-host srcmap table
#     OUT   the output folder
#     FINAL the final unit table, for the in_final column; column 8 is array_pos in G8
#     SMPS  optionally a comma-separated list of individuals, for a pilot run
#     SAMD  optionally the SAM folder, defaulting to D/sam
#   Where an MM block carries no terminating flag it is read as equivalent to '.', per the SAM specification (the fibertools output is A+a,/T-a,; jasmine and ONT state '.' explicitly and are unaffected).
#   CpG: the v11 code unchanged (C+m, minus 1 on the reverse strand, ML >= 128) -> METH_READUNIT and METH_SITE, with the hit rate against the canonical C and G.
#   m6A:
#     · A+a, the A of the original read, and T-a, the T of the original read, that is an A on the opposite strand, which occurs only in PacBio data.
#     · The flag '.' means a base absent from the list is 0, so the denominator is every base of that kind; the '?' of C+m means the list only.
#     · Being a single-base modification there is no minus-1 correction as there is for CpG.
#     · mod_strand is the strand the m6A lies on in the canonical frame: '+' is a canonical A position and '-' a canonical T position.
#         A+a is '+' on a forward read and '-' on a reverse read; T-a is '-' on a forward read and '+' on a reverse one.
#     · The hit rate is the fraction of positions where the canonical base equals the expected one, A or T; one position out and it falls to about 0.3.
#   Both ARRAY_ONLY and junction reads were realigned to the canonical window by ao_g8cut, so the offset is 0.
import io
import os
import re
import sys
import gzip
import glob
import bisect
import subprocess
import collections

NL = chr(10)
TAB = chr(9)
H = os.path.expanduser
D = H(os.environ["D"])
CFG = os.environ["CFG"]
LABEL = os.environ["LABEL"]
OUT = H(os.environ["OUT"])
FINAL = os.environ.get("FINAL", "")
SMPS = set(x for x in os.environ.get("SMPS", "").split(",") if x)
SAMD = H(os.environ.get("SAMD", os.path.join(D, "sam")))   # the SAM folder, defaulting to D/sam; the fibertools recall version uses pbft_v12/sam
MAPF = H("~/MAP_V11.tsv")
HOST = os.uname()[1]
CIG = re.compile(r"(\d+)([MIDNSHP=X])")
CMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
MLCUT = 128
OFF = 0
os.makedirs(os.path.join(OUT, "bam", "ref"), exist_ok=True)


def rc(s):
    return "".join(CMP.get(c, "N") for c in reversed(s))


def rd(p):
    return "".join(x.strip() for x in io.open(p, encoding="utf8").read().split(NL)[1:])


# -- the canonical sequence, built the same way as in mkbam_ao_v11
SEQF = {}
for L in io.open(os.path.join(D, "g8in", "hap_seq_node_v7.tsv"), encoding="utf8"):
    f = L.rstrip(NL).split(TAB)
    if len(f) >= 3:
        SEQF[(f[0], f[1])] = f[2]
FLIP = {}
f = io.open(os.path.join(D, "SUM_G8.tsv"), encoding="utf8")
f.readline()
for L in f:
    r = L.rstrip(NL).split(TAB)
    FLIP[(r[0], r[1])] = (r[4] == "Y")

# -- the unit boundaries, as in meth_units_ao_v11
U = collections.defaultdict(list)
f = io.open(MAPF, encoding="utf8")
ix = dict((n, i) for i, n in enumerate(f.readline().rstrip(NL).split(TAB)))
for L in f:
    r = L.rstrip(NL).split(TAB)
    if r[ix["seg_type"]] != "UNIT":
        continue
    U[(r[0], r[1])].append((int(r[ix["unit_start"]]), int(r[ix["unit_end"]]), int(r[ix["array_pos"]])))
for k in U:
    U[k].sort()
STARTS = dict((k, [x[0] for x in v]) for k, v in U.items())

# -- the OK rows of the final unit table, as (sample, read, slot, array_pos)
FIN = set()
if FINAL:
    for fp in glob.glob(H(FINAL)):
        for L in io.open(fp, encoding="utf8"):
            r = L.rstrip(NL).split(TAB)
            if len(r) > 13 and r[4] == CFG and r[13] == "OK":
                FIN.add((r[0], r[1], r[3], r[7]))                  # r[7] is array_pos in G8 and r[8] the G6 index (see the comments in fix_relabel_only)


# -- 1. the BAM and the reference, exactly as in mkbam_ao_v11, with only the output folder new
def mkbam(sp, smp):
    recs = [L.rstrip(NL) for L in io.open(sp, encoding="utf8") if L.strip() and not L.startswith("@")]
    nzr = sum(1 for x in recs if (TAB + "ZR:Z:") in x)
    recs = [TAB.join(g for g in x.split(TAB) if not g.startswith("ZR:Z:")) for x in recs]
    used = sorted(set(x.split(TAB, 3)[2] for x in recs))
    hdr = ["@HD" + TAB + "VN:1.6" + TAB + "SO:unsorted"]
    refs = {}
    fa = io.open(os.path.join(OUT, "bam", "ref", "%s.%s.G8v3ref.fa" % (smp, CFG)), "w", encoding="utf8", newline=NL)
    for nm in used:
        s_, sl = nm.rsplit("_h", 1)
        if (s_, sl) not in SEQF:
            return None, "NO_REFSEQ"
        a = rd(os.path.join(D, "cmp4d_v7", SEQF[(s_, sl)])).upper()
        a = rc(a) if FLIP.get((s_, sl)) else a
        refs[nm] = a
        hdr.append("@SQ" + TAB + "SN:" + nm + TAB + "LN:%d" % len(a))
        fa.write(">" + nm + NL + a + NL)
    fa.close()
    hdr.append("@CO" + TAB + "Fiber-seq placed reads (AOJV) re-aligned to haplotype canonical (G8v3 orientation) by ao_g8cut_v12; SEQ/QUAL/MM/ML untouched")
    tmp = os.path.join(OUT, "bam", "tmp_%s.%s.sam" % (smp, CFG))
    io.open(tmp, "w", encoding="utf8", newline=NL).write(NL.join(hdr + recs) + NL)
    out = os.path.join(OUT, "bam", "%s.%s.%s.G8v3.bam" % (smp, CFG, LABEL))
    r1 = subprocess.run(["samtools", "sort", "-@2", "-o", out, tmp], capture_output=True, text=True)
    r2 = subprocess.run(["samtools", "index", out], capture_output=True, text=True) if r1.returncode == 0 else r1
    os.remove(tmp)
    if r1.returncode or r2.returncode:
        return None, "SORT_FAIL"
    return (out, refs), "ZR%d rec%d" % (nzr, len(recs))


def blocks(seq, flag, mm, ml):
    """(block name, [(SEQ coordinate, value)]); MM is relative to the original read orientation. A '.' block fills bases outside the list with 0."""
    orig = rc(seq) if (flag & 16) else seq
    vals = [int(x) for x in ml.split(",") if x != ""]
    res = []
    base = 0
    for blk in [p for p in mm.split(";") if p]:
        f2 = blk.split(",")
        deltas = [int(x) for x in f2[1:]]
        head = f2[0]
        nm = head.rstrip(".?")
        b0 = nm[0]
        if nm in ("C+m", "A+a", "T-a"):
            idx = [i for i, c in enumerate(orig) if c == b0]
            got = {}
            j = -1
            for t, d in enumerate(deltas):
                j += d + 1
                if j >= len(idx):
                    break
                got[j] = vals[base + t] if base + t < len(vals) else 0
            if nm == "C+m":
                items = sorted(got.items())                       # C+m: the list only
            elif head.endswith(".") or head[-1] not in ".?":
                items = [(j, got.get(j, 0)) for j in range(len(idx))]   # the denominator is every base of that kind
            else:
                items = sorted(got.items())
            out = []
            for j, v in items:
                op = idx[j]
                sp = (len(seq) - 1 - op) if (flag & 16) else op
                out.append((sp, v))
            res.append((nm, head[-1] if head[-1] in ".?" else "", out))
        base += len(deltas)
    return res


RU = gzip.open(os.path.join(OUT, "METH_READUNIT_%s.tsv.gz" % HOST), "wt", encoding="utf8", newline=NL)
RU.write(TAB.join(["sample", "data_source", "cfg", "kind", "slot", "array_pos", "read_id",
                   "strand", "n_cpg", "n_meth", "beta", "in_final", "calls"]) + NL)
SI = gzip.open(os.path.join(OUT, "METH_SITE_%s.tsv.gz" % HOST), "wt", encoding="utf8", newline=NL)
SI.write(TAB.join(["sample", "data_source", "slot", "array_pos", "unit_pos", "n", "n_meth", "beta"]) + NL)
AU = gzip.open(os.path.join(OUT, "M6A_READUNIT_%s.tsv.gz" % HOST), "wt", encoding="utf8", newline=NL)
AU.write(TAB.join(["sample", "data_source", "cfg", "kind", "slot", "array_pos", "read_id", "strand",
                   "mod_strand", "n_a", "n_m6a", "frac", "in_final", "calls"]) + NL)
AS = gzip.open(os.path.join(OUT, "M6A_SITE_%s.tsv.gz" % HOST), "wt", encoding="utf8", newline=NL)
AS.write(TAB.join(["sample", "data_source", "slot", "array_pos", "unit_pos", "mod_strand", "n", "n_m6a", "frac"]) + NL)
LOG = io.open(os.path.join(OUT, "LOG_%s.tsv" % HOST), "w", encoding="utf8", newline=NL)
LOG.write(TAB.join(["host", "cfg", "sample", "bam", "read", "no_mm", "mm_empty", "no_m6a_block", "read_unit", "in_final",
                    "cpg_calls", "cpg_hit", "cpg_hit_noshift",
                    "m6a_calls_plus", "m6a_hit_plus", "m6a_calls_minus", "m6a_hit_minus", "m6a_hit_shift1"]) + NL)

TOT = collections.Counter()
for sp_ in sorted(glob.glob(os.path.join(SAMD, "*.%s.sam" % CFG))):
    smp = os.path.basename(sp_).split(".")[0]
    if SMPS and smp not in SMPS:
        continue
    got, msg = mkbam(sp_, smp)
    if got is None:
        LOG.write(TAB.join([HOST, CFG, smp, msg]) + NL)
        continue
    b, refs = got
    T = collections.Counter()
    SITE = collections.Counter()
    SITEM = collections.Counter()
    ASITE = collections.Counter()
    ASITEM = collections.Counter()
    pr = subprocess.Popen(["samtools", "view", b], stdout=subprocess.PIPE, text=True)
    for L in pr.stdout:
        r = L.rstrip(NL).split(TAB)
        if len(r) < 11:
            continue
        flag = int(r[1])
        if flag & 0x900 or r[2] == "*":
            continue
        sn, pos, cig, seq = r[2], int(r[3]), r[5], r[9]
        sl = sn.split("_h")[-1]
        key = (smp, sl)
        if key not in U:
            T["no_map"] += 1
            continue
        mm = ml = None
        for x in r[11:]:
            if x.startswith("MM:Z:"):
                mm = x[5:]
            elif x.startswith("ML:B:C,"):
                ml = x[7:]
        T["read"] += 1
        if not mm:
            T["no_mm"] += 1
            continue
        if not ml:
            T["mm_empty"] += 1                    # an MM block present but with no calls, the ML being empty, is dropped rather than filled with zeros
            continue
        if "A+a" not in mm and "T-a" not in mm:
            T["no_m6a_block"] += 1                # the raw data carries no m6A calls; the CpG part continues
        ref = refs[sn]
        m = {}
        rp = pos - 1
        q = 0
        for g in CIG.finditer(cig):
            k, op = int(g.group(1)), g.group(2)
            if op == "S":
                q += k
            elif op == "H":
                pass
            elif op in "M=X":
                for t in range(k):
                    m[q + t] = rp + t
                rp += k
                q += k
            elif op == "I":
                q += k
            elif op in "DN":
                rp += k
        st = "-" if (flag & 16) else "+"
        S = STARTS[key]
        per = collections.defaultdict(list)
        aper = collections.defaultdict(list)
        for nm, fl, calls in blocks(seq, flag, mm, ml):
            if nm == "C+m":
                for sp, v in calls:
                    if sp not in m:
                        continue
                    cp = m[sp] - OFF
                    # the test: the position before the correction
                    if 0 <= cp < len(ref) and ref[cp] == ("G" if flag & 16 else "C"):
                        T["cpg_hit_noshift"] += 1
                    if flag & 16:
                        cp -= 1
                    T["cpg_calls"] += 1
                    if 0 <= cp < len(ref) - 1 and ref[cp:cp + 2] == "CG":
                        T["cpg_hit"] += 1
                    p1 = cp + 1
                    if p1 < 1:
                        continue
                    i = bisect.bisect_right(S, p1) - 1
                    if i < 0:
                        continue
                    s0, e0, ap = U[key][i]
                    if p1 > e0:
                        continue
                    per[ap].append((p1 - s0 + 1, v))
            else:
                ms = "+" if ((nm == "A+a") != bool(flag & 16)) else "-"
                want = "A" if ms == "+" else "T"
                for sp, v in calls:
                    if sp not in m:
                        continue
                    cp = m[sp] - OFF
                    T["m6a_calls" + ms] += 1
                    if 0 <= cp < len(ref) and ref[cp] == want:
                        T["m6a_hit" + ms] += 1
                    if 0 <= cp + 1 < len(ref) and ref[cp + 1] == want:
                        T["m6a_hit_shift1"] += 1
                    p1 = cp + 1
                    if p1 < 1:
                        continue
                    i = bisect.bisect_right(S, p1) - 1
                    if i < 0:
                        continue
                    s0, e0, ap = U[key][i]
                    if p1 > e0:
                        continue
                    aper[(ap, ms)].append((p1 - s0 + 1, v))
        for ap in sorted(set(list(per.keys()) + [a for a, _ in aper.keys()])):
            T["read_unit"] += 1
            inf = "Y" if (smp, r[0], sl, str(ap)) in FIN else "N"
            if inf == "Y":
                T["in_final"] += 1
            lst = sorted(per.get(ap, []))
            if lst:
                nm_ = sum(1 for _, v in lst if v >= MLCUT)
                RU.write(TAB.join([smp, LABEL, CFG, "aojv", sl, str(ap), r[0], st, str(len(lst)), str(nm_),
                                   "%.4f" % (nm_ / float(len(lst))), inf,
                                   ",".join("%d:%d" % (a, v) for a, v in lst)]) + NL)
                for a, v in lst:
                    SITE[(sl, ap, a)] += 1
                    if v >= MLCUT:
                        SITEM[(sl, ap, a)] += 1
            for ms in "+-":
                al = sorted(aper.get((ap, ms), []))
                if not al:
                    continue
                na = sum(1 for _, v in al if v >= MLCUT)
                AU.write(TAB.join([smp, LABEL, CFG, "aojv", sl, str(ap), r[0], st, ms, str(len(al)), str(na),
                                   "%.4f" % (na / float(len(al))), inf,
                                   ",".join("%d:%d" % (a, v) for a, v in al)]) + NL)
                for a, v in al:
                    ASITE[(sl, ap, a, ms)] += 1
                    if v >= MLCUT:
                        ASITEM[(sl, ap, a, ms)] += 1
    pr.stdout.close()
    pr.wait()
    for k in sorted(SITE, key=lambda x: (int(x[0]), x[1], x[2])):
        n = SITE[k]
        nm_ = SITEM.get(k, 0)
        SI.write(TAB.join([smp, LABEL, k[0], str(k[1]), str(k[2]), str(n), str(nm_), "%.4f" % (nm_ / float(n))]) + NL)
    for k in sorted(ASITE, key=lambda x: (int(x[0]), x[1], x[2], x[3])):
        n = ASITE[k]
        na = ASITEM.get(k, 0)
        AS.write(TAB.join([smp, LABEL, k[0], str(k[1]), str(k[2]), k[3], str(n), str(na), "%.4f" % (na / float(n))]) + NL)

    def rt(a, b):
        return "%.4f" % (T[a] / float(T[b])) if T[b] else "NA"
    tot_m6a = T["m6a_calls+"] + T["m6a_calls-"]
    LOG.write(TAB.join([HOST, CFG, smp, os.path.basename(b), str(T["read"]), str(T["no_mm"]), str(T["mm_empty"]), str(T["no_m6a_block"]), str(T["read_unit"]),
                        str(T["in_final"]), str(T["cpg_calls"]), rt("cpg_hit", "cpg_calls"), rt("cpg_hit_noshift", "cpg_calls"),
                        str(T["m6a_calls+"]), rt("m6a_hit+", "m6a_calls+"), str(T["m6a_calls-"]), rt("m6a_hit-", "m6a_calls-"),
                        "%.4f" % (T["m6a_hit_shift1"] / float(tot_m6a)) if tot_m6a else "NA"]) + NL)
    LOG.flush()
    TOT.update(T)
    print("  %s %s read %d · read-unit %d · in_final %d" % (CFG, smp, T["read"], T["read_unit"], T["in_final"]), flush=True)

for fh in (RU, SI, AU, AS, LOG):
    fh.close()
print("DONE %s %s" % (HOST, CFG))
