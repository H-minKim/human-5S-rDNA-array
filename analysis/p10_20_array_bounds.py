# -*- coding: utf-8 -*-
#  For each P10 block, does the array reach non-repetitive flank at both ends?
#    The question: do the six arrays separate the flank correctly, and is the boundary between array and non-array really there?
#    For each P10 block, found in the assembly through a 10 kb window and cut with +-30 kb of slack, the unit intervals are
#    located again by the same landmark rule as p10_08, and the non-repetitive length remaining before the first unit and after the last, that is the flank, is measured.
#    The CHM13 flank 31-mer hits from P10_SEGMENT_POLARITY are attached as well, as head or tail against LEFT or RIGHT.
#    Read-only: 81 ref.fa folders and three tables. No BAM is opened.
#    Used by per_array_panels_fig04_G-L_figS4_G-J_figS6_A-U.py: the panels drawn along the array use only blocks
#      whose left_nonrepeat_bp_anyhit and right_nonrepeat_bp_anyhit are both at least 500 bp.
import collections
import glob
import hashlib
import io
import os
import statistics

def _need(_v):
    """The deposited copy has no default paths. Stop when a root variable is missing."""
    _x = os.environ.get(_v)
    if not _x:
        sys.exit("environment variable %s is not set; see README.md, section 'Environment variables'" % _v)
    return _x


P = os.path.join(_need("RDNA_WORK_ROOT"), "P10")
OUTD = os.path.join(_need("RDNA_WORK_ROOT"), "out", "P10")
OUT = os.path.join(_need("ARRAYMETH_HANDOFF_ROOT"), "02_output_tables", "P10", "P10_ARRAY_BOUNDS.tsv")
K = 31
#  the same values as p10_08_extract_frame.py
FAM_UNIT = {"5S": "5S", "D4Z4_4q": "D4Z4_4q", "D4Z4_10q": "D4Z4_10q", "DXZ4": "DXZ4",
            "RNU2": "RNU2", "SST1": "SST1_a"}
ULEN = {"5S": 2240, "D4Z4_4q": 3298, "D4Z4_10q": 3323, "DXZ4": 2972, "RNU2": 6171, "SST1": 2527}


def rc(s):
    return s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


UN, name = {}, None
for L in io.open(os.path.join(P, "data", "P10_UNITS_CHM13.fa"), encoding="utf8"):
    if L[0] == ">":
        name = L[1:].split("_unit2")[0]
        UN[name] = []
    else:
        UN[name].append(L.strip().upper())
UN = {k: "".join(v) for k, v in UN.items()}
KM = {}
for fam, un in FAM_UNIT.items():
    L = ULEN[fam]
    u = UN[un][:L]
    circ = u + u[:K - 1]
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


def landmarks(seq, fam):
    """The landmarks of p10_08 as they are, and also returns the landmark list and the first and last positions of any hit (unit fragments included)"""
    d, L = KM[fam], ULEN[fam]
    est = {"+": [], "-": []}
    anyhit = []
    for j in range(0, len(seq) - K + 1):
        h = d.get(seq[j:j + K])
        if h is None:
            continue
        anyhit.append(j)
        if h[0] == "+":
            est["+"].append(j - h[1])
        else:
            est["-"].append(j + h[1] + K - 1)
    units, nlm, skipped = [], 0, 0
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
        nlm += len(lm)
        for k in range(len(lm) - 1):
            g = lm[k + 1] - lm[k]
            if not (0.8 * L <= g <= 1.2 * L):
                skipped += 1
                continue
            if st == "+":
                units.append((lm[k], lm[k + 1], "+"))
            else:
                units.append((lm[k] + 1, lm[k + 1] + 1, "-"))
    units.sort()
    clean = []
    for u in units:
        if clean and u[0] < clean[-1][1]:
            clean.pop()
            continue
        clean.append(u)
    return clean, skipped, (min(anyhit), max(anyhit) + K) if anyhit else (None, None)


REG = {}
for i, L in enumerate(io.open(os.path.join(P, "out", "P10_SEGMENT_REGISTRY.tsv"), encoding="utf8")):
    if i == 0:
        continue
    f = L.rstrip("\n").split("\t")
    REG[(f[0], f[1], f[2])] = (f[15], f[16], f[3])        # final_label, status, span
POLT = {}
with io.open(os.path.join(P, "out", "P10_SEGMENT_POLARITY.tsv"), encoding="utf8") as fh:
    h = fh.readline().rstrip("\n").split("\t")
    for L in fh:
        r = dict(zip(h, L.rstrip("\n").split("\t")))
        POLT[(r["sample"], r["hap"], r["span"])] = r
USED = {}
with io.open(os.path.join(OUTD, "P10_FRAME_SEGS_ALL.tsv"), encoding="utf8") as fh:
    h = fh.readline().rstrip("\n").split("\t")
    for L in fh:
        r = dict(zip(h, L.rstrip("\n").split("\t")))
        USED[(r["sample"], r["hap"], r["seg_md5"])] = r

rows, seen = [], set()
for fd in sorted(glob.glob(os.path.join(P, "data", "place", "*"))):
    smp = os.path.basename(fd).split("_")[0]
    ref = os.path.join(fd, "ref.fa")
    if not os.path.exists(ref):
        continue
    SEG, nm, buf = {}, None, []
    for L in io.open(ref, encoding="utf8"):
        if L[0] == ">":
            if nm:
                SEG[nm] = "".join(buf)
            nm, buf = L[1:].strip(), []
        else:
            buf.append(L.strip())
    if nm:
        SEG[nm] = "".join(buf)
    for nm, s in SEG.items():
        hap = nm.split("_")[2]
        md = hashlib.md5(s.encode()).hexdigest()[:12]
        if (smp, hap, md) in seen:
            continue
        lab, status, span = REG.get((smp, hap, md), ("?", "?", "?"))
        if status != "ok" or lab not in ULEN:
            continue
        seen.add((smp, hap, md))
        units, skipped, (h0, h1) = landmarks(s.upper(), lab)
        if len(units) < 2:
            continue
        a, b = units[0][0], units[-1][1]
        po = POLT.get((smp, hap, span), {})
        us = USED.get((smp, hap, md), {})
        rows.append((smp, hap, lab, md, span, len(s), len(units), skipped, a, len(s) - b,
                     h0, len(s) - h1, po.get("head_LEFT", "NA"), po.get("head_RIGHT", "NA"),
                     po.get("tail_LEFT", "NA"), po.get("tail_RIGHT", "NA"), po.get("orientation", "NA"),
                     us.get("polarity_used", "NA"), us.get("n_units", "NA")))
with io.open(OUT, "w", encoding="utf8", newline="\n") as fo:
    fo.write("sample\thap\tfamily\tseg_md5\tspan\tseg_len\tn_units\tskipped_gaps\tleft_flank_bp\tright_flank_bp\t"
             "left_nonrepeat_bp_anyhit\tright_nonrepeat_bp_anyhit\thead_LEFT\thead_RIGHT\ttail_LEFT\ttail_RIGHT\t"
             "orientation\tpolarity_used_p10_08\tn_units_p10_08\n")
    for r in rows:
        fo.write("\t".join(str(x) for x in r) + "\n")
print("조각 %d · 저장 %s" % (len(rows), OUT))
