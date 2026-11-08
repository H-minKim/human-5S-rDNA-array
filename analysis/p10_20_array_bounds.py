# -*- coding: utf-8 -*-
#  p10_20 — P10 조각마다 array 의 양 끝이 flank(비반복 서열)에 닿았나 (2026-10-30)
#    사용자 물음 「array 6 개는 flank 를 정확히 분리한거 맞아? … array 와 array 가 아닌것의 경계가 제대로 있어?」
#    P10 조각(조립본에서 10 kb 창으로 찾아 ±30 kb 여유로 자른 것)마다 p10_08 과 같은 랜드마크 규칙으로
#    유닛 구간을 다시 찾고, 첫 유닛 앞 · 마지막 유닛 뒤에 남은 비반복 길이(= flank)를 잰다.
#    + P10_SEGMENT_POLARITY 의 CHM13 flank 31-mer 적중(head/tail × LEFT/RIGHT)을 붙인다.
#    읽기만 한다 (ref.fa 81 폴더 · 표 셋). BAM 은 안 연다.
#    첫 판 = `04_result_records\2026-10-30_P10_array_boundary_check\array_bounds16_2026-10-30.py` (출력 경로만 다르다)
#    쓰는 곳 : per_array_panels_fig04_G-L_figS4_G-J_figS6_A-U.py — array 를 따라 그린 칸은
#      `left_nonrepeat_bp_anyhit` · `right_nonrepeat_bp_anyhit` 둘 다 500 bp 이상인 조각만 쓴다.
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
        sys.exit("environment variable %s is not set; see README.md section Figures after 30 October 2026" % _v)
    return _x


P = os.path.join(_need("RDNA_WORK_ROOT"), "P10")
OUTD = os.path.join(_need("RDNA_WORK_ROOT"), "out", "P10")
OUT = os.path.join(_need("ARRAYMETH_HANDOFF_ROOT"), "02_output_tables", "P10", "P10_ARRAY_BOUNDS.tsv")
K = 31
#  ↓ p10_08_extract_frame.py 와 같은 값
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
    """p10_08 의 landmarks 그대로 + 랜드마크 목록과 아무 적중이나 처음·끝 자리(유닛 조각 포함)도 돌려준다"""
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
