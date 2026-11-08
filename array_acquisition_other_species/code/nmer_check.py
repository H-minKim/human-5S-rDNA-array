# -*- coding: utf-8 -*-
# 유닛 하나에 5S 유전자가 몇 개 들어가나 (n-mer 판정) — cmsearch 유전자 끝 간격의 차례를 본다
#   단량체 : 간격이 한 부류 (붉은털원숭이 기각 array 처럼 한 값에 몰림)
#   n 량체 : 간격이 두세 부류로 규칙적으로 번갈아 난다
#            (제브라피시 chr4 56.53 Mb = 164/313 교대 -> 유닛 477 에 유전자 2 개
#             붉은털원숭이 chr1 20.4 Mb = 2,9xx/4,4xx 교대 -> 유닛 7,440 에 유전자 2 개)
#   불규칙 : 어느 k 로도 차례가 맞지 않음 (유닛 경계 없음)
#   판정은 직접 측정 주기(옵셋 사다리)와 대조해 확정한다
#   산출 : NMER.tsv
# MIT License · Copyright (c) 2026 Hungmin Kim
import io
import os
import glob
import statistics
import collections

NL = chr(10)
TAB = chr(9)
BASE = "$RDNA_OTHERSP_DIR/data_v2"
LOC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "lit")
TOL = 0.08     # 같은 부류로 볼 상대 차이
MINFRAC = 0.80  # 차례가 맞아야 하는 비율


def classes(gaps):
    cs = []
    for g in sorted(gaps):
        for c in cs:
            if abs(g - c[0]) <= TOL * max(g, c[0]):
                c.append(g)
                break
        else:
            cs.append([g])
    cs.sort(key=lambda c: -len(c))
    return cs


def label(gaps, cs):
    rep = [statistics.median(c) for c in cs]
    out = []
    for g in gaps:
        out.append(min(range(len(rep)), key=lambda i: abs(g - rep[i]) / float(max(g, rep[i]))))
    return out


ST = {}
f = io.open(os.path.join(BASE, "00_catalog", "5S_array_structure.tsv"), encoding="utf8")
hd = f.readline().rstrip(NL).split(TAB)
ix = dict((n, i) for i, n in enumerate(hd))
for L in f:
    r = L.rstrip(NL).split(TAB)
    ST["%s__%s_%s_%s" % (r[ix["species"]], r[ix["chrom"]], r[ix["array_start"]], r[ix["array_end"]])] = r
P = {}
f = io.open(os.path.join(BASE, "00_catalog", "REFUNIT_PARAMS.tsv"), encoding="utf8")
hd = f.readline().rstrip(NL).split(TAB)
ipr = dict((n, i) for i, n in enumerate(hd))
for L in f:
    r = L.rstrip(NL).split(TAB)
    P[r[1]] = r
D = {}
f = io.open(os.path.join(LOC, "PERIOD_FINAL.tsv"), encoding="utf8")
hd = f.readline().rstrip(NL).split(TAB)
idr = dict((n, i) for i, n in enumerate(hd))
for L in f:
    r = L.rstrip(NL).split(TAB)
    D[r[0]] = r
HIT = {}
for p in glob.glob(os.path.join(BASE, "cms_tbl", "*.tbl")):
    sp = os.path.basename(p).replace(".tbl", "")
    H = collections.defaultdict(list)
    for L in io.open(p, encoding="utf8"):
        if L.startswith("#"):
            continue
        ff = L.split()
        if len(ff) < 9:
            continue
        a, b = int(ff[7]), int(ff[8])
        H[ff[0]].append((min(a, b), max(a, b)))
    for c in H:
        H[c].sort()
    HIT[sp] = H
w = io.open(os.path.join(LOC, "NMER.tsv"), "w", encoding="utf8", newline=NL)
w.write(TAB.join(["array_id", "copies", "n_gap", "n_class", "class_reps", "genes_per_unit", "cycle_frac",
                  "unit_len_from_gaps", "g8_refunit_len", "period_direct", "repeat_frac", "verdict"]) + NL)
print("%-50s %5s %5s %5s %4s %6s %7s %7s %6s  %s" % ("array", "copy", "간격", "부류", "k", "차례맞", "유닛(간격)", "G8유닛", "직접", "판정"))
cnt = collections.Counter()
for aid in sorted(D):
    sp = aid.split("__")[0]
    if sp not in HIT or aid not in ST:
        continue
    r = ST[aid]
    c, s0, e1 = r[ix["chrom"]], int(r[ix["array_start"]]), int(r[ix["array_end"]])
    hs = [h for h in HIT[sp].get(c, []) if s0 <= h[0] and h[1] <= e1]
    ends = sorted(h[1] for h in hs)
    gaps = [ends[i + 1] - ends[i] for i in range(len(ends) - 1)]
    if len(gaps) < 4:
        continue
    cs = classes(gaps)
    lab = label(gaps, cs)
    reps = [int(statistics.median(x)) for x in cs]
    nk, cyc, ul = 0, 0.0, 0
    if len(cs[0]) >= MINFRAC * len(gaps):
        nk, cyc, ul = 1, len(cs[0]) / float(len(gaps)), int(statistics.median(cs[0]))
    else:
        for k in (2, 3, 4):
            ok = sum(1 for i in range(len(lab) - k) if lab[i] == lab[i + k])
            fr = ok / float(max(1, len(lab) - k))
            if fr >= MINFRAC:
                sums = [sum(gaps[i:i + k]) for i in range(0, len(gaps) - k + 1, k)]
                nk, cyc, ul = k, fr, int(statistics.median(sums))
                break
    g8 = P[aid][ipr["refunit_len"]] if aid in P else "NA"
    pdv = int(D[aid][idr["period_direct"]])
    frc = float(D[aid][idr["frac_of_slice"]])
    if nk == 0:
        v = "불규칙 (유닛 경계 없음)"
    elif nk == 1:
        v = "단량체" if (g8 != "NA" and abs(ul - int(g8)) <= max(0.10 * ul, 25)) else "단량체 · G8 유닛과 어긋남"
    else:
        v = "%d 량체 (유전자 %d 개/유닛)" % (nk, nk)
        if g8 != "NA" and abs(ul - int(g8)) <= max(0.10 * ul, 25):
            v += " · G8 유닛과 같음"
    if pdv and abs(pdv - ul) > max(0.12 * max(ul, 1), 40):
        v += " · 직접주기(%d)와 어긋남" % pdv
    cnt[v.split(" ·")[0].split(" (")[0]] += 1
    w.write(TAB.join([aid, r[ix["copies"]], str(len(gaps)), str(len(cs)),
                      ",".join(str(x) for x in reps[:4]), str(nk), "%.2f" % cyc, str(ul), g8,
                      str(pdv), "%.3f" % frc, v]) + NL)
    if nk != 1 or "어긋남" in v:
        print("%-50s %5s %5d %5d %4d %6.2f %7d %7s %6d  %s" % (aid.split("__")[1], r[ix["copies"]], len(gaps),
                                                               len(cs), nk, cyc, ul, g8, pdv, v))
w.close()
print()
for k, v in cnt.most_common():
    print("   %-22s %d" % (k, v))
print("표 : lit/NMER.tsv")
