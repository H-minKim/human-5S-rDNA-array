# -*- coding: utf-8 -*-
# 각 array 에서 직접 반복을 확인한다 (갭 허용 정렬)
#   대표 유닛을 자로 array 조각에 정렬해 : 몇 번 나오는가 · array 의 몇 %를 덮는가 ·
#   연속 출현 시작점 사이 간격(= 주기, 인델 있어도 흔들리지 않음) · 유닛끼리 일치도
#   공개 반복단위가 있는 계통은 그것으로도 같이 재서 내 유닛과 대조한다
#   산출 : TILE_CHECK.tsv
# MIT License · Copyright (c) 2026 Hungmin Kim
import io
import os
import sys
import glob
import subprocess
import statistics

NL = chr(10)
TAB = chr(9)
H = os.path.expanduser("~/lit5s")
RULER = {"human": "X12811.1", "chimp": "X12811.1", "bonobo": "X12811.1", "gorilla": "X12811.1",
         "orang_sumatran": "X12811.1", "orang_bornean": "X12811.1", "siamang": "X12811.1",
         "rhesus": "X12811.1", "arabidopsis": "AY130738.1", "zebrafish": "AF213516.1",
         "mouse": "D43630.1", "zfinch": "AF419701.3"}


def rdfa(p):
    d, k = {}, None
    for L in io.open(p, encoding="utf8"):
        L = L.strip()
        if not L:
            continue
        if L[0] == ">":
            k = L[1:].split()[0]
            d[k] = []
        else:
            d[k].append(L.upper())
    return dict((a, "".join(b)) for a, b in d.items())


ALL = rdfa(os.path.join(H, "ALL_lit.fa"))


def hsps(qfa, sfa):
    r = subprocess.run(["blastn", "-task", "blastn", "-word_size", "14", "-dust", "no", "-evalue", "1e-10",
                        "-max_target_seqs", "5", "-query", qfa, "-subject", sfa,
                        "-outfmt", "6 pident length qstart qend sstart send bitscore"],
                       capture_output=True, text=True)
    out = []
    for L in r.stdout.splitlines():
        f = L.split(TAB)
        out.append((float(f[0]), int(f[1]), int(f[2]), int(f[3]), int(f[4]), int(f[5]), float(f[6])))
    return out


def tile(qfa, sfa, slen):
    hs = hsps(qfa, sfa)
    hs.sort(key=lambda x: -x[-1])
    taken, pid = [], []
    for p, aln, qs, qe, ss, se, bs in hs:
        a, b = min(ss, se), max(ss, se)
        if any(not (b < x or a > y) for x, y in taken):
            continue
        taken.append((a, b))
        pid.append(p)
    if not taken:
        return 0, 0.0, [], 0.0
    cov = sum(b - a + 1 for a, b in taken) / float(slen)
    st = sorted(x for x, y in taken)
    gaps = [st[i + 1] - st[i] for i in range(len(st) - 1)]
    return len(taken), cov, gaps, statistics.median(sorted(pid))


w = io.open(os.path.join(H, "TILE_CHECK.tsv"), "w", encoding="utf8", newline=NL)
w.write(TAB.join(["array_id", "span", "slice_len", "my_unit_len",
                  "U_copies", "U_cov", "U_spacing_med", "U_spacing_p10", "U_spacing_p90", "U_pident_med",
                  "ruler", "R_copies", "R_cov", "R_spacing_med", "R_pident_med", "verdict"]) + NL)
print("%-50s %6s %5s %6s %5s %8s %6s %6s  %s" % ("array_id", "unit", "copy", "덮힘", "주기", "p10-p90", "일치", "R주기", "판정"))
for sfa in sorted(glob.glob(os.path.join(H, "slices", "*.arr.fa"))):
    aid = os.path.basename(sfa).replace(".arr.fa", "")
    hd = io.open(sfa, encoding="utf8").readline()
    kv = dict(x.split("=", 1) for x in hd.split() if "=" in x)
    span, slen = int(kv["span"]), int(kv["slice_len"])
    ufa = os.path.join(H, "units", aid + ".unit.fa")
    ul = len(list(rdfa(ufa).values())[0])
    n, cov, gaps, pid = tile(ufa, sfa, slen)
    gs = sorted(gaps)
    gm = int(statistics.median(gs)) if gs else 0
    g10 = gs[int(0.1 * (len(gs) - 1))] if gs else 0
    g90 = gs[int(0.9 * (len(gs) - 1))] if gs else 0
    sp = aid.split("__")[0]
    key = [k for k in RULER if sp.startswith(k)]
    rn, rc, rg, rp = "NA", 0, 0.0, 0
    rpid = 0.0
    if key and RULER[key[0]] in ALL:
        rn = RULER[key[0]]
        qq = os.path.join(H, "rq.fa")
        io.open(qq, "w", encoding="utf8", newline=NL).write(">" + rn + NL + ALL[rn] + NL)
        rc, rcov, rgaps, rpid = tile(qq, sfa, slen)
        rgs = sorted(rgaps)
        rg = int(statistics.median(rgs)) if rgs else 0
    else:
        rcov = 0.0
    if n == 0:
        v = "유닛이 array 에 안 붙음"
    elif n < 3:
        v = "출현 %d 회 — 반복 판정 불가" % n
    elif abs(gm - ul) <= max(0.10 * ul, 20) and cov >= 0.80:
        v = "OK 직접 확인 (연속 타일링)"
    elif abs(gm - ul) <= max(0.10 * ul, 20):
        v = "주기 일치 · 덮힘 %.2f" % cov
    elif gm < 0.75 * ul:
        v = "내 유닛이 주기(%d)보다 길다 — 덩어리" % gm
    else:
        v = "주기 %d != 유닛 %d" % (gm, ul)
    print("%-50s %6d %5d %6.2f %6d %4d-%-4d %6.1f %6d  %s" % (aid, ul, n, cov, gm, g10, g90, pid, rg, v))
    w.write(TAB.join([aid, str(span), str(slen), str(ul), str(n), "%.3f" % cov, str(gm), str(g10), str(g90),
                      "%.1f" % pid, rn, str(rc), "%.3f" % rcov, str(rg), "%.1f" % rpid, v]) + NL)
    sys.stdout.flush()
w.close()
print("표 : ~/lit5s/TILE_CHECK.tsv")
