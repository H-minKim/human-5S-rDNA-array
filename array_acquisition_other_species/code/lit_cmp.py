# -*- coding: utf-8 -*-
# 논문·GenBank 5S 반복단위 vs 내가 실측한 대표 유닛 대조
#   입력 : ~/lit5s/ALL_lit.fa (공개 서열) · ~/lit5s/units/<array_id>.unit.fa (실측 대표 유닛 77)
#   1) 공개 서열 -> 실측 유닛 blastn : 어느 array 유닛과 가장 닮았나 (일치도 · 덮힘)
#   2) 실측 유닛 안에 공개 반복단위가 몇 번 들어가는가 (겹치지 않는 HSP 수)
#      -> 2 이상이면 내 "유닛" 은 여러 유닛이 한 덩어리로 잡힌 것
#   산출 : LIT_BEST.tsv · LIT_TILING.tsv
# MIT License · Copyright (c) 2026 Hungmin Kim
import io
import os
import glob
import subprocess
import collections

NL = chr(10)
TAB = chr(9)
H = os.path.expanduser("~/lit5s")


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


LIT = rdfa(os.path.join(H, "ALL_lit.fa"))
LITTITLE = {}
for L in io.open(os.path.join(H, "ALL_lit.fa"), encoding="utf8"):
    if L.startswith(">"):
        p = L[1:].strip().split(None, 1)
        LITTITLE[p[0]] = (p[1] if len(p) > 1 else "")[:70]
UNIT = {}
for p in sorted(glob.glob(os.path.join(H, "units", "*.unit.fa"))):
    aid = os.path.basename(p).replace(".unit.fa", "")
    UNIT[aid] = list(rdfa(p).values())[0]
with io.open(os.path.join(H, "units_all.fa"), "w", encoding="utf8", newline=NL) as g:
    for a in sorted(UNIT):
        g.write(">" + a + NL + UNIT[a] + NL)


def blast(q, s, extra=""):
    r = subprocess.run("blastn -task blastn -word_size 11 -dust no -evalue 1e-5 -query %s -subject %s %s "
                       "-outfmt '6 qseqid sseqid pident length qstart qend sstart send bitscore'"
                       % (q, s, extra), shell=True, capture_output=True, text=True)
    out = []
    for L in r.stdout.splitlines():
        f = L.split(TAB)
        out.append((f[0], f[1], float(f[2]), int(f[3]), int(f[4]), int(f[5]), int(f[6]), int(f[7]), float(f[8])))
    return out


rows = blast(os.path.join(H, "ALL_lit.fa"), os.path.join(H, "units_all.fa"))
best = {}
for q, s, pid, aln, qs, qe, ss, se, bs in rows:
    k = (q, s)
    if k not in best or bs > best[k][-1]:
        best[k] = (pid, aln, qs, qe, ss, se, bs)
w = io.open(os.path.join(H, "LIT_BEST.tsv"), "w", encoding="utf8", newline=NL)
w.write(TAB.join(["lit_acc", "lit_len", "lit_title", "best_array_id", "unit_len", "pident",
                  "aln_len", "cov_of_lit", "cov_of_unit", "n_arrays_hit"]) + NL)
print("=== 1. 공개 반복단위가 가장 닮은 실측 유닛")
for q in sorted(LIT, key=lambda x: -len(LIT[x])):
    hits = [(s, v) for (qq, s), v in best.items() if qq == q]
    if not hits:
        print("  %-10s %5d bp  실측 유닛 어디에도 안 붙음" % (q, len(LIT[q])))
        w.write(TAB.join([q, str(len(LIT[q])), LITTITLE.get(q, ""), "NOHIT", "", "", "", "", "", "0"]) + NL)
        continue
    hits.sort(key=lambda x: -x[1][-1])
    s, v = hits[0]
    pid, aln, qs, qe, ss, se, bs = v
    cq = aln / float(len(LIT[q]))
    cu = aln / float(len(UNIT[s]))
    print("  %-10s %5d bp -> %-52s 유닛 %5d bp · 일치 %.1f%% · 정렬 %4d (공개의 %.2f · 유닛의 %.2f) · 붙은 array %d"
          % (q, len(LIT[q]), s, len(UNIT[s]), pid, aln, cq, cu, len(hits)))
    w.write(TAB.join([q, str(len(LIT[q])), LITTITLE.get(q, ""), s, str(len(UNIT[s])), "%.2f" % pid,
                      str(aln), "%.3f" % cq, "%.3f" % cu, str(len(hits))]) + NL)
w.close()

# 2) 타일링 : 실측 유닛 안에 공개 단위가 몇 번 들어가나
GRP = {"human": "human_CHM13v2", "macaca_fas": "rhesus_T2T_MMU8v2", "zebrafish": "zebrafish_GRCz12tu",
       "arabidopsis": "arabidopsis_ColCC", "chicken": "zfinch_v03_dip", "rat": "mouse_C57BL6J_T2T"}
ACC2GRP = {}
for g in GRP:
    for k in rdfa(os.path.join(H, g + ".lit.fa")) if os.path.exists(os.path.join(H, g + ".lit.fa")) else []:
        ACC2GRP[k] = g
if not ACC2GRP:  # 묶음 파일이 없으면 ALL 에서 유추 불가 -> 전부 검사
    ACC2GRP = dict((k, None) for k in LIT)
t = io.open(os.path.join(H, "LIT_TILING.tsv"), "w", encoding="utf8", newline=NL)
t.write(TAB.join(["array_id", "unit_len", "lit_acc", "lit_len", "n_nonoverlap_hsp", "covered_bp",
                  "cov_frac", "pident_med", "hsp_spans"]) + NL)
print(NL + "=== 2. 실측 유닛 안에 공개 반복단위가 몇 번 들어가나 (2 이상 = 여러 유닛 덩어리)")
for aid in sorted(UNIT, key=lambda a: -len(UNIT[a])):
    sp = aid.split("__")[0]
    accs = [a for a, g in ACC2GRP.items() if g and GRP.get(g) == sp]
    if not accs:
        continue
    ufa = os.path.join(H, "one.fa")
    io.open(ufa, "w", encoding="utf8", newline=NL).write(">" + aid + NL + UNIT[aid] + NL)
    for acc in accs:
        qfa = os.path.join(H, "q.fa")
        io.open(qfa, "w", encoding="utf8", newline=NL).write(">" + acc + NL + LIT[acc] + NL)
        hs = blast(qfa, ufa)
        hs.sort(key=lambda x: -x[-1])
        taken, spans, pids = [], [], []
        for q, s, pid, aln, qs, qe, ss, se, bs in hs:
            a, b = min(ss, se), max(ss, se)
            if any(not (b < x or a > y) for x, y in taken):
                continue
            taken.append((a, b))
            spans.append("%d-%d" % (a, b))
            pids.append(pid)
        if not taken:
            continue
        cov = sum(b - a + 1 for a, b in taken)
        pids.sort()
        t.write(TAB.join([aid, str(len(UNIT[aid])), acc, str(len(LIT[acc])), str(len(taken)), str(cov),
                          "%.3f" % (cov / float(len(UNIT[aid]))), "%.1f" % pids[len(pids) // 2],
                          ",".join(sorted(spans, key=lambda z: int(z.split("-")[0])))]) + NL)
t.close()
print("표 : ~/lit5s/LIT_BEST.tsv · ~/lit5s/LIT_TILING.tsv")
