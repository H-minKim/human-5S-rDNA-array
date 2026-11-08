# -*- coding: utf-8 -*-
# 5S array list, built from cmsearch --tblout
#   hit 사이 간격 10 kb 이하 = 같은 array · 기록 : 염색체 · 시작 · 끝 · copy 수 · 단위 주기(이웃 hit 시작 간격 중앙값) · 우세 가닥
#   ±100 kb 안의 N 연속(gap) · array 안의 N · 판정
#   판정 : array 안과 ±100 kb 에 N 없음 = 조립완료(COMPLETE) / 있으면 GAP_NEAR 또는 GAP_INSIDE
#   copy 수는 hit 수로 센다 (bp/주기 환산 금지)
#   인자 : 종 · 어셈블리 이름 · tblout · fasta(.fai 필요) · 출력 tsv(덧붙임)
import io
import os
import re
import sys
import collections

NL = chr(10)
TAB = chr(9)
sp, asm, tbl, fa, out = sys.argv[1:6]
MINCOPY = int(os.environ.get("MINCOPY", "3"))
GAP = 10000
PAD = 100000
H = collections.defaultdict(list)
for L in io.open(tbl, encoding="utf8"):
    if L.startswith("#"):
        continue
    f = L.split()
    if len(f) < 16:
        continue
    chrom, a, b, strand, score = f[0], int(f[7]), int(f[8]), f[9], float(f[14])
    H[chrom].append((min(a, b), max(a, b), strand, score))
FAI = {}
for L in io.open(fa + ".fai", encoding="utf8"):
    f = L.split(TAB)
    FAI[f[0]] = (int(f[1]), int(f[2]), int(f[3]), int(f[4]))


def fetch(chrom, s, e):
    # 1 기반 포함 좌표를 .fai 로 직접 읽는다 (samtools 없이)
    ln, off, lb, lw = FAI[chrom]
    s = max(1, s)
    e = min(ln, e)
    with open(fa, "rb") as fh:
        def pos(p):
            return off + (p - 1) // lb * lw + (p - 1) % lb
        fh.seek(pos(s))
        raw = fh.read(pos(e) - pos(s) + 1).decode("ascii")
    return raw.replace(NL, "").replace(chr(13), "")


NRUN = re.compile("[Nn]+")
rows = []
for chrom, hs in H.items():
    hs.sort()
    grp = [[hs[0]]]
    for h in hs[1:]:
        if h[0] - grp[-1][-1][1] <= GAP:
            grp[-1].append(h)
        else:
            grp.append([h])
    for g in grp:
        if len(g) < MINCOPY:
            continue
        s, e = g[0][0], max(x[1] for x in g)
        st = [x[0] for x in g]
        per = sorted(st[i + 1] - st[i] for i in range(len(st) - 1))
        med = per[len(per) // 2] if per else 0
        c = collections.Counter(x[2] for x in g)
        strand = c.most_common(1)[0][0]
        ln = FAI[chrom][0]
        inside = sum(len(m.group()) for m in NRUN.finditer(fetch(chrom, s, e)))
        left = fetch(chrom, s - PAD, s - 1) if s > 1 else ""
        right = fetch(chrom, e + 1, e + PAD) if e < ln else ""
        near = sum(len(m.group()) for m in NRUN.finditer(left)) + sum(len(m.group()) for m in NRUN.finditer(right))
        status = "GAP_INSIDE" if inside else ("GAP_NEAR" if near else "COMPLETE")
        edge = "CHROM_EDGE" if (s <= PAD or e >= ln - PAD) else "-"
        rows.append([sp, asm, chrom, str(s), str(e), str(e - s + 1), str(len(g)), str(med),
                     "%s(%d/%d)" % (strand, c[strand], len(g)), str(inside), str(near), edge, status,
                     "%.1f" % (sum(x[3] for x in g) / len(g))])
rows.sort(key=lambda r: -int(r[6]))
new = not os.path.exists(out)
o = io.open(out, "a", encoding="utf8", newline=NL)
if new:
    o.write(TAB.join(["species", "assembly", "chrom", "start", "end", "span_bp", "n_copies", "period_median",
                      "strand", "N_inside", "N_within_100kb", "chrom_edge", "status", "mean_bitscore"]) + NL)
for r in rows:
    o.write(TAB.join(r) + NL)
o.close()
tot = sum(len(v) for v in H.values())
print("%s %s : cmsearch hit %d · array(≥%d copy) %d · array 안 copy %d · 흩어진 hit %d"
      % (sp, asm, tot, MINCOPY, len(rows), sum(int(r[6]) for r in rows), tot - sum(int(r[6]) for r in rows)))
for r in rows[:8]:
    print("  %s:%s-%s  %s bp · copy %s · 주기 %s · %s · N안 %s · ±100kb N %s · %s %s" % (r[2], r[3], r[4], r[5], r[6], r[7], r[8], r[9], r[10], r[12], r[11]))
