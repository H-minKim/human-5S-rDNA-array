# -*- coding: utf-8 -*-
# 금화조 5S array 의 PacBio CpG 메틸 (제출본 v0.3 bigwig · 좌표 = GCA_048771995.1 / GCA_048772025.1)
#   PacBio 집계값(자리마다 %)이다 — ONT read 단위가 아니다. read 단위는 pod5 재콜링(GPU) 뒤에 가능
#   유닛 경계 = cmsearch 5S hit 시작점 (이웃 hit 시작 사이 = 유닛 하나) — bp/주기 환산 금지
#   산출 : zfinch_5S_units_meth.tsv (유닛마다 CpG 수 · 평균 %) · zfinch_5S_array_meth_summary.tsv (array 안 vs 양옆 5 kb)
#   인자 : 작업 폴더(~/xsp) · 최소 copy 수(기본 10)
import io
import os
import sys
import subprocess
import collections

NL = chr(10)
TAB = chr(9)
X = sys.argv[1]
MINC = int(sys.argv[2]) if len(sys.argv) > 2 else 10
Z = os.path.join(X, "zfinch")
BW2BED = os.path.join(X, "tools", "bigWigToBedGraph")
FLANK = 5000
ASM = {"mat": ("bTaeGut7_mat_GCA_048771995.1", "sub_mat.bw", "zfinch_v03_mat"),
       "pat": ("bTaeGut7_pat_GCA_048772025.1", "sub_pat.bw", "zfinch_v03_pat")}


def gb2bw(h):
    M = {}
    for L in io.open(os.path.join(Z, "len_match_%s.tsv" % h), encoding="utf8"):
        f = L.rstrip(NL).split(TAB)
        if len(f) == 3 and f[2] != "-" and "," not in f[2]:
            M[f[2]] = f[0]
    return M


def hits(h):
    H = collections.defaultdict(list)
    for L in io.open(os.path.join(X, "cms", ASM[h][2] + ".tbl"), encoding="utf8"):
        if L.startswith("#"):
            continue
        f = L.split()
        if len(f) >= 16:
            a, b = int(f[7]), int(f[8])
            H[f[0]].append((min(a, b), max(a, b), f[9]))
    for k in H:
        H[k].sort()
    return H


def bedgraph(bw, chrom, s, e):
    r = subprocess.run([BW2BED, "-chrom=" + chrom, "-start=%d" % max(0, s), "-end=%d" % e, bw, "/dev/stdout"],
                       capture_output=True, text=True)
    out = []
    for L in r.stdout.splitlines():
        f = L.split()
        if len(f) == 4:
            out.append((int(f[1]), int(f[2]), float(f[3])))
    return out


U = io.open(os.path.join(Z, "zfinch_5S_units_meth.tsv"), "w", encoding="utf8", newline=NL)
U.write(TAB.join(["hap", "assembly", "chrom_genbank", "chrom_bigwig", "array_start", "array_end", "strand", "unit_idx",
                  "unit_start", "unit_end", "n_cpg", "mean_meth_pct"]) + NL)
S = io.open(os.path.join(Z, "zfinch_5S_array_meth_summary.tsv"), "w", encoding="utf8", newline=NL)
S.write(TAB.join(["hap", "assembly", "chrom_genbank", "chrom_bigwig", "array_start", "array_end", "n_copies", "status",
                  "cpg_in_array", "mean_in_array", "frac_units_ge50", "cpg_left5kb", "mean_left5kb", "cpg_right5kb", "mean_right5kb"]) + NL)
cat = [L.rstrip(NL).split(TAB) for L in io.open(os.path.join(Z, "catalog_v03.tsv"), encoding="utf8")]
hd = dict((n, i) for i, n in enumerate(cat[0]))
for h in ("mat", "pat"):
    asm, bw, _ = ASM[h]
    M = gb2bw(h)
    H = hits(h)
    for r in cat[1:]:
        if r[hd["assembly"]] != asm or int(r[hd["n_copies"]]) < MINC:
            continue
        chrom, s, e = r[hd["chrom"]], int(r[hd["start"]]), int(r[hd["end"]])
        bwc = M.get(chrom)
        if not bwc:
            S.write(TAB.join([h, asm, chrom, "-", str(s), str(e), r[hd["n_copies"]], r[hd["status"]]] + ["NA"] * 7) + NL)
            continue
        bg = bedgraph(os.path.join(Z, bw), bwc, s - FLANK - 1, e + FLANK)
        st = [x[0] for x in H[chrom] if s <= x[0] <= e]
        units = [(st[i], st[i + 1] - 1) for i in range(len(st) - 1)]
        um = []
        j = 0
        for i, (a, b) in enumerate(units):
            v = [x[2] for x in bg if a - 1 <= x[0] < b]
            m = sum(v) / len(v) if v else float("nan")
            um.append(m)
            U.write(TAB.join([h, asm, chrom, bwc, str(s), str(e), r[hd["strand"]], str(i + 1), str(a), str(b), str(len(v)),
                              "%.2f" % m if v else "NA"]) + NL)
        ins = [x[2] for x in bg if s - 1 <= x[0] < e]
        lf = [x[2] for x in bg if x[0] < s - 1]
        rt = [x[2] for x in bg if x[0] >= e]
        good = [m for m in um if m == m]
        S.write(TAB.join([h, asm, chrom, bwc, str(s), str(e), r[hd["n_copies"]], r[hd["status"]], str(len(ins)),
                          "%.2f" % (sum(ins) / len(ins)) if ins else "NA",
                          "%.3f" % (sum(1 for m in good if m >= 50) / float(len(good))) if good else "NA",
                          str(len(lf)), "%.2f" % (sum(lf) / len(lf)) if lf else "NA",
                          str(len(rt)), "%.2f" % (sum(rt) / len(rt)) if rt else "NA"]) + NL)
U.close()
S.close()
for L in io.open(os.path.join(Z, "zfinch_5S_array_meth_summary.tsv"), encoding="utf8"):
    print(L.rstrip(NL))
