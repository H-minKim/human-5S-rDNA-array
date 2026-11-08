# -*- coding: utf-8 -*-
# 제브라 short-read 메틸 — chr18(체세포형 array 1) vs chr4(모계형 array 묶음) 을 시료(조직·시기)별로
#   read 단위 : WGBS_PERREAD_input.tsv (= 12_zebrafish_devseries 의 ADULT_PERREAD, adult_stats.py 입력과 같은 파일)
#   array 전체 = 코딩부 + 스페이서(같은 read 는 합침). flank 는 쓰지 않음
#   차이 = chr18 − chr4, 두 군 read 독립 부트스트랩(정렬 전 차이)
#   cite : chr18 쪽 read < 10 인용 금지 · 10–29 약함 · ≥30 가능 
#   Stage III 난모세포(GSE175951)는 저자 가공 CpG 표라 read 가 없음 → CpG 자리 수 합산값만 따로(구간 없음)
# MIT License · Copyright (c) 2026 Hungmin Kim
import collections
import csv
import io
import numpy as np

rng = np.random.default_rng(0)
NB = 2000
R = "$RDNA_OTHERSP_DIR/data_v2/"
D = R + "13_array_level_diff/"
T = "\t"
SRC = {"알_Jiang": ("알(난자)", "MBT 전", "PRJNA188516 Jiang2013"),
       "16세포": ("16세포", "MBT 전", "PRJNA188516 Jiang2013"),
       "1k세포": ("1k세포", "MBT 직전/중", "PRJNA188516 Jiang2013"),
       "Germring": ("Germring 6hpf", "MBT 후", "PRJNA188516 Jiang2013"),
       "정자_Jiang": ("정자", "생식세포", "PRJNA188516 Jiang2013"),
       "정소": ("정소", "성체 생식선", "PRJNA188516 Jiang2013"),
       "알_Potok": ("알(난자)", "MBT 전", "PRJNA193544 Potok2013"),
       "sphere": ("sphere ≈4hpf", "MBT 후", "PRJNA193544 Potok2013"),
       "정자_Potok": ("정자", "생식세포", "PRJNA193544 Potok2013"),
       "근육_Potok": ("근육", "성체 체세포", "PRJNA193544 Potok2013"),
       "배아24hpf": ("배아 24hpf(EM-seq)", "MBT 후", "PRJNA628650"),
       "장": ("장", "성체 체세포", "PRJNA628650"),
       "심장": ("심장", "성체 체세포", "PRJNA628650"),
       "early_Control_1to28": ("생식선 체세포 1–28dpf", "생식선", "GSE122695 Ortega-Recalde2019"),
       "early_GSC_1to28": ("생식줄기세포 1–28dpf", "생식세포", "GSE122695 Ortega-Recalde2019"),
       "F_Control_35to70": ("난소 체세포 35–70dpf", "생식선", "GSE122695 Ortega-Recalde2019"),
       "F_GSC_35to70": ("난소 생식세포 35–70dpf", "생식세포", "GSE122695 Ortega-Recalde2019"),
       "M_Control_35to70": ("정소 체세포 35–70dpf", "생식선", "GSE122695 Ortega-Recalde2019"),
       "M_GSC_35to70": ("정소 생식세포 35–70dpf", "생식세포", "GSE122695 Ortega-Recalde2019")}

d = collections.defaultdict(lambda: collections.defaultdict(lambda: {"코딩부": [0, 0], "스페이서": [0, 0]}))
for r in csv.DictReader(io.open(D + "WGBS_PERREAD_input.tsv", encoding="utf8"), delimiter=T):
    if r["part"] not in ("코딩부", "스페이서"):
        continue
    x = d[(r["group"], r["branch"])][r["read"]][r["part"]]
    x[0] += int(r["n_cpg"])
    x[1] += int(r["n_meth"])


def arr(g, br, part):
    rows = []
    for v in d[(g, br)].values():
        if part == "전체":
            c, m = v["코딩부"][0] + v["스페이서"][0], v["코딩부"][1] + v["스페이서"][1]
        else:
            c, m = v[part]
        if c:
            rows.append((c, m))
    return np.array(rows, dtype=float).reshape(-1, 2)


def boot(a):
    ix = rng.integers(0, len(a), size=(NB, len(a)))
    return a[ix, 1].sum(1) / np.maximum(a[ix, 0].sum(1), 1)


def cite(n):
    return "인용금지" if n < 10 else ("약함" if n < 30 else "가능")


o = io.open(D + "WGBS_CHR18_VS_CHR4.tsv", "w", encoding="utf8", newline="\n")
o.write(T.join(["sample", "label", "window", "source", "part", "n_read_chr18", "n_read_chr4", "chr18", "chr4",
                "diff_chr18_minus_chr4", "lo95", "hi95", "cite"]) + "\n")
for g, (lab, win, src) in SRC.items():
    for part in ("전체", "코딩부", "스페이서"):
        a, b = arr(g, "chr18", part), arr(g, "chr4", part)
        if len(a) == 0 or len(b) == 0:
            o.write(T.join([g, lab, win, src, part, str(len(a)), str(len(b)), "-", "-", "-", "-", "-", "read없음"]) + "\n")
            continue
        fa, fb = a[:, 1].sum() / a[:, 0].sum(), b[:, 1].sum() / b[:, 0].sum()
        if len(a) >= 2:
            dd = boot(a) - boot(b)
            lo, hi = "%+.3f" % np.percentile(dd, 2.5), "%+.3f" % np.percentile(dd, 97.5)
        else:
            lo = hi = "-"
        o.write(T.join([g, lab, win, src, part, str(len(a)), str(len(b)), "%.3f" % fa, "%.3f" % fb, "%+.3f" % (fa - fb),
                        lo, hi, cite(len(a))]) + "\n")
# Stage III 난모세포 : CpG 호출 합산(자리 단위 · read 없음 · 구간 없음)
st = collections.defaultdict(lambda: [0, 0])
for r in csv.DictReader(io.open(R + "12_zebrafish_devseries/oocyte_GSE175951/STEM_5S.tsv", encoding="utf8"), delimiter=T):
    if r["part"] in ("코딩부", "스페이서"):
        for p in (r["part"], "전체"):
            k = (r["sample"], r["locus"], p)
            st[k][0] += int(r["n_call"])
            st[k][1] += int(r["n_meth"])
for s in ("oocyte_StageIII_rep1", "oocyte_StageIII_rep2"):
    for part in ("전체", "코딩부", "스페이서"):
        a, b = st[(s, "chr18_somatic", part)], st[(s, "chr4_maternal", part)]
        fa, fb = a[1] / a[0], b[1] / b[0]
        o.write(T.join([s, "Stage III 난모세포", "난자 형성기(Stage III)", "GSE175951(저자 가공 CpG 표)", part,
                        "호출%d" % a[0], "호출%d" % b[0], "%.3f" % fa, "%.3f" % fb, "%+.3f" % (fa - fb), "-", "-",
                        "자리단위·read없음"]) + "\n")
o.close()
print("완료")
