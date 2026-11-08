# -*- coding: utf-8 -*-
# Fig 3 (G) 의 입력표를 만든다 — FIG3_G_TYPEDIST.tsv
#   (2026-11-04 에 새로 넣은 칸 · 앞 이름 r96_make_table.py)
#
#   무엇을 내나
#     점 하나 = 한 하플로타입. 같은 array-type 군 구성원까지의 평균 거리와
#     다른 군까지의 평균 거리를 한 줄에 함께 적는다.
#
#   코호트 : H3 = H1 ∩ H2 (변이로 묶은 군과 메틸을 함께 쓰는 주장이므로)
#   관문   : L2_HAP_UNIT 을 n_read >= 5 로 거른 뒤, 다섯 구역 모두에 값이 있는 하플로타입만
#   거리   : 두 하플로타입의 다섯 구역 beta 의 |delta| 평균
#   제외   : 같은 개체의 짝은 거리 계산에서 뺀다 (비독립)
#   군     : 무리 크기 >= 5 이고 무리 번호가 0 이 아닌 것만
#   난수   : default_rng(20261104) · 군 라벨 치환 2,000 회 (표에는 쓰지 않고 화면에만 찍는다)
#
#   입력 넷 (절대 경로는 아래 ROOT 와 GRP 로 모은다)
#     ① <ROOT>/02_tables/TIERS_H0H3_SLOTS_2026-10-05.tsv   층 (H3 를 고르는 데 쓴다)
#     ② <GRP>   5S-rDNA-paper-en/tables/X12_HAP_GROUPS.tsv
#               Fig 2M 의 군 (single linkage · Jaccard >= 0.90 · 무리 크기 >= 5)
#     ③ <ROOT>/02_tables/L2_HAP_UNIT.parquet               하플로타입x유닛 메틸 (열 일곱 개만 읽는다)
#     ④ <ROOT>/02_tables/POP_COHORT.tsv                    초집단과 유닛 수
#
#   내는 곳 : <FIG>/FIG3_G_TYPEDIST.tsv
#             환경변수 FIG3G_OUT 을 주면 그 자리에 쓴다 (검산용 · 발행 입력표를 안 건드리고 재현)
#
#   이 표를 읽는 글 : fig3_main_draw_en.py 의 (G) 칸 · fig3_sourcedata_en.py 의 G_type_distance 시트
#   그림이 찍는 p 는 짝 Wilcoxon 이다 (치환 p 가 아니다).
import hashlib
import io
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

sys.dont_write_bytecode = True
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
S = ""
NL = chr(10)
ROOT = os.environ["DATA_ROOT"]
FIG = ROOT + "/03_figures/paper_Fig"
GRP = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "tables", "X12_HAP_GROUPS.tsv")
md5 = lambda p: hashlib.md5(io.open(p, "rb").read()).hexdigest()[:12]
rng = np.random.default_rng(20261104)

TI = pd.read_csv(ROOT + "/02_tables/TIERS_H0H3_SLOTS_2026-10-05.tsv", sep="\t", dtype=str)
TI["hap"] = TI["sample"].astype(str) + "_" + TI.slot.astype(str)
H3 = set(TI[TI.tier == "H3"].hap)

G = pd.read_csv(GRP, sep="\t")
G.columns = ["clu", "hap", "sample", "slot", "clu_size", "n_token"]
big = G[(G.clu_size >= 5) & (G.clu != 0)]

L2 = pd.read_parquet(ROOT + "/02_tables/L2_HAP_UNIT.parquet",
                     columns=["sample", "slot", "ap", "ad", "N", "n_read", "b_mean"])
L2 = L2[L2.n_read >= 5].copy()
L2["hap"] = L2["sample"].astype(str) + "_" + L2.slot.astype(str)
Z = [("P1_3", (L2.ap <= L2.ad) & (L2.ap <= 3)),
     ("P4_10", (L2.ap <= L2.ad) & (L2.ap >= 4) & (L2.ap <= 10)),
     ("INT", ((L2.ap <= L2.ad) & (L2.ap > 10)) | ((L2.ad < L2.ap) & (L2.ad > 21))),
     ("D11_21", (L2.ad < L2.ap) & (L2.ad >= 11) & (L2.ad <= 21)),
     ("D1_10", (L2.ad < L2.ap) & (L2.ad <= 10))]
V = pd.DataFrame(index=sorted(L2.hap.unique()))
for nm, msk in Z:
    V[nm] = L2[msk].groupby("hap").b_mean.mean()
V = V.dropna()
P = pd.read_csv(ROOT + "/02_tables/POP_COHORT.tsv", sep="\t")
P["hap"] = P["sample"].astype(str) + "_" + P.slot.astype(str)
P = P.set_index("hap")[["sample", "sp", "N"]]
D = big.set_index("hap").join(V, how="inner").join(P[["sp", "N"]], how="inner")
D = D.rename(columns={"sample": "ind"})
D = D[[h in H3 for h in D.index]]
COLS = ["P1_3", "P4_10", "INT", "D11_21", "D1_10"]
X = D[COLS].to_numpy(float)
clu = D.clu.astype(int).to_numpy()
ind = D.ind.astype(str).to_numpy()
hap = np.array(D.index)
n = len(hap)
print("%s코호트 H3 · 하플로타입 %d · 개체 %d · 군 %d%s"
      % (S, n, len(set(ind)), len(set(clu)), S))

DIST = np.abs(X[:, None, :] - X[None, :, :]).mean(axis=2)
OK = (ind[:, None] != ind[None, :])
np.fill_diagonal(OK, False)
same = clu[:, None] == clu[None, :]
rows = []
for i in range(n):
    w = DIST[i][same[i] & OK[i]]
    b = DIST[i][(~same[i]) & OK[i]]
    if len(w) == 0 or len(b) == 0:
        continue
    rows.append(dict(hap=hap[i], 개체=ind[i], 무리=clu[i], 무리크기=int((clu == clu[i]).sum()),
                     초집단=D.sp.iloc[i], 유닛수=int(D.N.iloc[i]),
                     n_군안=len(w), 군안거리=round(float(w.mean()), 5),
                     n_군밖=len(b), 군밖거리=round(float(b.mean()), 5),
                     차=round(float(w.mean() - b.mean()), 5)))
T = pd.DataFrame(rows)
print("  점(하플로타입) %d · 군 안 쪽이 더 가까운 것 %d (%.1f%%)"
      % (len(T), int((T.차 < 0).sum()), 100.0 * (T.차 < 0).mean()))
st, pw = wilcoxon(T.군안거리, T.군밖거리)
print("  짝 윌콕슨 p = %.3g · 짝 안 중앙 차 %+.4f"
      % (pw, float(np.median(T.차))))
print("  중앙값 : 군 안 %.4f · 군 밖 %.4f"
      % (T.군안거리.median(), T.군밖거리.median()))

# 군 라벨 치환 귀무 (하플로타입 수준 통계로)
obs = float(np.median(T.차))
cnt = 0
NP = 2000
idx = {h: k for k, h in enumerate(hap)}
for _ in range(NP):
    p = clu.copy()
    rng.shuffle(p)
    sm = p[:, None] == p[None, :]
    dif = []
    for i in range(n):
        w = DIST[i][sm[i] & OK[i]]
        b = DIST[i][(~sm[i]) & OK[i]]
        if len(w) and len(b):
            dif.append(w.mean() - b.mean())
    if np.median(dif) <= obs:
        cnt += 1
print("  군 라벨 치환 %d 회 : p = %.4f (바닥 %.4f)" % (NP, (cnt + 1) / (NP + 1), 1 / (NP + 1)))

OUT = os.environ.get("FIG3G_OUT") or os.path.join(FIG, "FIG3_G_TYPEDIST.tsv")
HDR = [
    "#  Fig 3 (G) 입력 — array-type 군 안 ↔ 군 밖 메틸 거리 (2026-11-04)",
    "#  코호트 : H3 = H1 ∩ H2 (변이 참조 ∩ 메틸 참조) · 하플로타입 %d · 개체 %d · 군 14" % (
        len(T), T["개체"].nunique()),
    "#  군 : X12_HAP_GROUPS.tsv (Fig 2M · single linkage Jaccard >= 0.90 · 무리 크기 >= 5)",
    "#  메틸 : L2_HAP_UNIT.parquet (n_read >= 5) 의 다섯 구역 평균 beta",
    "#  거리 : 두 하플로타입의 다섯 구역 |delta beta| 평균 · 같은 개체 짝은 뺀다",
    "#  점 하나 = 한 하플로타입 (군 안 평균 거리 · 군 밖 평균 거리)",
]
#  출처 머리글은 표 안에 쓰지 않는다 — 이 폴더의 그림 입력표 넷(FIG3_G_TYPEDIST ·
#  FIG5V2_PAIRZONE5 · FIG8_PAIRS · NORM_VS_CANCER_SAMPLE)은 모두 `#` 줄이 없고,
#  읽는 글 둘이 `comment="#"` 없이 읽는다
#  (`fig3_main_draw_en.py:103` 과 `fig3_sourcedata_en.py:36` 이 그냥 `pd.read_csv(sep="\t")`).
#  ⇒ `#` 줄을 넣으면 그림이 그것을 열 이름으로 읽어 깨진다 (2026-11-05 실측으로 잡았다).
#  그래서 머리글은 곁 파일 로 따로 쓴다. 표는 열 이름 줄부터 시작한다.
with io.open(OUT, "w", encoding="utf-8", newline=NL) as f:
    T.to_csv(f, sep="\t", index=False, lineterminator=NL)
with io.open(OUT + ".about.txt", "w", encoding="utf-8", newline=NL) as f:
    f.write(NL.join(HDR).replace("#  ", "") + NL)
print()
print("%s표 %s · %s · %d 줄%s" % (S, os.path.basename(OUT), md5(OUT), len(T), S))
print("   곡 파일 : %s.about.txt" % os.path.basename(OUT))
print(S + "끝" + S)
