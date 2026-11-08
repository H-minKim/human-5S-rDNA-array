# -*- coding: utf-8 -*-
r"""Source data Excel files (English) for Figure 2 and Supplementary Figure 2. of 2026-10-20:
   'put the result part of the summary table into Excel form (in English), and when submitting the paper, provide the source Excel of every data item'
   - Values: the values actually drawn by the published panel figures (English edition), i.e. the JSON that srcdata.py read from the Figure just before saving (SRCDATA/<panel>__en.json)
   - Panel letter, panel content, denominator, aggregation unit and cut: from the Fig 2 layout table (Korean-named file, in Korean), translated to English with the dictionary below (stops if a term is missing from the dictionary, and the multiset of numbers in the original and in the English must be the same)
   - Output: <OUTD>/Fig2_SourceData.xlsx , FigS2a_SourceData.xlsx , FigS2b_SourceData.xlsx (a README sheet plus one sheet per panel), with a check that no Korean is left
   - the supplement is split in two (layout table SFig2 a-i -> FigS2a, j-s -> FigS2b, sheet names FigS2aA ... FigS2bS)
     - the location of the drawn-value JSON can be changed with the environment variable SRCDATA_DIR (default HERE/SRCDATA)
   - 2026-11-01 user 'make the S2 naming a single S2': the supplement is one figure over two pages (Supplementary Figure S2, panels A-Q) =>
     output = Fig2_SourceData.xlsx , FigS2_SourceData.xlsx (sheets FigS2A ... FigS2Q). The page split (a-g = page 1) appears only in the README sheet"""
import io, os, re, sys, json, glob, hashlib, collections, datetime
sys.stdout.reconfigure(encoding="utf-8")
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment, PatternFill
from openpyxl.utils import get_column_letter
HERE = os.path.dirname(os.path.abspath(__file__))
P = (__import__("os").environ["RDNA_DATA_ROOT"] + r"/for_paper/00_paper_figures_distal_layout")
SRC = os.environ.get("SRCDATA_DIR") or os.path.join(HERE, "SRCDATA")
OUTD = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "SOURCEDATA_OUT")
os.makedirs(OUTD, exist_ok=True)
md5 = lambda q: hashlib.md5(open(q, "rb").read()).hexdigest()
HAN = re.compile("[가-힣]")

# ---------------- Korean to English dictionary for the layout table (2026-10-20, moved by ) ----------------
TR = {
 # (second) main Fig 2 A: blue line = all substitutions plus deletions and insertions outside blocks, orange = deletions and insertions inside blocks (old keys are kept to reproduce old editions)
 "개별 자리(bp)마다, 배열 전체 평균 — 치환과 블록 밖 결실 삽입은 염기 단위, 반복 블록 안 결실 삽입은 블록 폭": "Per position (bp), averaged over the whole array; substitutions, and insertions or deletions outside repeat blocks, per base; insertions or deletions within repeat blocks drawn as block width",
 "유닛 (치환과 블록 밖 결실 삽입은 자리별 · 블록 안 결실 삽입은 블록별 · 한 유닛은 한 번만)": "Unit (per position for substitutions and for insertions or deletions outside repeat blocks, per block for insertions or deletions within them; each unit counted once)",
 # layout table text after the terminal-type definition (10 or more of 12 diagnostic variants) (old keys are kept to reproduce old editions)
 "말단형 사이 간격 분포 (간격 1 은 3 건뿐, 0.0071)": "Distribution of gaps between terminal-type units (gap 1 occurs only 3 times, 0.0071)",
 "말단형 보유 패턴 × 집단, 패턴 69 개 전부": "Terminal-unit carrier pattern by population, all 69 patterns",
 "집단 표기 하플로타입 188 (미표기 2 제외)": "188 haplotypes with a population label (2 unlabelled excluded)",
 "primary 564 (개체 305 · H1 으로 다시 만들지 않은 세대) · 이웃 짝 말단형-보통 1,229 · 보통-보통 원위 6,660 · 가운데 6,254 · 근위 7,892 · 말단형-말단형 3 (5 쌍 미만이라 그리지 않음) ‖ 2026-10-21 : 옛 판 'H1 557' 은 틀렸다 — ZB_neighbour_divergence 는 FIGTAB_H1_출처 에 '564 세대 그대로' · Z_B 주석 말단형 없음 372 + 보유 192 = 564 · load('primary') 실측 564/305":
     "primary 564 (305 individuals; generation not rebuilt on H1); neighbour pairs: terminal-type vs ordinary 1,229; ordinary vs ordinary, distal 6,660, central 6,254, proximal 7,892; terminal-type vs terminal-type 3 (fewer than 5 pairs, not drawn); the cohort is the 564 generation, which ZB_neighbour_divergence uses unchanged: Z_B notes give no terminal type 372 + carriers 192 = 564, and load('primary') measures 564/305",
 # panel content (panel title)
 "개별 자리(bp)마다, 배열 전체 평균 — 블록 밖은 염기 단위, 반복 블록은 블록 폭": "Per position (bp), averaged over the whole array; outside repeat blocks per base, repeat blocks drawn as block width",
 "array 위치별 변이량 (다수형 기준)": "Variant load by array position (majority-unit reference)",
 "기원별 누적 (진단 변이 분리 없음, 검은 선은 총 변이 수)": "Stacked by origin class (diagnostic variants not separated; black line = total variant count)",
 "전부, 유닛당 변이 구성 (11 종 + 기타)": "All units and single-unit variants only: variant composition per unit (11 classes + other)",
 "구간 나눔표 — 근위에서 원위까지": "Zoning map from the proximal to the distal end",
 "깊이 2 에서 시작한다 — 변이 봉우리도 같은 깊이에 있다 (겹쳐 그림)": "Terminal units start at depth 2; the variant peak sits at the same depth (overlaid)",
 "말단형 사이 간격 분포 (간격 1 은 6 건뿐, 0.0138)": "Distribution of gaps between terminal-type units (gap 1 occurs only 6 times, 0.0138)",
 "말단형↔비말단형만 안 나눈다 (간격 ≤2)": "Only terminal and non-terminal units do not exchange (gap ≤2)",
 "집단마다 크게 다르다, AFR 이 가장 낮다": "Carrier fraction differs strongly between populations; AFR is lowest",
 "닮은 상대를 가진 유닛의 몫 — array 내 ↔ 상동 ↔ 무관 (집단 4 단계)": "Fraction of units with a matching partner: within array, homologous, unrelated (4 population levels)",
 "거리별 쌍 평균 차이 + 전역 봉투 띠 (전체 · 발산 사분위 제외 · 곡선마다 자기 띠)": "Mean pairwise difference by separation, with global envelope bands (all pairs and divergent quartile removed; each curve has its own band)",
 "패치 길이 추정 (구간검열 최대가능도) + 교환 없는 바닥 187.7 bp": "Patch length estimate (maximum-likelihood tract length) with the no-exchange floor of 187.7 bp",
 "패치 길이 추정 (구간검열 최대가능도) + 교환 없는 바닥 190.1 bp": "Patch length estimate (maximum-likelihood tract length) with the no-exchange floor of 190.1 bp",
 "패치 길이 추정 (구간검열 최대가능도) + 교환 없는 바닥 191.9 bp": "Patch length estimate (maximum-likelihood tract length) with the no-exchange floor of 191.9 bp",
 "근위부, 영역(문턱 넘는 몫)": "Proximal side: excess over threshold by region",
 "원위부, 영역(문턱 넘는 몫)": "Distal side: excess over threshold by region",
 "array 영역별": "By array region",
 "길어질 때 느는 것은 가운데 구획뿐": "As the array lengthens, only the central compartment grows",
 "드문 층에서 가장 큰 신호 — 원위에서 CpG 탈아미노가 빠진다 (메틸 β 선 없음)": "Strongest signal in the rare stratum: CpG deamination is depleted at the distal end (no methylation beta line)",
 "무리마다 집단 구성 — 집단 못 찾음까지": "Population make-up of each array-type cluster, including haplotypes without a population label",
 "'그 집단에만 있는 토큰' 이 무엇인가 — 200 회 내려뽑기": "Population-private tokens, 200 downsampling rounds",
 "세포주에서도 같은 규칙 (붉은 점 = 말단형, 간격 1 없음)": "The same rule holds in cell lines (red points = terminal-type units; no gap of 1)",
 "영역마다 거리에 따라 — 그 영역이 통째로 같을 확률": "By region and distance: probability that the whole region is identical",
 "구역마다 흔한 순서, 층을 바꾸면 순서가 바뀐다": "Rank of variant classes in each zone; the order changes between strata",
 "왜 거기에 경계를 두었나 — 분할 이득과 벌점 문턱": "Why the boundaries sit where they do: split gain and the penalty threshold",
 "경계를 가르는 것은 거의 전부 다개체 공유형": "Variants separating boundaries are almost all multi-individual shared",
 "거의 항상 깊이 2 에서 시작하고 2-3 칸 간격 (HPRC 557 하플로타입)": "Terminal units almost always start at depth 2 and recur at gaps of 2-3 (HPRC, 557 haplotypes)",
 "말단형을 정의하는 변이 (상위 14) 말단형 묶음 첫 그림": "Variants that define the terminal-type unit (top 14)",
 "말단형만의 성질이다, 다른 소수형은 간격 1 이 오히려 많다": "The spacing is specific to terminal units; other minor types show more gaps of 1",
 "사람이 달라도 닮는다 (사용자가 말한 '말단형은 하플로간에 닮음')": "Terminal units are similar even between different people",
 "말단형 보유 패턴 × 집단, 패턴 70 개 전부": "Terminal-unit carrier pattern by population, all 70 patterns",
 "array 유형이란 무엇인가 : 몇 개의 변이 토큰을 array 전체에 깔고 있는 하플로타입 묶음": "What an array type is: a group of haplotypes that carry a few variant tokens across the whole array",
 "개별 변이별": "By individual variant",
 "삼자 전달을 유닛 자리까지 맞춰 본다 — 전달 5 짝은 447 유닛 중 436 이 자리까지 같다": "Trio transmission matched to unit position: in the 5 transmitted pairs, 436 of 447 units are identical including position",
 "말단형과 이웃 보통 유닛은 약 2% 다르다": "Terminal units differ from their neighbouring ordinary units by about 2%",
 "희귀 변이를 함께 가진 두 유닛은 다른 변이도 함께 가진다 — 다른 변이가 희귀할 때 · 흔할 때": "Two units sharing a rare variant also share other variants, when the other variant is rare and when it is common",
 # denominator
 "A 거리마다 유닛쌍 2,000 (값 눈금 0.0005) · B 가운데 창 0.2–0.8 하플로타입 534 · 유닛 30,268": "A: 2,000 unit pairs per distance (value resolution 0.0005); B: central window 0.2 to 0.8, 534 haplotypes, 30,268 units",
 "H1 557 · 이웃 짝 말단형-보통 1,250 · 보통-보통 원위 6,636 · 가운데 6,254 · 근위 7,892 · 말단형-말단형 6": "H1 set, 557 haplotypes; neighbouring pairs: terminal-ordinary 1,250, ordinary-ordinary distal 6,636, middle 6,254, proximal 7,892, terminal-terminal 6",
 "H1 557 · 전달 5 짝(유닛 447 · 부모와 유닛 수 같음) · 비전달 5 · 무관 150 (짝마다 정렬 최선 하나)": "H1 set, 557 haplotypes; transmitted 5 pairs (447 units, same unit count as the parent), non-transmitted 5, unrelated 150 (best alignment per pair)",
 "H1 557 가운데 창(0.2-0.8 · 창 유닛 23 이상) 하플로타입 533 · 희귀 변이 A(코호트 유닛 < 1%) 를 함께 가진 쌍 108,529 · 같은 거리에서 한쪽만 가진 쌍 95,804": "H1 set (557); central window (0.2 to 0.8, at least 23 window units), 533 haplotypes; pairs sharing rare variant A (cohort unit frequency < 1%) 108,529; pairs at the same distance where only one unit carries A 95,804",
 "H1 557 슬롯 (개체 301)": "H1 set, 557 haplotypes (301 individuals)",
 "H1 557 슬롯 (개체 301) · 기원 f 분모 555": "H1 set, 557 haplotypes (301 individuals); origin f denominator 555",
 "거리마다 유닛쌍 4,000 · 사건 14,144~16,016 (X4) · 바닥은 무관 쌍 4,000 (X9)": "4,000 unit pairs per distance; 14,144 to 16,016 events (table X4); floor from 4,000 unrelated pairs (table X9)",
 "내려뽑기 200 회 · 검정 가능 토큰 중앙 1,372": "200 downsampling rounds; median 1,372 testable tokens",
 "무리 13 (하플로타입 5 개 이상) · 하플로타입 143": "13 clusters (5 or more haplotypes); 143 haplotypes",
 "세포주 (그림에 적힘)": "Cell lines (named in the panel)",
 "유닛 45,383 (interior · 하플로타입 546 / H1 557 중) ‖ 비교 집합은 모두 하플로타입 하나 크기로 맞춤(5 회 평균)": "45,383 interior units (546 of the 557 H1 haplotypes); every comparison set is matched to the size of one haplotype (mean of 5 draws)",
 "집단 표기 하플로타입 189 (미표기 2 제외)": "189 haplotypes with a population label (2 unlabelled excluded)",
 "하플로 542 (S18 basis B · pen 1.0) ‖ 구역은 깊이 프로파일 (MINHAP 100 · DMAX 40)": "542 haplotypes (table S18, criterion B, penalty 1.0); zones from the depth profile (minimum 100 haplotypes, maximum depth 40)",
 "하플로타입 510 · 개체 280 (무리 30 · 하플로타입 5 개 이상 13 을 그림)": "510 haplotypes, 280 individuals (30 clusters; the 13 with 5 or more haplotypes are drawn)",
 "하플로타입 534 · 유닛 30,268 (가운데 창 0.2–0.8)": "534 haplotypes, 30,268 units (central window 0.2 to 0.8)",
 # denominator text revised in the layout table on 2026-10-21 (supplement B)
 "하플로타입 534 (가운데 창 0.2–0.8 · Q2b·XB2) · 거리마다 유닛쌍 22,280~29,756 (XB2 · 거리 1–15) ‖ 2026-10-21 : 옛 판 '유닛 30,268' 은 X5·X6(하플로타입 533)의 수라 지웠다 · XB2 쌍 수로 거꾸로 세면 30,290 (쌍 세는 법 가정 · 미확정)":
     "534 haplotypes (central window 0.2 to 0.8; tables Q2b and XB2); 22,280 to 29,756 unit pairs per separation (XB2; separations 1 to 15); counting back from the XB2 pair numbers gives 30,290 units, which assumes a pair-counting rule and is not confirmed",
 "거리마다 유닛쌍 4,000 · 사건 13,727~16,016 (X4 · 거리 1·2·3·5·8·12·18·25 · 최소는 거리 2) · 바닥은 무관 쌍 4,000 · 사건 19,384 (X9) ‖ 2026-10-21 : 옛 판 14,144 는 거리 1 값이었다":
     "4,000 unit pairs per separation; 13,727 to 16,016 events (X4; separations 1, 2, 3, 5, 8, 12, 18 and 25; the minimum is at separation 2); floor from 4,000 unrelated pairs, 19,384 events (X9)",
 "A 거리마다 유닛쌍 2,000 (값 눈금 0.0005 · P08_X2.py NP_=2000) · B 가운데 창 0.2–0.8 하플로타입 534 (Q2b) ‖ 2026-10-21 : 옛 판 '유닛 30,268' 은 X5·X6(하플로타입 533)의 수라 지웠다":
     "A: 2,000 unit pairs per separation (value resolution 0.0005; P08_X2.py NP_=2000); B: central window 0.2 to 0.8, 534 haplotypes (Q2b)",
 "primary 564 (개체 305 · H1 으로 다시 만들지 않은 세대) · 이웃 짝 말단형-보통 1,250 · 보통-보통 원위 6,636 · 가운데 6,254 · 근위 7,892 · 말단형-말단형 6 ‖ 2026-10-21 : 옛 판 'H1 557' 은 틀렸다 — ZB_neighbour_divergence 는 FIGTAB_H1_출처 에 '564 세대 그대로' · Z_B 주석 말단형 없음 372 + 보유 192 = 564 · load('primary') 실측 564/305":
     "primary 564 (305 individuals; generation not rebuilt on H1); neighbour pairs: terminal-type vs ordinary 1,250; ordinary vs ordinary, distal 6,636, central 6,254, proximal 7,892; terminal-type vs terminal-type 6; the cohort is the 564 generation, which ZB_neighbour_divergence uses unchanged: Z_B notes give no terminal type 372 + carriers 192 = 564, and load('primary') measures 564/305",
 "하플로타입 542 → 깊이 프로파일 (MINHAP 100 · DMAX 40)": "542 haplotypes, summarised as a depth profile (minimum 100 haplotypes, maximum depth 40)",
 "하플로타입 549 · 개체 297": "549 haplotypes, 297 individuals",
 # aggregation unit
 "A 하플로타입 (깊이별 경계 수) · B 구역": "A: haplotype (boundary counts by depth); B: zone",
 "같은 array 안 유닛쌍 (거리 k · 거리마다 N_k 22,280~29,756)": "Unit pairs within the same array (separation k; 22,280 to 29,756 pairs per separation)",
 "같은 array 안 유닛쌍 (거리 k · 거리마다 N_k 22,223~29,657)": "Unit pairs within the same array (separation k; 22,223 to 29,657 pairs per separation)",
 "같은 array 안 유닛쌍 (거리 k · 거리마다 N_k 22,593~30,167)": "Unit pairs within the same array (separation k; 22,593 to 30,167 pairs per separation)",
 "경계 1 개 (경계를 가르는 변이 수 · 귀무 = 같은 하플로타입의 경계 아닌 무작위 자리)": "1 boundary (number of variants separating the boundary; null = random non-boundary positions in the same haplotype)",
 "깊이 (자르는 자리)": "Depth (cut position)",
 "말단형 유닛쌍 (간격)": "Pairs of terminal-type units (gap)",
 "무리 (집단 구성 비율)": "Cluster (population composition)",
 # 4-zone windows, permutation test, S2aF by zone (layout table MK_TAB_1029)
 "무리 14 (하플로타입 5 개 이상) · 하플로타입 295 · 코호트 하플로타입 510 · 개체 280": "14 clusters (5 or more haplotypes); 295 haplotypes; cohort 510 haplotypes, 280 individuals",
 "내려뽑기 200 회 · 검정 가능 토큰 중앙 1,372 · 치환 1,000 회 P 0.001 · 하플로타입 503 · 개체 276": "200 downsampling rounds; median 1,372 testable tokens; 1,000 permutations of population labels between individuals, P 0.001; 503 haplotypes, 276 individuals",
 "하플로타입 510 · 개체 280 (무리 27 · 하플로타입 5 개 이상 14 를 그림)": "510 haplotypes, 280 individuals (27 clusters; the 14 with 5 or more haplotypes are drawn)",
 "무리 × 토큰 (값 = 서명으로 잡힌 구역 수 ÷ 4)": "Cluster by token (value = number of zones in which the token separates the cluster, divided by 4)",
 "경계 1 개 (경계를 가르는 변이 수 · 귀무 = 같은 하플로타입의 경계 아닌 무작위 자리) · 경계 자리 = 구역 4 (근위 1-3 · 4-17 · 가운데 · 원위 1-10)":
     "1 boundary (number of variants separating the boundary; null = random non-boundary positions in the same haplotype); boundary position = 4 zones (proximal 1-3, 4-17, central, distal 1-10)",
 "무리 × 토큰 (값 = 서명으로 잡힌 구역 창 수 ÷ 9)": "Cluster by token (value = number of zone windows in which the token separates the cluster, divided by 9)",
 "변이 (상위 14)": "Variant (top 14)",
 "유닛": "Unit",
 "유닛 (51,801)": "Unit (51,801)",
 "유닛 (array 길이 N 층 × array 안 자리 · 값표 '근위다움 · N<=72 · P1 … (유닛 146…)')": "Unit (array-length stratum N by position in the array)",
 "유닛 (array 길이 N 층 × array 안 자리 · 값표와 같음)": "Unit (array-length stratum N by position in the array)",
 "유닛 (구역 × 층)": "Unit (zone by stratum)",
 "유닛 (닮은 상대 있음/없음)": "Unit (with or without a matching partner)",
 "유닛 (블록 밖은 자리별 · 블록 안은 블록별 · 한 유닛은 한 번만)": "Unit (per position outside repeat blocks, per block inside them; each unit counted once)",
 "유닛 (짝지은 자리)": "Unit (paired position)",
 "유닛 쌍 (구간은 하플로타입 cluster bootstrap 500)": "Unit pair (intervals from a haplotype cluster bootstrap with 500 replicates)",
 "유닛 자리 (61 칸 · 가운데는 양쪽 깊이 >30 압축)": "Unit position (61 positions; the middle, depth >30 from both ends, is compressed)",
 "유닛쌍": "Unit pair",
 "유닛쌍 (간격)": "Unit pair (gap)",
 "유닛쌍 × 영역 (영역이 통째로 같은가)": "Unit pair by region (whether the whole region is identical)",
 "유닛쌍 사이 구간 (구간검열 · 마커 간격 105.5 bp 대리)": "Interval between a unit pair (marker spacing 105.5 bp as proxy)",
 "유닛쌍 사이 구간 (구간검열 · 마커 간격 105.7 bp 대리)": "Interval between a unit pair (marker spacing 105.7 bp as proxy)",
 "유닛쌍 사이 구간 (구간검열 · 마커 간격 105.9 bp 대리)": "Interval between a unit pair (marker spacing 105.9 bp as proxy)",
 "이웃한 두 유닛": "Two neighbouring units",
 "집단 (내려뽑기 회차별 토큰 수)": "Population (token count per downsampling round)",
 "하플로타입": "Haplotype",
 "하플로타입 (N 칸별 · 값표 'N칸 1 | … 하플로 112')": "Haplotype (by array-length bin N)",
 "하플로타입 (N 칸별 × 개별 변이 · 값표 '비례 개별변이 2205TC … 하플로 111')": "Haplotype (by array-length bin N and individual variant)",
 "하플로타입 (깊이별 보유) · 겹친 선은 유닛 (깊이별 평균 변이 수)": "Haplotype (carrier status by depth); the overlaid line is per unit (mean variant count by depth)",
 "하플로타입 (패턴별)": "Haplotype (by pattern)",
 "하플로타입 · 개체 둘 다": "Both haplotype and individual",
 "하플로타입쌍": "Haplotype pair",
 # scale, cut, direction
 "다른 변이 B 의 공유 = 두 유닛 모두 같은 대립 · 후보 = 둘 중 하나라도 가진 B": "Sharing of another variant B = both units carry the same allele; candidates = B carried by at least one of the two units",
 "가르는 토큰 : 무리안 보유율 ≥0.80 & 무리안−무리밖 ≥0.30": "Separating token: in-cluster carrier rate ≥0.80 and in-cluster minus out-of-cluster ≥0.30",
 "공유형 뺌 (f < 0.10)": "Multi-individual shared type removed (f < 0.10)",
 "희귀 = 코호트 유닛 빈도 < 0.01 (A · 패널 A 의 B) · 흔함 = ≥ 0.10 (패널 B 의 B) · 100 bp 이하 구간 제외": "Rare = cohort unit frequency < 0.01 (A, and B in panel A); common = ≥ 0.10 (B in panel B); distances of 100 bp or less excluded",
 "array_pos 1 = 원위": "array_pos 1 = distal end",
 "array_pos 1 = 원위 (U84C=자리 2205 실측 : 근위끝 0.1961 · 원위끝 0.0000)": "array_pos 1 = distal end (U84C at position 2205, measured: proximal end 0.1961, distal end 0.0000)",
 "L2.ap = 근위에서 센 깊이": "Depth counted from the proximal end (L2.ap)",
 "해당없음": "Not applicable",
}
# numbers deliberately dropped in translation (example numbers inside value-table names, value-table citations in parentheses in the original), only those with a reason for each panel
DROP_OK = {
 "유닛 (array 길이 N 층 × array 안 자리 · 값표 '근위다움 · N<=72 · P1 … (유닛 146…)')": {"72", "1", "146"},
 "하플로타입 (N 칸별 · 값표 'N칸 1 | … 하플로 112')": {"1", "112"},
 "하플로타입 (N 칸별 × 개별 변이 · 값표 '비례 개별변이 2205TC … 하플로 111')": {"2205", "111"},
 "말단형을 정의하는 변이 (상위 14) 말단형 묶음 첫 그림": set(),
}
NUM = re.compile(r"\d+(?:,\d{3})*(?:\.\d+)?")


def nums(s): return collections.Counter(x.replace(",", "") for x in NUM.findall(s))


def en(s):
    s = s.strip()
    if not HAN.search(s): return s
    if s not in TR: raise SystemExit("사전에 없음 : %r" % s)
    t = TR[s]
    a, b = nums(s), nums(t)
    miss = a - b
    miss = collections.Counter({k: v for k, v in miss.items() if k not in DROP_OK.get(s, set())})
    extra = b - a
    if miss or extra: raise SystemExit("숫자 어긋남 %r\n  원문 %s\n  영어 %s\n  빠짐 %s 더함 %s" % (s, dict(a), dict(b), dict(miss), dict(extra)))
    assert not HAN.search(t), t
    return t


def rows(p):
    L = [l.rstrip("\n").split("\t") for l in io.open(p, encoding="utf8") if "\t" in l and not l.startswith("#")]
    return {c: i for i, c in enumerate(L[0])}, L[1:]


H, R = rows(os.path.join(P, "tables", "FIG2_layout_2026-10-09.tsv"))
PANELS = {"Fig2": [], "SFig2": []}
for r in R:
    g = r[H["그림"]]
    if g in PANELS and r[H["본문칸그림"]].strip():
        PANELS[g].append(r)
for k in PANELS: PANELS[k].sort(key=lambda r: r[H["칸"]])
assert [len(PANELS[k]) for k in ("Fig2", "SFig2")] == [14, 17], {k: len(v) for k, v in PANELS.items()}
# supplement is one figure S2, pages = two composites (page 1 a-g = fig2_suppa, page 2 h-q = fig2_suppb, new layout of 2026-10-29)
PAGE1 = [r[H["칸"]].upper() for r in PANELS["SFig2"] if r[H["칸"]] <= "g"]
PAGE2 = [r[H["칸"]].upper() for r in PANELS["SFig2"] if r[H["칸"]] > "g"]
assert (len(PAGE1), len(PAGE2)) == (7, 10), (PAGE1, PAGE2)
TITLE = {"Fig2": ("Figure 2", "Fig2_SourceData.xlsx", [("fig2_main_en", "")], "Fig2"),
         "SFig2": ("Supplementary Figure S2", "FigS2_SourceData.xlsx",
                   [("fig2_suppa_en", " (page 1, panels %s-%s)" % (PAGE1[0], PAGE1[-1])), ("fig2_suppb_en", " (page 2, panels %s-%s)" % (PAGE2[0], PAGE2[-1]))], "FigS2")}
B = Font(bold=True)
FILL = PATTERN = PatternFill("solid", fgColor="EAF2F8")
CODE_MD5 = md5(next(q for q in (os.path.join(HERE, "_fix", "srcdata.py"), os.path.join(HERE, "srcdata.py")) if os.path.exists(q)))[:10]
report = []
for fig, rs in PANELS.items():
    ftitle, fname, comp, spre = TITLE[fig]
    wb = Workbook(); ws = wb.active; ws.title = "README"
    ws.append(["Source data for %s" % ftitle]); ws["A1"].font = Font(bold=True, size=13)
    ws.append(["Each panel sheet lists the values drawn in that panel. The values were read directly from the plotted objects of the published panel "
               "(English version) immediately before it was saved, so they are exactly the numbers shown in the figure."])
    ws.append(["Series types: line (x, y); scatter points (x, y); bars / rectangles (position, height); filled band (x, lower and upper edge); "
               "line segments (e.g. error bars: start and end); reference line (horizontal or vertical guide; x or y in axes fraction 0 to 1); "
               "image (matrix); text in panel (value labels, cell numbers, annotations). Names come from the panel legend; an empty name means the object has no legend entry."])
    ws.append(["Panel letters follow the current figure. Composite: %s; extraction code srcdata.py md5 %s; generated %s."
               % (" and ".join("%s.png md5 %s%s" % (c, md5(os.path.join(P, "merged_Fig2_2026-10-10", c + ".png"))[:12], pg) for c, pg in comp),
                  CODE_MD5, datetime.date.today().isoformat())])
    ws.append([])
    hdr = ["Panel", "Sheet", "Panel content", "Denominator / cohort", "Unit of aggregation", "Definition (measure)", "Cut / threshold", "Orientation",
           "x-axis label", "y-axis label", "Legend entries", "Series blocks", "Data rows", "Panel image (English)", "Panel image md5 (png)"]
    ws.append(hdr)
    for c in ws[ws.max_row]: c.font = B; c.fill = FILL
    for r in rs:
        letter = r[H["칸"]].upper()
        cell = r[H["본문칸그림"]].replace(".svg", "")
        sheet = spre + letter
        J = json.load(open(os.path.join(SRC, cell + "__en.json"), encoding="utf8"))
        sh = wb.create_sheet(sheet)
        sh.append(["%s, panel %s" % (ftitle, letter)]); sh["A1"].font = Font(bold=True, size=12)
        sh.append(["Panel content", en(r[H["패널제목"]])])
        sh.append(["Denominator / cohort", en(r[H["분모"]])])
        sh.append(["Unit of aggregation", en(r[H["집계단위"]])])
        sh.append(["Panel image (English)", cell + ".png"])
        sh.append([])
        nblk = nrow = 0
        xl, yl, legs = [], [], []
        for a in J["axes"]:
            sub = (" (axis %d)" % a["axis_index"]) if len(J["axes"]) > 1 else ""
            if a["x_label"]: xl.append(a["x_label"])
            if a["y_label"]: yl.append(a["y_label"])
            legs += [t for _, t in a["legend"]]
            sh.append(["Axis%s" % sub, "x-axis: %s" % (a["x_label"] or "(none)"), "y-axis: %s" % (a["y_label"] or "(none)"),
                       "x range: %s to %s" % tuple(a["x_lim"]), "y range: %s to %s" % tuple(a["y_lim"])])
            sh.cell(sh.max_row, 1).font = B
            if a["title"]: sh.append(["Panel title", a["title"]])
            cat = [t for t in a["x_ticks"] if t[1] and not re.fullmatch(r"[−\-\d.,%]+", t[1])]
            if cat: sh.append(["x tick labels (position: label)"] + ["%s: %s" % (p, t) for p, t in cat])
            ycat = [t for t in a["y_ticks"] if t[1] and not re.fullmatch(r"[−\-\d.,%]+", t[1])]
            if ycat: sh.append(["y tick labels (position: label)"] + ["%s: %s" % (p, t) for p, t in ycat])
            # lines with the same kind, name, color, marker, style and column are grouped (like supplement K, with hundreds of lines drawn per haplotype)
            grp = collections.OrderedDict()
            for s in a["series"]:
                key = (s["type"], s.get("name", ""), s.get("color"), s.get("marker", ""), s.get("linestyle", ""), tuple(s["columns"]))
                grp.setdefault(key, []).append(s)
            for (typ, nm, col, mk, ls, cols), ss in grp.items():
                sh.append([])
                many = len(ss) > 1
                label = "Series: %s | type: %s%s%s" % (nm or "(no legend entry)", typ, (" | colour %s" % col) if col else "", (" | %d objects" % len(ss)) if many else "")
                sh.append([label]); sh.cell(sh.max_row, 1).font = B
                if typ.startswith("image") and ss[0].get("extent"): sh.append(["extent (left, right, bottom, top)"] + ss[0]["extent"])
                sh.append((["object"] if many else []) + list(cols))
                for c in sh[sh.max_row]: c.font = B; c.fill = FILL
                for oi, s in enumerate(ss, 1):
                    for row in s["rows"]:
                        sh.append(([oi] if many else []) + list(row)); nrow += 1
                nblk += 1
        for ci in range(1, min(sh.max_column, 12) + 1):
            sh.column_dimensions[get_column_letter(ci)].width = 16
        sh.column_dimensions["A"].width = 30
        ws.append([letter, sheet, en(r[H["패널제목"]]), en(r[H["분모"]]), en(r[H["집계단위"]]), en(r[H["자"]]), en(r[H["컷"]]), en(r[H["방향"]]),
                   " | ".join(dict.fromkeys(xl)), " | ".join(dict.fromkeys(yl)), " | ".join(dict.fromkeys(legs)), nblk, nrow, cell + ".png",
                   md5(os.path.join(P, "figures_en", cell + ".png"))[:12]])
        report.append((fig, letter, cell, nblk, nrow))
    for ci, w in enumerate([7, 10, 48, 48, 40, 30, 30, 22, 30, 30, 40, 8, 9, 26, 16], 1):
        ws.column_dimensions[get_column_letter(ci)].width = w
    for row in ws.iter_rows(min_row=6):
        for c in row: c.alignment = Alignment(wrap_text=True, vertical="top")
    # American spelling in every text cell (the manuscript uses American English)
    _US = [(r"\bcolour", "color"), (r"\bColour", "Color"), (r"\bneighbour", "neighbor"), (r"\bNeighbour", "Neighbor"),
           (r"\bgrey\b", "gray"), (r"\bGrey\b", "Gray"), (r"\bcentre\b", "center"), (r"\bCentre\b", "Center"),
           (r"\blabell(ed|ing)\b", r"label\1"), (r"\bnormalis", "normaliz"), (r"\bhomogenis", "homogeniz"), (r"\bbehaviour", "behavior"), (r"\btumour", "tumor")]
    for _sh in wb.worksheets:
        for _row in _sh.iter_rows():
            for _c in _row:
                if isinstance(_c.value, str):
                    _v = _c.value
                    for _p, _r in _US: _v = re.sub(_p, _r, _v)
                    if _v != _c.value: _c.value = _v
    out = os.path.join(OUTD, fname)
    wb.save(out)
    # check that no Korean is left (all sheets, all cells)
    from openpyxl import load_workbook
    wb2 = load_workbook(out, read_only=True)
    ko = [(s.title, c.coordinate if hasattr(c, "coordinate") else "?", c.value) for s in wb2.worksheets for row in s.iter_rows() for c in row
          if isinstance(c.value, str) and HAN.search(c.value)]
    print("■ %s · 시트 %d · 한글 칸 %d · md5 %s" % (fname, len(wb2.worksheets), len(ko), md5(out)[:10]))
    assert not ko, ko[:5]
for fig, letter, cell, nb, nr in report:
    print("   %-5s %s %-24s 묶음 %3d · 줄 %6d" % (fig, letter, cell, nb, nr))
