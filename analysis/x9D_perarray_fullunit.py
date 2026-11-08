# -*- coding: utf-8 -*-
# 완전유닛 자료를 array 단위로 뽑는다
#   v13 파이프라인(x9A)과 같은 여과·같은 판정 이지만, 덤프에 array_id·유닛번호·시작좌표를 싣고
#   자리별 프로파일을 array 단위로도 누적한다 → 대주제 A·B·C·E 를 완전유닛 기준으로 다시 그릴 수 있다.
#   (지금까지 A·B·C·E 는 느슨한 여과 n_calls >= 0.5 x array 중앙 을 썼다 — 문서 25 의 기준 불일치 ①)
#   산출 : tables/UNITBETA_PERARRAY_v13.tsv.gz · tables/PROFILE_BP_PERARRAY_v13.tsv.gz
#          (PROFILE_*_v13X 는 x9A 결과와 같아야 하는 대조용 이다 — 정본은 x9A 산출)
# 아래는 v13 설명 그대로 :
# 유닛 안 메틸 프로파일 v13 — 판정을 문턱에서 검정으로
#   v12 까지는 Ashman D·Sarle BC·소수성분 문턱 + 귀무 관문으로 '이봉' 을 정했다. 그 문턱은 유닛 수에 흔들렸다
#   (수마트라 과분산 귀무 통과율 전수 3.3% ↔ 3,000 부분표본 13.7%).
#   v13 은 묶음마다 베타-이항 1 성분 vs 2 성분 혼합 을 적합해 셋을 모두 요구한다 :
#     ① 모수 부트스트랩 우도비가 유의 (p <= 0.01 · 귀무 LR 분포는 적합된 1 성분에서 모의)
#     ② split-half 재현 (성분 평균 차 <= 간격/4 · 저몫 차 <= 0.05)
#     ③ 적합이 만드는 분자 수준 밀도의 봉우리 >= 2 개
#   셋 중 ①만 서면 보류 로 적고 분할은 단봉과 같게(가장 가까운 컷) 한다. 옛 판정은 '옛판정_v12' 열에 남긴다.
#   아래는 v12 설명 그대로 :
# 유닛 안 메틸 프로파일 v12 — v11 + 가장 가까운 컷 + 과분산 귀무
# A unimodal group is split with nearest_cut rather than with a hard cut on the ratio (ratio_split).
#              하드 컷은 컷 값 동점이 전부 HIGH 로 가 달성 LOW 가 목표 0.2191 보다 늘 작았다(중앙 0.175·최대 차 0.154).
#   바뀐 점 ② 과분산(베타-이항) 귀무 통과율 을 함께 적는다 — '두 상태' 가 아니라 '하나의 넓은 분포' 여도
#              3 조건이 통과하는지 보는 더 센 귀무. 판정에는 쓰지 않고 표에만 적는다(판정 규칙은 v11 그대로).
#   아래는 v11 설명 그대로 :
#   v10 에서 제브라 25/31 묶음이 '이봉' 으로 나왔는데, 실측해 보니 CpG 7~8 개짜리 묶음은
#   진짜 단봉이어도 3 조건이 99~100% 로 '이봉' 을 뱉는다 (표 BIMODAL_NULL_CHECK.tsv).
#   그래서 묶음마다 같은 n_cpg·같은 평균 β 의 단봉(이항) 모의 를 돌려 3 조건 통과율을 재고,
#   통과율 >= 5% 면 이봉으로 인정하지 않고 단봉으로 보아 비율 분할(f) 을 쓴다. 통과율을 표에 적는다.
#   (CpG 가 많은 유인원·금화조·애기장대는 통과율 0.0% 로 영향 없음)""
#   ① 묶음(array)마다 따로 이봉 판정 : Ashman D > 2 (가우스 2성분 EM) · Sarle BC > 0.555 · 소수 성분 >= 5%
# Definition check: human 5S (P10) gives D 3.002, BC 0.8311, w_minor 0.306, matching the specification values 3.00 / 0.83 / 0.306
#   ② 이봉 → 그 묶음 자신의 두 이항 성분(EM) 우도비 로 가른다 (p_low·p_high·동치 컷을 표에 적음)
#   ③ 단봉 → 사람 5S 의 LOW 비율 f 를 그 묶음 순위에 적용 (상보 분할 · 가운데 비우지 않음)
#      f 는 우리 사람 자료풀(P10 25,644 완전유닛) 에서 직접 냈다 : p_low 0.1580 · p_high 0.7867 · 컷 0.4610 · f = 0.2191
# (matches the specification values to the decimal; table HUMAN_RATIO_F.tsv)
#   ④ 금화조·제브라의 '고메틸 array vs 저메틸 array' 는 array 사이 차이다 — 합쳐서 이봉이라 부르지 않는다(묶음마다 따로 판정)
#   그 여과는 v5 에서 '흔한 자리 90% 덮음' 을 쓰려고 넣었던 것인데, 완전유닛 기준이 참조 CpG 90% 관측 으로 바뀐 뒤로는
#   같은 일을 두 번 하는 셈이고 유닛 길이가 들쭉날쭉한 종(침팬지)에서만 과하게 깎였다(완전유닛 679 → 여과 없이 4,666).
#   그래서 완전유닛 = n_cpg >= MINC 이고 그 유닛의 참조 CpG 중 90% 이상을 관측한 read-유닛, 하나로 간다.
#   · 묶음마다 x10 원 함수(density/call) 를 그 묶음의 완전유닛 β 에 돌려 이봉/단봉과 봉우리 둘을 얻는다
#   · 이봉 이면 우도비로 가른다 : d = n_meth*ln(a1/a0) + (n_cpg-n_meth)*ln((1-a1)/(1-a0)) · d>0 → hi
# - if unimodal, the top and bottom 25% by rank are used; the grouping column of the table records which was used
#   · 그 밖(하드콜 ML>=128 · 완전유닛=참조 CpG 90% 관측 · bp 자리 · 깊이가중)은 v7 과 같다
#   v6 은 값 문턱(q1c/q4c)을 썼는데, 유닛의 75% 넘게 β=1.00 인 묶음(제브라 177 bp · 금화조 chr2/chr4)에서
#   하위컷=상위컷=1.0 이 되어 전부 q1 으로 들어가고 q4 가 비었다(그림에서 하위 25% n = 전체 n).
#   → 정렬 순위로 아래 25% / 위 25% 를 정확히 잘라낸다(동점은 씨앗 0 의 무작위 순서).
# Why v4/v5 came out 0 (measured): the "90% coverage of common positions" criterion holds only in a canonical frame where positions are shared between units, as in the human case.
#     타종은 유닛마다 서열이 조금씩 달라 같은 CpG 가 다른 bp 자리에 오고, '흔한 자리' 가 유닛당 CpG 의 2~8 배가 된다 → 아무도 90% 를 못 덮는다.
#   그래서 같은 뜻(부분 read 걸러내기)을 그 유닛의 참조 CpG 수 로 바꾼다 :
#     완전유닛 = 그 유닛(참조 서열)이 가진 CpG 의 90% 이상 을 관측한 read-유닛 · 그리고 n_cpg >= MINC
#   ※ v4(x60)는 완전유닛 여과에서 모든 유닛이 탈락했다. 원인 : 사람 판의 자리 c 는 정본틀 이라 유닛끼리 자리가 공유되는데
#     우리 좌표는 유닛 길이가 제각각(고릴라 2,248/2,246/2,244/2,183…)이라 같은 CpG 가 여러 bp 자리로 번져
#     '흔한 자리' 가 유닛당 CpG 의 2~8 배가 됐다 → 어떤 유닛도 90% 를 못 덮는다.
#   그래서 먼저 g8 지도로 틀이 맞는 유닛 만 고른다 : unit_len 이 그 array 의 최빈값 ±2 이고 cov_start<=2 · cov_end>=unit_len-2
#   ① 하드콜 : CpG 마다 ML >= 128 → 0/1 (확률평균 안 씀)
#   ② 집계 단위 = bp 자리 그대로 (칸으로 미리 묶지 않는다)
#   ③ β = 깊이가중 : 자리마다 메틸수 합 / 호출수 합 (유닛별 평균의 평균 금지)
#   ④ 자리마다 최소 깊이 관문 (그릴 때 ±HALF bp 창 안 호출 >= MINCALL)  ← 그리는 쪽 x61 에서
#   ⑤ 미끄럼 평균은 그릴 때만 (bp 로 고른 깊이가중 창)                ← x61
#   완전유닛 여과 (이걸 빼면 가짜 저메틸 봉우리) :
#       (가) n_cpg >= MINC  · (나) '흔한 자리'(자리별 관측수 >= 0.10 x 최대 자리 관측수)의 90% 이상을 덮는 유닛만
#       ※ 타종은 유닛당 CpG 가 7~190 으로 종마다 달라 MINC 를 그 묶음 중앙 CpG 의 절반(최대 20)으로 둔다 — 표에 적는다.
#   무리 : all · q1(하위 25%) · q4(상위 25%)  ← 순위 절단(자료가 찾은 두 무리가 아님)
#   좌표 틀 : g8 유닛틀(유전자가 유닛 끝) · array strand 로 뒤집음 · 유전자 구간을 같은 틀로 옮겨 적는다
#   산출 : tables/PROFILE_BP_v13X.tsv.gz (묶음 x 무리 x bp 자리 : 호출수·메틸수·beta)
#          tables/PROFILE_META_v13X.tsv  (완전유닛 수·버린 수·흔한자리·컷·Ashman D·유전자 구간)
import io, os, sys, glob, re, struct, collections, math
from scipy import stats
import numpy as np
import pandas as pd

import os as _os_need, sys as _sys_need


def _need(k):
    """No defaults are set: if an environment variable is missing the script stops instead of silently reading the wrong location."""
    v = _os_need.environ.get(k)
    if not v:
        _sys_need.exit("environment variable %s is not set (the deposited copy has no default; "
                       "see README.md, section 'How to run it')" % k)
    return v


NEW = _need("RDNA_OUT_DIR")
B = _need("RDNA_DATA_ROOT")
OLD = _need("RDNA_PREV_OUT_DIR")
sys.path.insert(0, NEW + "/code")
import arch  # noqa
import bbmix  # 베타-이항 1/2 성분 가중 우도 적합(원판 xa1 과 값 동일 — 표 BBMIX_CHECK.tsv)
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
LOG = io.open(NEW + "/logs/X9D_PERARRAY.txt", "w", encoding="utf8", newline="\n")
REASON = "완전유닛 array 단위 덤프(A·B·C·E 완전유닛판 자료)(2026-10-04)"
MLTHR, COVMIN = 128, 0.9
FRATIO = 0.2191
NSIM, NCAP = 100, 3000   # 옛 귀무 모의(기록용) : 묶음당 모의 횟수 · 유닛 표본 상한(씨앗 0)
# ---- v13 판정 규칙 : 문턱(D·BC)은 기록만 하고, 판정은 아래 셋을 모두 요구한다
BOOT = 99          # 부트스트랩 벌 수 → p 의 바닥은 1/(BOOT+1) = 0.01
LRMIN = 30.0       # 이보다 작은 LR 은 부트스트랩 없이 '비유의' (실측된 귀무 LR 최대가 15.5 였다)
PALPHA = 0.01      # 유의 문턱 (= p 바닥)
WTOL = 0.05        # split-half 에서 저몫 차 허용
SEPDIV = 4.0       # split-half 에서 성분 평균 차 허용 = (두 성분 간격)/SEPDIV
# ---- 파일럿(xa5·xa6)에서 드러난 구멍을 막는 관문 — 문턱은 실측 분포에서 골랐다(표 V13_RULE_V2.tsv)
NULLMAX_OK = 100.0  # 귀무 LR 최대가 이보다 크면 적합이 폭주한 것 → 검정 불가(실측 6e8 까지 나왔다)
SEPRES = 2.0        # 두 성분 간격 >= SEPRES x (1/CpG중앙) — 양자화가 만든 성분 제외
MINW = 0.05 # lower bound on the minor-component weight, so that a 1.7% tail is not counted as a peak
VALLEY = 0.6        # 골깊이비 상한 — 0.9 를 넘는 '어깨' 를 두 봉우리로 세지 않기 위해
NSPLIT, STARTS = 5, 5   # split-half 를 5 벌 (한 벌이면 최적화 국소해 때문에 판정이 뒤집힌다 — 실측)


def em_two_binomial(kn, iters=200): # two-component binomial EM
    p0, p1, w = 0.20, 0.80, 0.40
    for _ in range(iters):
        a0 = b0 = a1 = b1 = sw = 0.0
        for k, n in kn:
            l0 = k * math.log(p0) + (n - k) * math.log(1 - p0) + math.log(max(w, 1e-12))
            l1 = k * math.log(p1) + (n - k) * math.log(1 - p1) + math.log(max(1 - w, 1e-12))
            m = max(l0, l1)
            e0, e1 = math.exp(l0 - m), math.exp(l1 - m)
            r0 = e0 / (e0 + e1)
            sw += r0; a0 += r0 * k; b0 += r0 * n; a1 += (1 - r0) * k; b1 += (1 - r0) * n
        w = sw / len(kn)
        p0 = a0 / b0 if b0 else p0
        p1 = a1 / b1 if b1 else p1
        if p0 > p1:
            p0, p1, w = p1, p0, 1 - w
    return p0, p1, w


def lr_split(kn, p0, p1):
    cut = (math.log((1 - p0) / (1 - p1)) / (math.log(p1 / p0) + math.log((1 - p0) / (1 - p1))))
    low = {i for i, (k, n) in enumerate(kn)
           if k * math.log(p1 / p0) + (n - k) * math.log((1 - p1) / (1 - p0)) <= 0}
    return low, cut


def ratio_split(betas, f):
    """Hard-cut split on the ratio. Kept for comparison; nearest_cut below is used instead."""
    b = sorted(betas)
    i = max(1, min(len(b) - 1, int(round(len(b) * f))))
    return b[i]


def nearest_cut(betas, f):
    """v12 : 달성 비율이 목표 f 에 가장 가까운 컷
       A hard cut on the ratio sends every tie at the cut value to HIGH, so the attained LOW fraction was always below target
       (v11 실측 : 달성 중앙 0.175 · 목표 0.2191 · |차| 최대 0.154).
       서로 다른 β 값을 전부 컷 후보로 두고 (β < c) 의 몫이 f 에 가장 가까운 c 를 고른다.
       동점은 쪼개지 않는다(값이 같은 유닛은 늘 같은 무리). 동률이면 작은 컷을 골라 재현되게 한다."""
    b = np.asarray(betas, dtype=float)
    uq = np.unique(b)
    fr = np.array([float((b < c).mean()) for c in uq])
    j = int(np.argmin(np.abs(fr - f)))
    return float(uq[j]), float(fr[j])


def betabinom_null_rate(kk, nn, nsim, rng, iters=60):
    """과분산 귀무 — '두 상태' 가 아니라 하나의 넓은(연속) 분포 여도 3 조건을 통과하는가?
       단봉 이항 귀무는 표본오차만 담아 분산이 너무 작다(생물 변이가 0). 그래서 같은 (k,n) 에
       적률법으로 베타-이항 을 맞춰 모의한다. 반환 : (통과율, α, β, 쓸 수 있나)
       α<1 또는 β<1 이면 베타 자체가 U 자(이봉)라 '단봉 귀무' 가 아니다 → 쓸 수 없다고 표시한다."""
    p = kk / nn
    m = float(kk.sum() / nn.sum())
    s2 = float(np.var(p, ddof=1)) if len(p) > 1 else 0.0
    inv = float(np.mean(1.0 / nn))
    denom = m * (1 - m) * (1 - inv)
    rho = (s2 - m * (1 - m) * inv) / denom if denom > 1e-12 else 0.0
    rho = float(min(max(rho, 1e-6), 0.999))
    t = 1.0 / rho - 1.0
    a, bta = m * t, (1 - m) * t
    usable = bool((a >= 1.0) and (bta >= 1.0))
    hits = 0
    for _ in range(nsim):
        pp = rng.beta(a, bta, size=len(nn))
        sim = rng.binomial(nn.astype(int), pp) / nn
        if len(sim) > 4 and sim.std() > 1e-9:
            g1s = float(stats.skew(sim, bias=False)); kus = float(stats.kurtosis(sim, bias=False))
            ns = len(sim)
            BCs = (g1s ** 2 + 1) / (kus + 3 * ((ns - 1) ** 2) / ((ns - 2) * (ns - 3)))
            ms = np.array([np.quantile(sim, .2), np.quantile(sim, .8)]); sss = np.array([sim.std(), sim.std()]) + 1e-6
            wws = np.array([.5, .5])
            for _i in range(iters):
                pdfs = wws / (sss * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((sim[:, None] - ms) / sss) ** 2)
                rrs = pdfs / np.clip(pdfs.sum(1, keepdims=True), 1e-300, None)
                nks = rrs.sum(0) + 1e-12
                ms = (rrs * sim[:, None]).sum(0) / nks
                sss = np.sqrt((rrs * (sim[:, None] - ms) ** 2).sum(0) / nks) + 1e-6
                wws = nks / ns
            ADs = float(np.sqrt(2) * abs(ms[0] - ms[1]) / np.sqrt(sss[0] ** 2 + sss[1] ** 2))
            hits += (ADs > 2 and BCs > 0.555 and min(wws) >= 0.05)
    return hits / float(nsim), float(a), float(bta), usable


def P(s=""):
    print(s)
    LOG.write(s + "\n")
    LOG.flush()


src = io.open(OLD + "/code/x11_coding_readunit.py", encoding="utf8").read()
L11 = 'LOG = io.open(OUT + ("/X11_CODING_PILOT.txt" if PILOT else "/X11_CODING_BUILD.txt"), "a", encoding="utf8", newline="\\n")'
head = src.split("STAT = collections.Counter()")[0].replace(L11, "LOG = io.StringIO()")
G = {"__name__": "x11hdr"}
exec(compile(head, OLD + "/code/x11_coding_readunit.py", "exec"), G)
NIBCOMP, load_species, bam_records, parse_tags, read_fa = (
    G["NIBCOMP"], G["load_species"], G["bam_records"], G["parse_tags"], G["read_fa"])
# ---- x10 의 판정 함수(density/call)를 원문 그대로 가져온다 (내가 옮겨 적지 않는다)
_x10 = io.open(NEW + "/code/x10_readunit_per_array.py", encoding="utf8").read()
_hdr = _x10.split("def density(")[0]
_body = "def density(" + _x10.split("def density(")[1].split('P("=" * 104)')[0]
_ns = {"__name__": "x10fn"}
_pre = "import numpy as np\n" + "\n".join([l for l in _hdr.split("\n") if re.match(r"^SEP, TROUGH, MINOR *=", l)])
exec(compile(_pre + "\n" + _body, "x10_fn", "exec"), _ns)
x10_call = _ns["call"]

BC = pd.read_csv(NEW + "/tables/BONOBO_READ_BASECALLER.tsv", sep="\t")
B637 = set(BC.loc[BC.basecaller.astype(str).str.contains("6.3.7"), "read_id"])
MT3 = pd.read_csv(NEW + "/tables/PERARRAY_META_v3.tsv", sep="\t")      # 묶음 정의(= x52 가 쓴 것)를 그대로 쓴다
GRP_OF = {}
for t in MT3.itertuples(index=False):
    for a in str(t.array목록).split("|"):
        GRP_OF[a] = t.묶음
GSPEC = {t.묶음: dict(종=t.종, 유닛길이=int(t.유닛길이)) for t in MT3.itertuples(index=False)}
SPK = [("고릴라", "gorilla_mGorGor1"), ("보노보", "bonobo_mPanPan1"), ("침팬지", "chimp_mPanTro3"),
       ("보르네오오랑우탄", "orang_bornean_mPonPyg2"), ("수마트라오랑우탄", "orang_sumatran_mPonAbe1"),
       ("시아망", "siamang_mSymSyn1"), ("애기장대", "arabidopsis_ColCC"),
       ("금화조", "zfinch_v03_dip"), ("제브라피시", "zebrafish_GRCz12tu")]


def ashman(x):
    """2성분 정규 혼합 EM 으로 Ashman D = sqrt(2)|m1-m2| / sqrt(s1^2+s2^2)"""
    x = np.asarray(x, float)
    if len(x) < 50 or x.std() < 1e-6:
        return np.nan
    m = np.array([np.quantile(x, 0.2), np.quantile(x, 0.8)])
    s = np.array([x.std(), x.std()]) + 1e-6
    w = np.array([0.5, 0.5])
    for _ in range(200):
        pdf = w[None, :] / (s[None, :] * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((x[:, None] - m[None, :]) / s[None, :]) ** 2)
        tot = pdf.sum(axis=1, keepdims=True)
        tot[tot == 0] = 1e-300
        r = pdf / tot
        nk = r.sum(axis=0) + 1e-12
        m = (r * x[:, None]).sum(axis=0) / nk
        s = np.sqrt((r * (x[:, None] - m[None, :]) ** 2).sum(axis=0) / nk) + 1e-6
        w = nk / len(x)
    return float(np.sqrt(2) * abs(m[0] - m[1]) / np.sqrt(s[0] ** 2 + s[1] ** 2))


ROWS, META, UROWS, AROWS = [], [], [], [] # UROWS = per unit (CpG count, methylated count, beta, group), for the "cut line on the beta distribution" panel
for ko, sp in SPK:
    Asp = load_species(sp)
    STR, FRAMEOK, FRAMEN = {}, {}, {}
    for f_ in sorted(glob.glob(B + "/02_consensus_map/%s/*.g8map.tsv" % sp)):
        if not re.search(r"_[0-9]+_[0-9]+[.]g8map[.]tsv$", f_.replace("\\", "/")):
            continue
        g_ = pd.read_csv(f_, sep="\t")
        aid_ = g_.array_id.iloc[0]
        STR[aid_] = g_.strand.iloc[0]
        # 틀이 맞는 유닛만 : 길이가 그 array 최빈값 ±2 이고 지도 덮임이 유닛 전체인 것
        modal = int(g_.unit_len.mode().iloc[0])
        ok_ = g_[(g_.unit_len.sub(modal).abs() <= 2) & (g_.cov_start <= 2) & (g_.cov_end >= g_.unit_len - 2)]
        FRAMEOK[aid_] = set(int(x) for x in ok_.array_pos)
        FRAMEN[aid_] = (len(ok_), len(g_), modal)
    REFCPG = {}                      # (array, 유닛index) -> 그 유닛 참조 서열의 CpG 수
    for a in [x for x in Asp if x in GRP_OF]:
        v = Asp[a]; sq = v["seq"]
        cg = np.flatnonzero((sq[:-1] == 67) & (sq[1:] == 71)) + 1 + (v["lo1"] - 1)   # 1기반 게놈 좌표
        ui_ = np.searchsorted(v["us"], cg, side="right") - 1
        ok2 = (ui_ >= 0) & (cg <= v["ue"][np.clip(ui_, 0, None)])
        cnt_ = collections.Counter(ui_[ok2].tolist())
        REFCPG[a] = cnt_
    keep_aids = [a for a in Asp if a in GRP_OF]
    bychrom = collections.defaultdict(list)
    for a in keep_aids:
        bychrom[Asp[a]["chrom"]].append(a)
    # 묶음마다 read-유닛 모으기 : (n_cpg, n_meth, 자리 비트, 메틸 비트)
    UNITS = collections.defaultdict(list)
    for bam in sorted(glob.glob(B + "/04_bam/%s/*.5S.bam" % sp)):
        for refs, r in bam_records(bam):
            refid, pos0, lrn, mapq, _bin, ncig, flag, lseq = struct.unpack("<iiBBHHHi", r[:20])
            if refid < 0 or (flag & 0x904) or refs[refid] not in bychrom:
                continue
            rid = r[32:32 + lrn - 1].decode()
            if sp == "bonobo_mPanPan1" and rid in B637:
                continue
            p = 32 + lrn
            cig = np.frombuffer(r[p:p + 4 * ncig], dtype="<u4"); p += 4 * ncig
            sq = r[p:p + (lseq + 1) // 2]; p += (lseq + 1) // 2 + lseq
            T = parse_tags(r, p)
            if "MM" not in T or "ML" not in T:
                continue
            ops, lens = cig & 15, (cig >> 4).astype(np.int64)
            if (ops == 5).any():
                continue
            q2r = np.full(lseq, -1, dtype=np.int64)
            q, rr = 0, pos0
            for o, n in zip(ops.tolist(), lens.tolist()):
                if o in (0, 7, 8):
                    q2r[q:q + n] = np.arange(rr, rr + n); q += n; rr += n
                elif o in (1, 4):
                    q += n
                elif o in (2, 3):
                    rr += n
            bb = np.frombuffer(sq, dtype=np.uint8)
            nib = np.empty(2 * len(bb), dtype=np.uint8); nib[0::2] = bb >> 4; nib[1::2] = bb & 15
            codes = nib[:lseq]
            rev = bool(flag & 16)
            ocodes = NIBCOMP[codes][::-1] if rev else codes
            cpos = np.flatnonzero(ocodes == 2)
            ml = T["ML"]; mi = 0
            m_idx = m_ml = None
            for ent in T["MM"].rstrip(";").split(";"):
                if not ent:
                    continue
                parts = ent.split(","); hd2 = parts[0]; mode = "?"
                if hd2[-1] in "?.":
                    mode = hd2[-1]; hd2 = hd2[:-1]
                base, st_, mc = hd2[0], hd2[1], hd2[2:]
                d = np.array([int(x) for x in parts[1:]], dtype=np.int64)
                nc = len(mc) if not mc.isdigit() else 1
                vals = ml[mi:mi + nc * len(d)].reshape(-1, nc) if len(d) else np.zeros((0, nc), np.uint8)
                mi += nc * len(d)
                if base != "C" or st_ != "+" or "m" not in mc or len(d) == 0:
                    continue
                k = np.cumsum(d + 1) - 1
                if k[-1] >= len(cpos):
                    break
                m_idx, m_ml = cpos[k], vals[:, mc.index("m")]
                if mode == ".":
                    full = np.zeros(len(cpos), dtype=np.uint8); full[k] = m_ml
                    m_idx, m_ml = cpos, full
            if m_idx is None:
                continue
            sidx = (lseq - 1 - m_idx) if rev else m_idx
            ref0 = q2r[sidx]
            ok = ref0 >= 0
            ref0, mlv = ref0[ok], m_ml[ok]
            for aid in bychrom[refs[refid]]:
                v, s_ = Asp[aid], Asp[aid]["seq"]
                key = GRP_OF[aid]
                L = GSPEC[key]["유닛길이"]
                k0 = ref0 - (v["lo1"] - 1)
                inr = (k0 >= 1) & (k0 < len(s_) - 1)
                if not inr.any():
                    continue
                kk, vv, rp = k0[inr], mlv[inr], ref0[inr] + 1
                cg = (s_[kk] == 71) & (s_[kk - 1] == 67) if rev else (s_[kk] == 67) & (s_[kk + 1] == 71)
                vv, rp = vv[cg], rp[cg]
                if not len(rp):
                    continue
                call = (vv >= MLTHR).astype(np.uint8)          # ① 하드콜
                ui = np.searchsorted(v["us"], rp, side="right") - 1
                good = (ui >= 0) & (rp <= v["ue"][np.clip(ui, 0, None)])
                if not good.any():
                    continue
                ui2, call2, rp2 = ui[good], call[good], rp[good]
                inu = (rp2 - v["us"][ui2] + 1) if STR[aid] == "+" else (v["ue"][ui2] - rp2 + 1)
                for u in np.unique(ui2):
                    w = ui2 == u
                    pos = np.clip(inu[w].astype(np.int64), 1, L)
                    cl = call2[w]
                    bm_c = np.zeros(L + 2, dtype=bool); bm_m = np.zeros(L + 2, dtype=bool)
                    bm_c[pos] = True
                    bm_m[pos[cl.astype(bool)]] = True
                    UNITS[key].append((int(w.sum()), int(cl.sum()), np.packbits(bm_c), np.packbits(bm_m), int(REFCPG[aid].get(int(u), 0)), aid, int(u), int(v["us"][u])))
    # ---------------- 묶음마다 완전유닛 여과 → 자리별 누적
    for key, lst in UNITS.items():
        L = GSPEC[key]["유닛길이"]
        nbits = L + 2
        TOT = np.zeros(nbits, dtype=np.int64)
        for nc, nm, bc, bmm, nref, _a, _u, _st in lst:
            TOT += np.unpackbits(bc)[:nbits]
        mx = TOT.max() if TOT.size else 0
        CPGB = np.flatnonzero(TOT >= 0.10 * mx)                 # 흔한 자리 = 0.10 x 최대 자리 관측수
        med_cpg = float(np.median([x[0] for x in lst])) if lst else 0.0
        med_ref = float(np.median([x[4] for x in lst])) if lst else 0.0
        MINC = int(min(20, max(3, round(0.5 * med_cpg))))       # 종마다 CpG 수가 달라 적응(표에 적음)
        kept, keptinfo = [], []      # keptinfo = (array_id, 유닛번호, 시작좌표) — 완전유닛판 A·B·C·E 용
        drop_nc = drop_cov = 0
        for nc, nm, bc, bmm, nref, aid_, u_, st_ in lst:
            if nc < MINC:
                drop_nc += 1
                continue
            if nref > 0 and nc / float(nref) < COVMIN:        # 완전유닛 = 참조 CpG 의 90% 이상 관측
                drop_cov += 1
                continue
            bits = np.unpackbits(bc)[:nbits].astype(bool)
            kept.append((nc, nm, bits, np.unpackbits(bmm)[:nbits].astype(bool)))
            keptinfo.append((aid_, u_, st_))
        if not kept:
            P("  %s : 완전유닛 0 (버림 nc %d · cov %d) — 건너뜀" % (key, drop_nc, drop_cov))
            continue
        bvals = np.array([nm / float(nc) for nc, nm, _b, _m in kept])
        bb = np.sort(bvals)
        q1c, q4c = float(bb[len(bb) // 4]), float(bb[(3 * len(bb)) // 4])
        # ---- Bimodality rule, three conditions (Gaussian EM Ashman D, Sarle BC, minor component)
        nunit = len(bvals)
        g1 = float(stats.skew(bvals, bias=False)); ku = float(stats.kurtosis(bvals, bias=False))
        BC = (g1 ** 2 + 1) / (ku + 3 * ((nunit - 1) ** 2) / ((nunit - 2) * (nunit - 3))) if nunit > 4 else np.nan
        mm = np.array([np.quantile(bvals, .2), np.quantile(bvals, .8)])
        ss = np.array([bvals.std(), bvals.std()]) + 1e-6
        ww = np.array([.5, .5])
        for _ in range(300):
            pdf = ww / (ss * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((bvals[:, None] - mm) / ss) ** 2)
            rr = pdf / np.clip(pdf.sum(1, keepdims=True), 1e-300, None)
            nk = rr.sum(0) + 1e-12
            mm = (rr * bvals[:, None]).sum(0) / nk
            ss = np.sqrt((rr * (bvals[:, None] - mm) ** 2).sum(0) / nk) + 1e-6
            ww = nk / nunit
        AD = float(np.sqrt(2) * abs(mm[0] - mm[1]) / np.sqrt(ss[0] ** 2 + ss[1] ** 2))
        wmin = float(min(ww))
        cond3 = (AD > 2) and (BC > 0.555) and (wmin >= 0.05)
        # ---- 귀무 모의 : 같은 n_cpg·같은 평균으로 단봉(이항)을 만들어 3 조건 통과율
        kk_ = np.array([nm for _nc, nm, _b, _m in kept], dtype=float)
        nn_ = np.array([nc for nc, _nm, _b, _m in kept], dtype=float)
        rngs = np.random.default_rng(0)
        if len(nn_) > NCAP:
            sel_ = rngs.choice(len(nn_), NCAP, replace=False)
            kk_, nn_ = kk_[sel_], nn_[sel_]
        pbar = float(kk_.sum() / nn_.sum())
        hits = 0
        for _ in range(NSIM):
            sim = rngs.binomial(nn_.astype(int), pbar) / nn_
            if len(sim) > 4 and sim.std() > 1e-9:
                g1s = float(stats.skew(sim, bias=False)); kus = float(stats.kurtosis(sim, bias=False))
                ns = len(sim)
                BCs = (g1s ** 2 + 1) / (kus + 3 * ((ns - 1) ** 2) / ((ns - 2) * (ns - 3)))
                ms = np.array([np.quantile(sim, .2), np.quantile(sim, .8)]); sss = np.array([sim.std(), sim.std()]) + 1e-6
                wws = np.array([.5, .5])
                for _i in range(60):
                    pdfs = wws / (sss * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((sim[:, None] - ms) / sss) ** 2)
                    rrs = pdfs / np.clip(pdfs.sum(1, keepdims=True), 1e-300, None)
                    nks = rrs.sum(0) + 1e-12
                    ms = (rrs * sim[:, None]).sum(0) / nks
                    sss = np.sqrt((rrs * (sim[:, None] - ms) ** 2).sum(0) / nks) + 1e-6
                    wws = nks / ns
                ADs = float(np.sqrt(2) * abs(ms[0] - ms[1]) / np.sqrt(sss[0] ** 2 + sss[1] ** 2))
                hits += (ADs > 2 and BCs > 0.555 and min(wws) >= 0.05)
        nullrate = hits / float(NSIM)
        gate11 = cond3 and (nullrate < 0.05)       # v11·v12 의 옛 판정 — 기록만 한다(대조용)
        # ---- 과분산 귀무 (적는 값) : '하나의 넓은 분포' 여도 3 조건을 통과하나?
        bbrate, bba, bbb, bbok = betabinom_null_rate(kk_, nn_, NSIM, np.random.default_rng(1))
        # ================= v13 판정 : 문턱이 아니라 검정 =================
        #   ① 베타-이항 1 성분 vs 2 성분 혼합의 모수 부트스트랩 우도비 가 유의할 것 (p <= PALPHA)
        #   ② split-half 재현 : 무작위 반반으로 갈라 적합한 두 성분이 (간격/4) 안에서 맞고, 저몫 차 <= WTOL
        #   ③ 적합이 만드는 분자 수준 밀도의 봉우리가 2 개 이상 ('성분이 둘' != '봉우리가 둘')
        kk_all = np.array([nm for _nc, nm, _b, _m in kept], dtype=float)
        nn_all = np.array([nc for nc, _nm, _b, _m in kept], dtype=float)
        ku, nu, cu = bbmix.pack(kk_all, nn_all)
        LRobs, _lg1, _lg2, ab1, prm = bbmix.lr_stat(ku, nu, cu)
        if LRobs >= LRMIN:                          # 작은 LR 은 부트스트랩 없이 비유의(귀무 LR 최대 실측 15.5)
            rngb = np.random.default_rng(12345)
            nl = []
            for _ in range(BOOT):
                _ks, _ns, _cs = bbmix.boot_null(nn_all, ab1[0], ab1[1], rngb)
                nl.append(bbmix.lr_stat(_ks, _ns, _cs, starts=1)[0])
            nl = np.array(nl)
            pval = (np.sum(nl >= LRobs) + 1.0) / (BOOT + 1.0)
            nullmax = float(nl.max())
        else:
            pval, nullmax = 1.0, np.nan
        sep = float(prm["mu_high"] - prm["mu_low"])
        minw = float(min(prm["w_low"], 1 - prm["w_low"]))
        npk, vratio, pk0, pk1 = bbmix.mixture_peaks(prm)
        sepres = sep * med_cpg                       # 두 성분 간격이 분해능(1/CpG)의 몇 배인가
        # ---- 관문을 차례로 보고 탈락 사유별로 이름을 다르게 붙인다(뭉뚱그린 '보류' 금지)
        sig = bool(pval <= PALPHA)
        dlo = dhi = dw = np.nan
        npass = -1
        if not (np.isfinite(nullmax) and nullmax < NULLMAX_OK) and LRobs >= LRMIN:
            verdict = "검정불가(귀무 적합 폭주)"
        elif not sig:
            verdict = "단봉(2성분 유의하지 않음)"
        elif sepres < SEPRES:
            verdict = "판정불가(분해능 부족 · 간격이 1/CpG 의 %.2f 배)" % sepres
        elif npk < 2 or (np.isfinite(vratio) and vratio > VALLEY):
            verdict = "단봉(봉우리 하나%s)" % ("" if npk < 2 else " · 어깨 골깊이비 %.2f" % vratio)
        elif minw < MINW:
            verdict = "단봉+소수꼬리(소수 성분 %.3f)" % minw
        else:
            dlos, dhis, dws, npass = [], [], [], 0
            for _s in range(NSPLIT):               # 재현은 5 벌
                rh = np.random.default_rng(100 + _s)
                ix_ = rh.permutation(len(nn_all))
                hs = []
                for h_ in (ix_[:len(ix_) // 2], ix_[len(ix_) // 2:]):
                    _kh, _nh, _ch = bbmix.pack(kk_all[h_], nn_all[h_])
                    hs.append(bbmix.fit2(_kh, _nh, _ch, starts=STARTS)[1])
                a_ = abs(hs[0]["mu_low"] - hs[1]["mu_low"]); b_ = abs(hs[0]["mu_high"] - hs[1]["mu_high"])
                c_ = abs(hs[0]["w_low"] - hs[1]["w_low"])
                dlos.append(a_); dhis.append(b_); dws.append(c_)
                npass += int((a_ <= sep / SEPDIV) and (b_ <= sep / SEPDIV) and (c_ <= WTOL))
            dlo, dhi, dw = max(dlos), max(dhis), max(dws)
            verdict = "이봉" if npass == NSPLIT else "보류(재현 %d/%d)" % (npass, NSPLIT)
        is_bi = (verdict == "이봉")
        repro = (npass == NSPLIT)
        lab = np.full(nunit, "MID", dtype=object)
        if is_bi:
            kn = [(nm, nc) for nc, nm, _b, _m in kept]
            p0, p1, wlow = em_two_binomial(kn)
            lowset, cutv = lr_split(kn, p0, p1)
            for i in range(nunit):
                lab[i] = "q1" if i in lowset else "q4"
            method = ("우도비(이봉 · p %.4f/%.4f · 컷 %.4f · LRT p<=%.3f · LR %.0f)"
                      % (p0, p1, cutv, pval, LRobs))
        else:
            p0 = p1 = np.nan
            cutv, achv = nearest_cut(bvals, FRATIO)      # v12 : 달성이 목표에 가장 가까운 컷
            lab = np.where(bvals < cutv, "q1", "q4").astype(object)
            method = ("사람 5S 비율 목표 f=%.4f · 가장 가까운 컷 %.4f (달성 %.4f · %s)"
                      % (FRATIO, cutv, achv, verdict))
        acc = {g: np.zeros((2, nbits)) for g in ("all", "q1", "q4")}
        accA = collections.defaultdict(lambda: np.zeros((2, nbits)))   # array 단위 자리별 누적(완전유닛판 B·C·E 용)
        cnt = collections.Counter()
        cntA = collections.Counter()
        for (nc, nm, bits, mb), lb, (aid_, u_, st_) in zip(kept, lab, keptinfo):
            UROWS.append((GSPEC[key]["종"], key, aid_, int(u_), int(st_), int(nc), int(nm),
                          round(nm / float(nc), 6), lb))
            gs = ["all"] + ([lb] if lb in ("q1", "q4") else [])
            for g in gs:
                acc[g][0] += bits
                acc[g][1] += mb
                cnt[g] += 1
                accA[(aid_, g)][0] += bits
                accA[(aid_, g)][1] += mb
                cntA[(aid_, g)] += 1
        for g in ("all", "q1", "q4"):
            n, m = acc[g]
            for p_ in np.flatnonzero(n > 0):
                ROWS.append((GSPEC[key]["종"], key, g, L, int(p_), int(n[p_]), int(m[p_]), round(m[p_] / n[p_], 6)))
        for (aid_, g), arr in accA.items():                 # array 단위 자리별 프로파일
            n, m = arr
            for p_ in np.flatnonzero(n > 0):
                AROWS.append((GSPEC[key]["종"], key, aid_, g, L, int(p_), int(n[p_]), int(m[p_]),
                              round(m[p_] / n[p_], 6), int(cntA[(aid_, g)])))
        # 유전자 구간(같은 틀)
        aid0 = [a for a in GRP_OF if GRP_OF[a] == key][0]
        v0 = Asp[aid0] if aid0 in Asp else None
        gw = []
        if v0 is not None:
            for i in range(len(v0["us"])):
                for (x, y) in v0["genes"][i]:
                    gw.append((x - v0["us"][i] + 1, y - v0["us"][i] + 1) if STR[aid0] == "+"
                              else (v0["ue"][i] - y + 1, v0["ue"][i] - x + 1))
        fr = [FRAMEN[a] for a in GRP_OF if GRP_OF[a] == key and a in FRAMEN]
        META.append(dict(종=GSPEC[key]["종"], 묶음=key, 유닛길이=L, 원유닛=len(lst), 완전유닛=len(kept),
                         틀맞춤유닛=sum(x[0] for x in fr), 지도유닛=sum(x[1] for x in fr),
                         버림_CpG적음=drop_nc, 버림_덮임부족=drop_cov, MINC=MINC, 흔한자리=int(len(CPGB)),
                         최대자리관측=int(mx), CpG중앙=med_cpg, 참조CpG중앙=med_ref, q1컷=round(q1c, 4), q4컷=round(q4c, 4), 무리방법=method,
                         판정=verdict.split("(")[0].strip(), 판정상세=verdict,   # 판정 = 이봉/단봉/단봉+소수꼬리/보류/판정불가/검정불가
                         LRT_LR=round(float(LRobs), 2), LRT_귀무최대=(round(nullmax, 2) if np.isfinite(nullmax) else np.nan),
                         LRT_p=round(float(pval), 4), LRT_B=BOOT,
                         성분_저=round(prm["mu_low"], 4), 성분_고=round(prm["mu_high"], 4), 성분_저몫=round(prm["w_low"], 4),
                         소수성분몫=round(minw, 4), 간격배_분해능=round(float(sepres), 3),
                         반반_저차최대=(round(float(dlo), 4) if np.isfinite(dlo) else np.nan),
                         반반_고차최대=(round(float(dhi), 4) if np.isfinite(dhi) else np.nan),
                         반반_몫차최대=(round(float(dw), 4) if np.isfinite(dw) else np.nan),
                         반반허용=round(float(sep / SEPDIV), 4),
                         재현벌=("%d/%d" % (npass, NSPLIT) if npass >= 0 else "—"),
                         적합봉수=int(npk), 골깊이비=(round(float(vratio), 4) if np.isfinite(vratio) else np.nan),
                         한성분_알파=round(float(ab1[0]), 3), 한성분_베타=round(float(ab1[1]), 3),
                         옛판정_v12=("이봉" if gate11 else "단봉"), three조건=("통과" if cond3 else "탈락"), 귀무통과율=round(nullrate, 3),
                         과분산귀무통과율=round(bbrate, 3), 과분산귀무_알파=round(bba, 3), 과분산귀무_베타=round(bbb, 3),
                         과분산귀무_쓸수있나=("예" if bbok else "아니오(α·β<1 이라 귀무가 U자)"),
                         분할방식=("우도비" if is_bi else "가장 가까운 컷"),
                         달성LOW비율=round(float((lab == "q1").mean()), 4),
                         AshmanD=round(AD, 3), BC=round(BC, 4), 소수성분=round(wmin, 4),
                         컷=round(float(cutv), 4), p_low=(round(p0, 4) if np.isfinite(p0) else np.nan), p_high=(round(p1, 4) if np.isfinite(p1) else np.nan),
                         
                         n_all=cnt["all"], n_q1=cnt["q1"], n_q4=cnt["q4"],
                         유전자시작=(float(np.median([a for a, _b in gw])) if gw else np.nan),
                         유전자끝=(float(np.median([b for _a, b in gw])) if gw else np.nan)))
    P("%-10s 묶음 %2d · 완전유닛 %s (버림 %s)"
      % (ko, sum(1 for m in META if m["종"] == ko),
         format(sum(m["완전유닛"] for m in META if m["종"] == ko), ","),
         format(sum(m["버림_CpG적음"] + m["버림_덮임부족"] for m in META if m["종"] == ko), ",")))

R = pd.DataFrame(ROWS, columns=["종", "묶음", "무리", "유닛길이", "위치", "호출수", "메틸수", "beta"])
M = pd.DataFrame(META)
TMP = os.environ.get("TEMP", ".") + "/x60"
os.makedirs(TMP, exist_ok=True)
R.to_csv(TMP + "/PROFILE_BP_v13X.tsv.gz", sep="\t", index=False, encoding="utf8", lineterminator="\n")
M.to_csv(TMP + "/PROFILE_META_v13X.tsv", sep="\t", index=False, encoding="utf8", lineterminator="\n")
U = pd.DataFrame(UROWS, columns=["종", "묶음", "array_id", "유닛번호", "시작좌표", "n_cpg", "n_meth", "beta", "무리"])
AP = pd.DataFrame(AROWS, columns=["종", "묶음", "array_id", "무리", "유닛길이", "위치", "호출수", "메틸수", "beta", "유닛수"])
AP.to_csv(TMP + "/PROFILE_BP_PERARRAY_v13.tsv.gz", sep="\t", index=False, encoding="utf8", lineterminator="\n")
arch.put(TMP + "/PROFILE_BP_PERARRAY_v13.tsv.gz", "tables/PROFILE_BP_PERARRAY_v13.tsv.gz", REASON)
U.to_csv(TMP + "/UNITBETA_PERARRAY_v13.tsv.gz", sep="\t", index=False, encoding="utf8", lineterminator="\n")
arch.put(TMP + "/PROFILE_BP_v13X.tsv.gz", "tables/PROFILE_BP_v13X.tsv.gz", REASON)
arch.put(TMP + "/PROFILE_META_v13X.tsv", "tables/PROFILE_META_v13X.tsv", REASON)
arch.put(TMP + "/UNITBETA_PERARRAY_v13.tsv.gz", "tables/UNITBETA_PERARRAY_v13.tsv.gz", REASON)
P("")
P("묶음 %d · 자리 행 %s" % (len(M), format(len(R), ",")))
P("")
P("유인원 — 유전자부 vs 스페이서 (완전유닛 · 깊이가중)")
if not len(M):
    P("  완전유닛이 하나도 없다 — 여과를 다시 볼 것"); LOG.close(); sys.exit(0)
for t in M[M.종.isin(["고릴라", "보노보", "침팬지", "보르네오오랑우탄", "수마트라오랑우탄", "시아망"])].itertuples(index=False):
    d = R[R.묶음 == t.묶음]
    for g in ("all", "q1", "q4"):
        dd = d[d.무리 == g]
        ing = dd[(dd.위치 >= t.유전자시작) & (dd.위치 <= t.유전자끝)]
        out = dd[dd.위치 < t.유전자시작]
        gb = ing.메틸수.sum() / max(ing.호출수.sum(), 1)
        sb = out.메틸수.sum() / max(out.호출수.sum(), 1)
        P("  %-28s %-4s 유전자 %.3f · 스페이서 %.3f · 차 %+.3f  (자리 %d)"
          % (t.묶음, g, gb, sb, gb - sb, len(dd)))
LOG.close()
