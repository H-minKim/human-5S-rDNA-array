# -*- coding: utf-8 -*-
# v94 — trios recomputed with the v11-add2 source definition; the earlier v79b definition is kept alongside
#  원천 정의 (v11-add2_trio\09_tests\fix1_test_GIAB_precise\08_by_zone_both_ends_aligned_transmitted_vs_control.txt) :
#    β = beta_call(자리 단위 합본의 ML>=128 비율) · read>=2 · 구역 = 원위1-10 · 원위11-20 · 근위11-20 · 근위1-10 양 끝에서 센 번호로 맞춤
#    판정축 = 구역 안 |Δβ|(자식 vs 전달 부모 하플로) 중앙  vs  |Δβ|(자식 vs 비전달 = 같은 부모의 다른 하플로) 중앙
#    HG003(부) → HG002 짝은 전달 검정에서 제외
#  짝 : v79b 확정 짝(유닛 변이 서명 자카드 + 배경 대조). 비전달 하플로 = 같은 부모의 다른 슬롯(없으면 대조 없음)
#  [0] 원천 표 재현 : 같은 짝(HG01250 · HG01251 · HG03815 · HG00556 · HG004)에서 원천 08 표 값과 대조
#  [1] 모든 확정 짝에서 원천 정의로 · 구역별 전달 vs 비전달 (짝 단위 윌콕슨)
# [2] the result of the earlier definition (v79b : beta_read_mean, read>=3, consensus zone, unrelated control) is placed alongside
#  입력 v11-추가2 `08_methylation\READMAP_METH_V11F2+extra2_trios.tsv` · `03_maps\G8MAP_trio.tsv` · 02_tables/TRIO_PAIR_CONF.tsv · TRIO_TRANSMIT_ZONE2.tsv
#  출력 02_tables/TRIO_SRCDEF_PAIRZONE.tsv · TRIO_SRCDEF_SUMMARY.tsv · 04_logs/V94_TRIO_SRCDEF.txt
import io
import numpy as np, pandas as pd
from scipy import stats
OUT = __import__("os").environ["RDNA_DATA_ROOT"] + "/v11_readvar"; T = OUT + "/02_tables"; TR = __import__("os").environ["RDNA_DATA_ROOT"] + "/for_paper/08_versions/v11-add2_trio"
LOG = io.open(OUT + "/04_logs/V94_TRIO_SRCDEF.txt", "w", encoding="utf8", newline="\n")
def P(s=""): print(s); LOG.write(str(s) + "\n"); LOG.flush()
ME = pd.read_csv(TR + "/08_methylation/READMAP_METH_V11F2+extra2_trios.tsv", sep="\t", dtype={"slot": str})
ME["pri"] = np.where(ME.dataset.astype(str) != "v11", 0, 1)
ME = ME.sort_values("pri").drop_duplicates(["sample", "slot", "array_pos"], keep="first")   # trio rows take precedence
ME = ME[ME.n_reads >= 2]
M = pd.read_csv(TR + "/03_maps/G8MAP_trio.tsv", sep="\t", dtype={"slot": str}, usecols=["sample", "slot", "N"])
NN = M.drop_duplicates(["sample", "slot"]).set_index(["sample", "slot"]).N.astype(int).to_dict()
P("자리 단위 메틸(read>=2) %s 행 · 유닛 수 지도 %d 슬롯" % (format(len(ME), ","), len(NN)))
ZONES = [("원위1-10", "d", 1, 10), ("원위11-20", "d", 11, 20), ("근위11-20", "p", 11, 20), ("근위1-10", "p", 1, 10)]
def prof(s, sl):
    g = ME[(ME["sample"] == s) & (ME.slot == sl)]
    N = NN.get((s, sl))
    if g.empty or N is None: return None
    d = {}
    for r in g.itertuples(index=False):
        ad = int(r.array_pos); ap = N - ad + 1
        d[("d", ad)] = r.beta_call; d[("p", ap)] = r.beta_call
    return d
def zd(c, x, side, lo, hi):
    if c is None or x is None: return np.nan, 0
    v = [abs(c[(side, k)] - x[(side, k)]) for k in range(lo, hi + 1) if (side, k) in c and (side, k) in x]
    return (float(np.median(v)) if v else np.nan), len(v)
def zb(c, side, lo, hi):
    if c is None: return np.nan
    v = [c[(side, k)] for k in range(lo, hi + 1) if (side, k) in c]
    return float(np.median(v)) if v else np.nan

C = pd.read_csv(T + "/TRIO_PAIR_CONF.tsv", sep="\t", dtype={"child_slot": str, "parent_slot": str})
conf = C[C.call.str.startswith("전달")].copy()
conf = conf[~((conf.child.str.lower() == "hg002") & (conf.parent.str.lower() == "hg003"))]
rows = []
for r in conf.itertuples(index=False):
    others = sorted(set(C[(C.parent == r.parent)].parent_slot) - {r.parent_slot})
    c = prof(r.child, r.child_slot); t = prof(r.parent, r.parent_slot)
    n = prof(r.parent, others[0]) if others else None
    for zn, side, lo, hi in ZONES:
        dt, kt = zd(c, t, side, lo, hi); dn, kn = zd(c, n, side, lo, hi)
        rows.append(dict(child=r.child, child_slot=r.child_slot, parent=r.parent, role=r.role, trans_slot=r.parent_slot,
                         nontrans_slot=others[0] if others else "-", call=r.call, zone=zn, n_trans=kt, d_trans=dt, n_non=kn, d_non=dn,
                         closer=(dt < dn) if (dt == dt and dn == dn) else np.nan,
                         b_child=zb(c, side, lo, hi), b_trans=zb(t, side, lo, hi), b_non=zb(n, side, lo, hi)))
Z = pd.DataFrame(rows)
Z.to_csv(T + "/TRIO_SRCDEF_PAIRZONE.tsv", sep="\t", index=False, float_format="%.4f")

P()
P("[0] 원천 08 표 재현 (같은 짝 · 같은 정의)")
SRC = {("HG01250", "원위1-10"): (0.025, 0.041), ("HG01250", "원위11-20"): (0.112, 0.335), ("HG01250", "근위11-20"): (0.036, 0.190), ("HG01250", "근위1-10"): (0.075, 0.489),
       ("HG01251", "원위1-10"): (0.015, 0.010), ("HG01251", "원위11-20"): (0.139, 0.070), ("HG01251", "근위11-20"): (0.035, 0.026), ("HG01251", "근위1-10"): (0.099, 0.106),
       ("HG03815", "원위1-10"): (0.051, 0.010), ("HG03815", "원위11-20"): (0.096, 0.291), ("HG03815", "근위11-20"): (0.135, 0.494), ("HG03815", "근위1-10"): (0.034, 0.041),
       ("HG00556", "원위1-10"): (0.019, 0.013), ("HG00556", "원위11-20"): (0.064, 0.308), ("HG00556", "근위11-20"): (0.097, 0.095), ("HG00556", "근위1-10"): (0.046, 0.103),
       ("hg004", "원위1-10"): (0.052, 0.055), ("hg004", "원위11-20"): (0.444, 0.209), ("hg004", "근위11-20"): (0.026, 0.073), ("hg004", "근위1-10"): (0.042, 0.550)}
ok = 0; tot = 0
for (par, zn), (st, sn) in SRC.items():
    x = Z[(Z.parent.str.lower() == par.lower()) & (Z.zone == zn)]
    x = x[x.call == "전달(확정)"] if (x.call == "전달(확정)").any() else x
    if x.empty: P("   %-8s %-9s 원천 %.3f/%.3f · 재현 자료 없음" % (par, zn, st, sn)); continue
    r = x.iloc[0]; tot += 1
    m = abs(r.d_trans - st) < 0.0015 and abs(r.d_non - sn) < 0.0015
    ok += int(m)
    P("   %-8s %-9s 원천 전달 %.3f · 비전달 %.3f | 재현 %.3f · %.3f (유닛 %d/%d) %s" % (par, zn, st, sn, r.d_trans, r.d_non, r.n_trans, r.n_non, "일치" if m else "다름"))
P("   → %d / %d 칸 일치 (소수 셋째 자리)" % (ok, tot))

P()
P("[1] 원천 정의 · 모든 확정 짝 (HG003→HG002 제외)")
for (ch, cs, par), g in Z.groupby(["child", "child_slot", "parent"], sort=False):
    P("   %s h%s ← %s h%s (%s · 비전달 h%s)" % (ch, cs, par, g.trans_slot.iloc[0], g.role.iloc[0], g.nontrans_slot.iloc[0]))
    for r in g.itertuples(index=False):
        P("      %-9s 전달 %s · 비전달 %s → %s   (자식/전달/비전달 β %s / %s / %s)"
          % (r.zone, "%.3f (%d)" % (r.d_trans, r.n_trans) if r.d_trans == r.d_trans else "-", "%.3f (%d)" % (r.d_non, r.n_non) if r.d_non == r.d_non else "-",
             ("예 (%.1f배)" % (r.d_non / r.d_trans) if r.d_trans > 0 else "예") if r.closer is True else ("아니오" if r.closer is False else "대조 없음"),
             *["%.2f" % v if v == v else "-" for v in (r.b_child, r.b_trans, r.b_non)]))
P()
P("   구역별 요약 (비전달 대조가 있는 짝만)")
S = []
for zn, g in Z.groupby("zone", sort=False):
    g = g.dropna(subset=["d_trans", "d_non"])
    if len(g) < 3: continue
    p = stats.wilcoxon(g.d_trans, g.d_non)[1] if len(g) >= 5 else np.nan
    S.append(dict(zone=zn, n_pair=len(g), med_trans=g.d_trans.median(), med_non=g.d_non.median(),
                  frac_closer=float((g.d_trans < g.d_non).mean()), wilcoxon_p=p))
    P("   %-9s 짝 %d · 전달 중앙 %.3f · 비전달 중앙 %.3f · 전달이 더 가까운 짝 %.2f · 윌콕슨 p %s"
      % (zn, len(g), g.d_trans.median(), g.d_non.median(), (g.d_trans < g.d_non).mean(), "%.3f" % p if p == p else "(짝<5)"))
S = pd.DataFrame(S)
P()
P("[2] the earlier definition (v79b, alongside) : beta_read_mean · read>=3 · consensus zone · unrelated control")
O = pd.read_csv(T + "/TRIO_TRANSMIT_ZONE2.tsv", sep="\t")
po = O.pivot_table(index="zone", columns="kind", values="d", aggfunc="median")
for z in po.index:
    P("   %-10s 전달 %.3f · 무관 %.3f" % (z, po.loc[z].get("전달", np.nan), po.loc[z].get("무관 대조", np.nan)))
S.to_csv(T + "/TRIO_SRCDEF_SUMMARY.tsv", sep="\t", index=False, float_format="%.4f")
LOG.close()
