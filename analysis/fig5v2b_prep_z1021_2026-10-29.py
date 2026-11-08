# -*- coding: utf-8 -*-
# 복제 — 구역 경계만 근위 17→10 · 원위 22→21 로
#   원본은 그대로 · 출력은 TABLE_DIR(모래상자)로 간다 · 셈·관문·분모는 하나도 안 바꿨다.
# Fig 5 (v2) 자료 — 구역을 정본 8 구역으로 다시 센다 (근위가 왼쪽)
#   ① 먼저 원 표(TRIO_SRCDEF_PAIRZONE)를 합본 READMAP 으로 재현 해 방향·셈법을 확인한다
#   ② 그 다음 같은 셈법으로 8 구역을 만든다
#   구역 정본 : Z1 근1-3 · Z2 근4-10 · Z3 근11-17 · ZI 내부 · Z4 원11-22 · Z5 원9-10 · Z6 원3-8 · Z7 원1-2
import io, os, sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np, pandas as pd
from scipy import stats
# 자료 뿌리는 환경변수 하나로 (2026-10-21 · 저장소 공개 준비)
#   안 주면 지금까지 쓰던 자리를 그대로 쓴다 — 돌아가는 것이 달라지지 않는다.
DATA_ROOT = os.environ["DATA_ROOT"]
T = DATA_ROOT + "/02_tables/"
# 표를 읽고 쓰는 자리 — 저장소에서는 표가 코드 옆에 없으므로 환경변수로 받는다 (2026-10-21)
HERE = os.environ.get("TABLE_DIR", os.path.dirname(os.path.abspath(__file__)))
LOG = io.open(os.path.join(HERE, "FIG5V2_PREP.txt"), "w", encoding="utf8", newline="\n")
def P(s=""): print(s, flush=True); LOG.write(str(s) + "\n"); LOG.flush()

PZ0 = pd.read_csv(T + "TRIO_SRCDEF_PAIRZONE.tsv", sep="\t")
pairs = PZ0[["child", "child_slot", "parent", "role", "trans_slot", "nontrans_slot", "call"]].drop_duplicates()
P("[짝] %d (자식 %d · 가족 %d)" % (len(pairs), pairs.child.nunique(), pairs.parent.nunique()))

ids = set(pairs.child) | set(pairs.parent)
R = pd.read_csv(DATA_ROOT + "/02_tables/READMAP_METH_ALL_v2_2026-09-27.tsv", sep="\t",
                usecols=["sample", "data_source", "slot", "array_pos", "n_reads", "beta_call"])
R = R[R["sample"].isin(ids)].copy()
P("[자리 단위 β] 합본에서 트리오 %d 개체 · 행 %s · 자료원 %s"
  % (R["sample"].nunique(), format(len(R), ","), sorted(R.data_source.unique())))
R["hap"] = R["sample"].astype(str) + "|" + R.slot.astype(str)
NM = R.groupby("hap").array_pos.max().to_dict()
R["N"] = R.hap.map(NM)

# ---------- ① 방향 확인 : 원위 끝이 더 꺼져 있어야 한다 ----------
a1 = R[R.array_pos == 1].beta_call.median()
aN = R[R.array_pos == R.N].beta_call.median()
P()
P("[방향] array_pos 1 의 β 중앙 %.3f · array_pos N 의 β 중앙 %.3f → %s"
  % (a1, aN, "array_pos 1 = 원위 (정본 규칙과 같음)" if a1 < aN else "반대 — 아래 셈법을 뒤집어야 한다"))
FLIP = a1 >= aN
R["ad"] = R.N - R.array_pos + 1 if FLIP else R.array_pos      # 원위 끝에서 센 번호
R["ap"] = R.array_pos if FLIP else R.N - R.array_pos + 1      # 근위 끝에서 센 번호

# ---------- 원 표 재현 (4 구역) ----------
def zone4(ap, ad):
    if ap <= 10: return "근위1-10"
    if ap <= 20: return "근위11-20"
    if ad <= 10: return "원위1-10"
    if ad <= 20: return "원위11-20"
    return None
def zone8(ap, ad):
    if ap <= ad:
        if ap <= 3: return "Z1 근위1-3"
        if ap <= 10: return "Z2 근위4-10"
        if ap <= 17: return "Z3 근위11-17"
        return "ZI 내부"
    if ad <= 2: return "Z7 원위1-2"
    if ad <= 8: return "Z6 원위3-8"
    if ad <= 10: return "Z5 원위9-10"
    if ad <= 22: return "Z4 원위11-22"
    return "ZI 내부"

R["z4"] = [zone4(a, d) for a, d in zip(R.ap, R.ad)]
def zone5(ap, ad):
    """공통 5 구역 (2026-10-16) = 정본 8 구역을 묶은 것
       Z1 / Z2+Z3 / ZI / Z4 / Z5+Z6+Z7. 새 경계를 만들지 않았다.
       main (E)(F)(G) · FigX26 과 같은 조건식 이어야 한다."""
    if ap <= ad:
        if ap <= 3:  return "근위 1-3"
        if ap <= 10: return "근위 4-10"
        return "내부"
    if ad <= 10: return "원위 1-10"
    if ad <= 21: return "원위 11-21"
    return "내부"

R["z8"] = [zone8(a, d) for a, d in zip(R.ap, R.ad)]
R["z5"] = [zone5(a, d) for a, d in zip(R.ap, R.ad)]
MINREAD = 2
Rq = R[R.n_reads >= MINREAD].copy()

def zone_beta(hap, zcol):
    v = Rq[Rq.hap == hap]
    return v.groupby(zcol).beta_call.median()

def pair_table(zcol, zorder):
    rows = []
    for _, p in pairs.iterrows():
        kid = "%s|%d" % (p.child, p.child_slot)
        tr = "%s|%s" % (p.parent, p.trans_slot)
        nt = None if str(p.nontrans_slot).strip() == "-" else "%s|%s" % (p.parent, p.nontrans_slot)
        bk, bt = zone_beta(kid, zcol), zone_beta(tr, zcol)
        bn = zone_beta(nt, zcol) if nt else None
        for z in zorder:
            if z not in bk.index or z not in bt.index: continue
            r = dict(child=p.child, child_slot=p.child_slot, parent=p.parent, role=p.role,
                     call=p.call, zone=z, b_child=bk[z], b_trans=bt[z],
                     d_trans=abs(bk[z] - bt[z]))
            if bn is not None and z in bn.index:
                r["b_non"] = bn[z]; r["d_non"] = abs(bk[z] - bn[z])
            rows.append(r)
    return pd.DataFrame(rows)

Z4ORD = ["근위1-10", "근위11-20", "원위11-20", "원위1-10"]
P4 = pair_table("z4", Z4ORD)
P()
P("[① 재현] 원 표(4 구역)와 합본 재계산 대조 — |Δβ| 중앙")
P("   %-12s %10s %10s %10s %10s" % ("구역", "원표 전달", "재현 전달", "원표 비전달", "재현 비전달"))
for z in Z4ORD:
    o = PZ0[PZ0.zone == z]; n = P4[P4.zone == z]
    P("   %-12s %10.4f %10.4f %10.4f %10.4f"
      % (z, o.d_trans.median(), n.d_trans.median(),
         o.d_non.median(), n.d_non.median() if "d_non" in n else np.nan))
P("   ⇒ 값이 가까우면 셈법이 같은 것이다 (중앙값 비교 · 짝 %d vs %d)"
  % (PZ0[["child", "child_slot", "parent"]].drop_duplicates().shape[0], len(pairs)))

# ---------- ② 8 구역 ----------
Z8ORD = ["Z1 근위1-3", "Z2 근위4-10", "Z3 근위11-17", "ZI 내부",
         "Z4 원위11-22", "Z5 원위9-10", "Z6 원위3-8", "Z7 원위1-2"]
P8 = pair_table("z8", Z8ORD)
P8.to_csv(os.path.join(HERE, "FIG5V2_PAIRZONE8.tsv"), sep="\t", index=False, float_format="%.5f")
rows = []
P()
P("[② 8 구역] 자식 ↔ 전달 부모 하플로  vs  자식 ↔ 같은 부모의 비전달 하플로")
P("   %-14s %5s %10s %10s %8s %9s" % ("구역", "짝", "전달 |Δβ|", "비전달 |Δβ|", "더 가까움", "윌콕슨 p"))
for z in Z8ORD:
    g = P8[(P8.zone == z)].dropna(subset=["d_trans", "d_non"])
    if len(g) < 3:
        P("   %-14s %5d  (짝이 모자라 검정 안 함)" % (z, len(g)))
        rows.append(dict(zone=z, n_pair=len(g), med_trans=np.nan, med_non=np.nan,
                         closer=np.nan, p=np.nan)); continue
    w = stats.wilcoxon(g.d_trans, g.d_non)
    rows.append(dict(zone=z, n_pair=len(g), med_trans=g.d_trans.median(), med_non=g.d_non.median(),
                     closer=int((g.d_trans < g.d_non).sum()), p=w.pvalue))
    P("   %-14s %5d %10.4f %10.4f %5d/%-3d %9.4f"
      % (z, len(g), g.d_trans.median(), g.d_non.median(),
         int((g.d_trans < g.d_non).sum()), len(g), w.pvalue))
SM = pd.DataFrame(rows)
SM.to_csv(os.path.join(HERE, "FIG5V2_SUMMARY8.tsv"), sep="\t", index=False, float_format="%.6g")

# ---------- ②-b 공통 5 구역 (supple J 가 쓰는 것) ----------
#   정본 8 구역을 묶은 것 : Z1 / Z2+Z3 / ZI / Z4 / Z5+Z6+Z7 — 새 경계를 만들지 않았다.
Z5ORD = ["근위 1-3", "근위 4-10", "내부", "원위 11-21", "원위 1-10"]
P5 = pair_table("z5", Z5ORD)
P5.to_csv(os.path.join(HERE, "FIG5V2_PAIRZONE5.tsv"), sep="\t", index=False, float_format="%.5f")
rows5 = []
P()
P("[②-b 5 구역] 같은 셈법 · 구역만 공통 5 개로")
P("   %-12s %5s %10s %10s %8s %9s" % ("구역", "짝", "전달 |Δβ|", "비전달 |Δβ|", "더 가까움", "윌콕슨 p"))
for z in Z5ORD:
    g = P5[(P5.zone == z)].dropna(subset=["d_trans", "d_non"])
    if len(g) < 3:
        rows5.append(dict(zone=z, n_pair=len(g), med_trans=np.nan, med_non=np.nan,
                          n_closer=np.nan, p=np.nan))
        P("   %-12s %5d  (짝이 모자라 검정 안 함)" % (z, len(g))); continue
    w = stats.wilcoxon(g.d_trans, g.d_non)
    rows5.append(dict(zone=z, n_pair=len(g), med_trans=g.d_trans.median(),
                      med_non=g.d_non.median(),
                      n_closer=int((g.d_trans < g.d_non).sum()), p=w.pvalue))
    P("   %-12s %5d %10.4f %10.4f %5d/%-3d %9.4f"
      % (z, len(g), g.d_trans.median(), g.d_non.median(),
         int((g.d_trans < g.d_non).sum()), len(g), w.pvalue))
SM5 = pd.DataFrame(rows5)
SM5.to_csv(os.path.join(HERE, "FIG5V2_ZONE5.tsv"), sep="\t", index=False, float_format="%.6g")
P("   [저장] FIG5V2_PAIRZONE5.tsv · FIG5V2_ZONE5.tsv")

# ---------- ③ flank 전달은 잴 수 있나 ----------
P()
P("[③ flank 전달] 0–2 kb / 2–50 kb 가 전달되는지")
FH = pd.read_csv(T + "FLANK_METH_HAP_A0_CLEAN.tsv", sep="\t")
have = sorted(set(FH["sample"]) & ids)
miss = sorted(ids - set(FH["sample"]))
P("   flank β 표(FLANK_METH_HAP_A0_CLEAN)에 있는 트리오 개체 : %s" % have)
P("   없는 개체 : %s" % miss)
P("   ⇒ 부모 쪽 flank β 가 아예 없다 (부모는 1kGP-ONT·GIAB 이고 flank 창 계산은 HPRC 조립본에만 돌렸다)")
P("      따라서 '전달 하플로 vs 비전달 하플로' 를 flank 에서 견주는 것은 지금 자료로 불가능 하다.")
P("      가능해지려면 부모 ONT 자료에 flank 창 파이프라인을 돌려야 한다 (노드 작업).")
pd.DataFrame(dict(있음=[", ".join(have)], 없음=[", ".join(miss)])).to_csv(
    os.path.join(HERE, "FIG5V2_FLANK_AVAIL.tsv"), sep="\t", index=False)
P()
P("[저장] FIG5V2_PAIRZONE8.tsv · FIG5V2_SUMMARY8.tsv · FIG5V2_FLANK_AVAIL.tsv")
LOG.close()
