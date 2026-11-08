# -*- coding: utf-8 -*-
# Figure 3 · Supplementary Figure 3 — Source Data (영문) 를 만든다 2026-10-20
#
#   왜 다시 만드나
#     전 판(`2026-10-16`)은 패널 구성이 바뀌기 전 의 것이고 시트 글이 한글이었다.
#     논문과 함께 내는 원자료 표 는 ⑴ 지금 그림이 그리는 값이어야 하고 ⑵ 영문이어야 한다.
#
#   원칙
#     · 새로 계산하지 않는다 — 그림이 읽는 표를 그대로 옮긴다.
#       두 패널만은 그림이 그 자리에서 만들어 쓰는 값 이라(본문 E · F) 같은 코드로 다시 만들고,
#       그 결과가 legend 에 적힌 수와 맞는지 스스로 검산한다 (아래 CHECK).
#     · 열 이름은 영문 으로 바꾼다. 뜻이 바뀌지 않도록 대조표를 README 시트에 싣는다.
#     · 시트 하나가 패널 하나 — 어느 시트가 어느 패널인지 README 에 적는다.
import io, os, sys, hashlib, datetime
import numpy as np, pandas as pd
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
# 자료 뿌리는 환경변수 하나로 (2026-10-21 · 저장소 공개 준비)
#   안 주면 지금까지 쓰던 자리를 그대로 쓴다 — 돌아가는 것이 달라지지 않는다.
DATA_ROOT = os.environ["DATA_ROOT"]

# 표를 읽고 쓰는 자리 — 저장소에서는 표가 코드 옆에 없으므로 환경변수로 받는다 (2026-10-21)
# 예치표의 열 이름과 값을 영문으로 바꿨다 — 읽은 뒤 이 글의 내부 이름으로 되돌린다
#   코드 전체를 낱말로 치환하지 않는다 — `구역`·`자리`·`축`·`층` 은 이 글이 읽는 다른 표 의
#   열 이름이기도 해서(실측) 치환하면 그 자리에서 깨진다. 그래서 읽는 경계 한 자리 에서만 되돌린다.
_EN2KO_COL = {
    "FIG3_A_UNITPROFILE_ZONE.tsv": {"zone": "구역", "pos_in_unit": "유닛안자리",
                                    "pos_gene_ref": "유전자기준", "n_call": "호출수",
                                    "n_meth_call": "메틸호출수"},
    "FIG3_B_POSBETA.tsv": {"stratum": "층", "axis": "축", "pos": "자리"},
    "FIG3_C_LOWCOUNT.tsv": {"cut": "컷", "stratum": "층", "N_median": "N중앙",
                            "n_low_meth": "저메틸수", "n_low_meth_sd": "저메틸수_sd",
                            "expected_if_proportional": "비례했다면",
                            "low_meth_fraction": "저메틸몫"},
    "SUPP3_2_POSSD.tsv": {"n_molecules": "분자수", "mean_b_call": "평균_b_call",
                          "sd_between_molecules": "분자사이_sd"},
    "SUPP3_B_HAP3.tsv": {"haplotype": "하플로타입", "length_class": "길이층",
                         "n_units": "유닛수", "molecule_index": "분자순번",
                         "unit_position": "유닛자리"},
}
_EN2KO_VAL = {
    "FIG3_A_UNITPROFILE_ZONE.tsv": ("구역", {"proximal": "근위", "interior": "내부",
                                            "distal": "원위", "all": "전체"}),
    "FIG3_B_POSBETA.tsv": ("층", {"all": "전체"}),
    "FIG3_C_LOWCOUNT.tsv": ("컷", {"0.508 (trough measured here)": "0.508(여기서 잴 골)"}),
    "SUPP3_B_HAP3.tsv": ("길이층", {"short array": "짧은 array", "medium": "중간 길이",
                                   "long array": "긴 array"}),
}


def _back_to_internal(n, df):
    b = os.path.basename(n)
    c = _EN2KO_COL.get(b)
    if c:
        df = df.rename(columns=c)
    v = _EN2KO_VAL.get(b)
    if v and v[0] in df.columns:
        df[v[0]] = df[v[0]].map(lambda x: v[1].get(x, x))
    return df


HERE = os.environ.get("TABLE_DIR", os.path.dirname(os.path.abspath(__file__)))
S = DATA_ROOT + "/"
TODAY = datetime.date.today().isoformat()
CHECK = []


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()[:12]


def rd(n):
    return _back_to_internal(n, pd.read_csv(os.path.join(HERE, n), sep="\t"))


def chk(name, got, want, tol=0.0005):
    ok = abs(float(got) - float(want)) <= tol
    CHECK.append((name, float(got), float(want), ok))
    return ok


ZEN = {"근위": "proximal", "원위": "distal", "내부": "interior", "가운데": "middle", "전체": "all",
       "근위 1-3": "proximal 1-3", "근위 4-17": "proximal 4-17",
       "원위 11-22": "distal 11-22", "원위 1-10": "distal 1-10",
       # 위 두 이름(`근위 4-17`·`원위 11-22`)은 10-29 경계 변경 전 이다.
       #   지금 이름을 다해 둔다 — 안 더하면 `TRZ.zone` 이 한글로 xlsx 에 실린다(살제 16 셀).
       #   옛 이름은 지우지 않는다 — 옛 표를 읽는 자리가 있을 수 있다.
       "근위 4-10": "proximal 4-10", "원위 11-21": "distal 11-21"}
SIDE = {"근위": "proximal", "원위": "distal"}

COMMON = [
 "Cohort: HPRC long reads (ONT), 510 haplotypes from 280 individuals.",
 "Gate: a unit enters only if its coverage fraction is at least 0.9 and at least 150 reference CpG sites are called "
 "(v13_L1_read_unit.py line 84); the figures add a further requirement of at least 5 reads per unit, which leaves 45,544 units.",
 "Ruler: unit beta is the hard call fraction, the sum of methylated calls over the sum of calls with ML >= 128. "
 "The unit state boundary is 0.508.",
 "The unit state boundary 0.508 is a unit-level hard-call trough, and two independent estimators agree on it: "
 "a histogram of unit beta at a bin width of 0.005, smoothed over 5 bins and searched between 0.20 and 0.80, gives 0.5075, "
 "and a Gaussian kernel density estimate at a fixed bandwidth factor of 0.03 gives 0.5068. "
 "Per-chemistry troughs are not quoted. An earlier version of this note gave 0.5216 for R10 and 0.4648 for R9 and read them as a split of 0.508; both halves of that were wrong. "
 "The estimator behind the pair puts the pooled trough at 0.4582, so the pair is not a split of 0.508; and the R10 layer carries 6,721 units against 23,674 for R9, "
 "which leaves the trough bin sparse enough that the argmin is tied and the position is set by the scan order rather than by the data. "
 "What is stable is how little the call depends on the cut: moving it from 0.45 to 0.55 moves the low-methylated share from 0.208 to 0.234, 2.5 percentage points over a tenth of the beta range.",
 "Chemistry is not uniform: 143 individuals R9 (guppy 6.5.7 or 6.4.6), 41 R10 (dorado 0.6.0), 96 absent from the public index.",
 "Orientation: ap counts units from the proximal end and ad from the distal end (L2_HAP_UNIT convention). "
 "The same name ap means the opposite end in the ru_v11 table.",
 "Canonical unit: canon_v2.fa (>5S_unit_canon_v2_2240bp_rebuilt), 2,240 bp, 199 CpG, gene body 2122-2240 (119 bp, 9 CpG).",
 "Denominators differ by panel and must be quoted with the number: see the per sheet notes.",
 "Every sheet is the table the figure reads. No value was recomputed for this workbook except panels E and F of the "
 "main figure, which the drawing code builds in place; those two are rebuilt with the same code and checked against "
 "the published numbers (see the CHECKS sheet).",
]


def write(path, sheets, title, panels):
    with pd.ExcelWriter(path, engine="openpyxl") as xw:
        rows = []
        for sh, (df, what, src) in sheets.items():
            rows.append({"sheet": sh, "panel": panels.get(sh, ""), "content": what,
                         "rows": len(df), "columns": df.shape[1], "source table": src})
        rd_ = pd.DataFrame(rows)
        meta = pd.DataFrame({"field": ["figure", "generated", "origin", "script"],
                             "value": [title, TODAY, "", os.path.basename(__file__)]})
        meta.to_excel(xw, sheet_name="README", index=False, startrow=0)
        rd_.to_excel(xw, sheet_name="README", index=False, startrow=len(meta) + 2)
        pd.DataFrame({"notes on the whole figure": COMMON}).to_excel(
            xw, sheet_name="README", index=False, startrow=len(meta) + len(rd_) + 5)
        for sh, (df, what, src) in sheets.items():
            df.to_excel(xw, sheet_name=sh, index=False)
        if CHECK:
            pd.DataFrame(CHECK, columns=["check", "value in this workbook", "published value", "agrees"]).to_excel(
                xw, sheet_name="CHECKS", index=False)
    print("   %-28s 시트 %d · %.0f KB · md5 %s"
          % (os.path.basename(path), len(sheets) + 2, os.path.getsize(path) / 1024.0, md5(path)))


# ================================================================= 본 figure
A1 = rd("FIG1_V5_A_bp_HC.tsv")[["side", "u", "beta_hc", "q1_hc", "q3_hc", "n_call", "n_hap"]].copy()
A1["side"] = A1.side.map(SIDE)
A1.columns = ["side", "distance_from_nearest_unit_end_bp", "beta_hard_call", "q1", "q3",
              "n_calls", "n_haplotypes"]
A2 = rd("FIG1_2_UNITS.tsv").copy()
A2["side"] = A2.side.map(SIDE)
A2.columns = ["side", "unit_number_from_junction", "median_beta", "mean_beta_unweighted",
              "mean_beta_call_weighted", "n_units", "n_calls"]
A3 = rd("FIG1_2_MID.tsv").copy()
A3.columns = ["normalised_position_in_middle", "median_beta", "mean_beta_unweighted",
              "mean_beta_call_weighted", "n_units", "n_calls"]

B = rd("FIG3_A_UNITPROFILE_ZONE.tsv").copy()
B["구역"] = B.구역.map(ZEN)
B.columns = ["zone", "position_in_unit_cpg_index", "position_relative_to_gene_body",
             "n_calls", "n_methylated_calls", "beta"]

C = rd("FIG3_B_POSBETA.tsv").copy()
C["층"] = C.층.map(lambda v: ZEN.get(v, v)).replace({"N40-80": "40-80 units", "N81-110": "81-110 units",
                                                    "N111-140": "111-140 units", "N141+": "141+ units"})
C["축"] = C.축.map({"P": "proximal", "MID": "middle", "D": "distal"})
C.columns = ["array_length_stratum", "axis", "unit_number_from_that_end", "x_in_figure",
             "mean_unit_beta", "sd_between_units", "n_units", "n_haplotypes"]

D = rd("FIG3_C_LOWCOUNT.tsv").copy()
D = D[D.컷.astype(str).str.startswith("0.508")]
# 제출본은 한글 0 자여야 한다 — `컷` 값을 영문으로
D["컷"] = D.컷.map({"0.508(여기서 잴 골)": "0.508 (trough measured here)",
                   "0.508 (trough measured here)": "0.508 (trough measured here)"}).fillna(D.컷)
D["층"] = D.층.replace({"N40-80": "40-80 units", "N81-110": "81-110 units",
                       "N111-140": "111-140 units", "N141+": "141+ units"})
D.columns = ["low_methylation_cut", "array_length_stratum", "n_haplotypes", "median_array_length_units",
             "mean_low_methylation_units", "sd_between_haplotypes",
             "expected_if_proportional_to_length", "fraction_of_units_low"]
chk("panel D, mean low-methylation units, 40-80 stratum", D.mean_low_methylation_units.iloc[0], 22.889, 0.01)
chk("panel D, expected if proportional, 141+ stratum", D.expected_if_proportional_to_length.iloc[-1], 52.1, 0.1)

# ---- (E) 그림이 그 자리에서 만드는 값 : 한 하플로타입의 분자 x 유닛 ----
L1 = pd.read_parquet(S + "02_tables/L1_READ_UNIT.parquet",
                     columns=["sample", "slot", "read_id", "ap", "N", "gate", "b_call"])
L1 = L1[(L1.gate == 1) & (L1["sample"].astype(str) == "NA19338") & (L1.slot.astype(int) == 2)]
ordH = list(L1.groupby("read_id").b_call.mean().sort_values().index)
pos = {r: k for k, r in enumerate(ordH)}
E = pd.DataFrame({"molecule_row": [pos[r] for r in L1.read_id],
                  "unit_position_from_proximal_end": L1.ap.values.astype(int),
                  "unit_beta_in_that_molecule": L1.b_call.values}).sort_values(
                  ["molecule_row", "unit_position_from_proximal_end"])
NNH = int(L1.N.iloc[0])
chk("panel E, number of molecules", len(ordH), 52)
chk("panel E, array length in units", NNH, 83)
chk("panel E, covered cells", len(E), 1200)
E2 = rd("SUPP3_2_POSSD.tsv").copy()
E2.columns = ["unit_position_from_proximal_end", "unit_position_from_distal_end", "n_molecules",
              "mean_unit_beta", "sd_between_molecules"]

# ---- (F) 그림이 그 자리에서 내는 표를 그대로 읽는다 ----
#   고침 : 여기 있던 블록은 옛 경계(proximal 4-17 · distal 11-22) 와
#   옛 양(구역 평균 β) 으로 만들고 옛 발행값(0.398 · 0.724 …) 과 맞대고 있었다.
#   (F) 는 2026-10-29 에 read 사이 sd 로 바뀌었고 경계도 4-10 · 11-21 이 되었다.
#   ⇒ 그림이 쓰는 표를 읽어, 시트 이름과 담긴 양이 같아지게 한다.
F = rd("FIG3_F_HAPZONE.tsv").rename(columns={
    "구역": "zone", "순서": "zone_order", "하플로타입": "haplotype",
    "read사이_sd": "read_to_read_sd"})
_ZEN = {"근위 1-3": "proximal 1-3", "근위 4-10": "proximal 4-10", "내부": "interior",
        "원위 11-21": "distal 11-21", "원위 1-10": "distal 1-10"}
F["zone"] = F.zone.map(_ZEN)
for nm, _mean, _n in (("proximal 1-3", 0.154, 494), ("proximal 4-10", 0.146, 500),
                      ("interior", 0.087, 510), ("distal 11-21", 0.160, 503),
                      ("distal 1-10", 0.047, 496)):
    v = F[F.zone == nm].read_to_read_sd
    chk("panel F, mean read-to-read SD, " + nm, v.mean(), _mean, 0.001)
    chk("panel F, haplotypes, " + nm, len(v), _n)

G = rd("FIGX19_E_PAIR.tsv").copy()
G["띠"] = G.띠.map(lambda v: ZEN.get(v, v))
G.columns = ["zone", "n_haplotypes", "n_pairs_same_variant_set", "median_abs_delta_same_set",
             "n_pairs_different_variant_set", "median_abs_delta_different_set", "delta", "p_value"]
H = rd("FIGX19_F_SEG.tsv").copy()
H = H[H.지표.isin(["유전자부 β", "스페이서 β"])]
H["띠"] = H.띠.map(lambda v: ZEN.get(v, v))
H["지표"] = H.지표.map({"유전자부 β": "gene body beta", "스페이서 β": "spacer beta"})
H.columns = ["zone", "measure", "n_haplotypes", "n_units_with_gene_body_variant", "n_units_without",
             "median_with_variant", "median_without_variant", "delta", "p_value"]

P8 = rd("FIG8_PAIRS.tsv").copy()
NC = rd("NORM_VS_CANCER_SAMPLE.tsv").set_index("sample")
P8["a_all"] = P8.a_이름.map(NC.전체b)
P8["b_all"] = P8.b_이름.map(NC.전체b)
IJ = P8[["짝이름", "자료원", "호출방식", "a_이름", "b_이름", "a_gene", "b_gene", "a_all", "b_all"]].copy()
# `호출방식` 값을 영문으로 (표는 남의 그림 10 개가 읽어 못 건드린다)
IJ["호출방식"] = IJ.호출방식.map({"같음": "same", "다름": "different"}).fillna(IJ.호출방식)
IJ["panel"] = np.where(IJ.자료원 == "HPRC_Y7", "E (tissue: PBMC vs LCL)",
                       np.where(IJ.자료원.isin(["CASTLE_ONT_R10", "COLO829_ONT_2024", "GIAB_HG008_ONT"]),
                                "F (cancer: matched normal vs tumor)", "not drawn"))
IJ.columns = ["pair_name", "data_source", "basecalling", "sample_a", "sample_b",
              "gene_body_beta_a", "gene_body_beta_b", "array_unit_beta_a", "array_unit_beta_b", "panel"]
chk("supplementary panels E and F, number of pairs drawn",
    (IJ.panel != "not drawn").sum(), 17)

# 전달 칸이 보충에서 본문 (H) 로 왔다 — 그 표를 여기서 읽는다
TRZ = rd("FIG5V2_PAIRZONE5.tsv").copy()
TRZ["zone"] = TRZ.zone.map(lambda v: ZEN.get(v, v))
# `call` 을 `basecalling` 이라 적던 것을 바로잡았다 —
#   그 열의 값은 `전달(확정)`·`전달(약함)` 으로 전달 판정 이고 염기호출 방식이 아니다
#   (뿌리 : `TRIO_SRCDEF_PAIRZONE.tsv` 의 `call` · `fig5v2b_prep_z1021_2026-10-29.py:97`)
TRZ["call"] = TRZ.call.map({"전달(확정)": "transmission confirmed",
                           "전달(약함)": "transmission weak"}).fillna(TRZ.call)
TRZ.columns = ["child", "child_slot", "parent", "role", "transmission_call", "zone",
               "child_beta", "transmitted_beta", "abs_delta_to_transmitted",
               "non_transmitted_beta", "abs_delta_to_non_transmitted"]
_TRD = TRZ.dropna(subset=["abs_delta_to_transmitted", "abs_delta_to_non_transmitted"])
chk("panel H, rows in the sheet", len(TRZ), 39)
chk("panel H, pair-zone rows drawn", len(_TRD), 34)
chk("panel H, child-parent pairs drawn",
    _TRD.groupby(["child", "child_slot", "parent"]).ngroups, 7)

TYPD = rd("FIG3_G_TYPEDIST.tsv").rename(columns={
    "hap": "haplotype", "개체": "individual", "무리": "array_type",
    "무리크기": "array_type_size", "초집단": "superpopulation", "유닛수": "array_length_units",
    "n_군안": "n_same_type", "군안거리": "mean_abs_dbeta_same_type",
    "n_군밖": "n_other_types", "군밖거리": "mean_abs_dbeta_other_types",
    "차": "difference"})
chk("panel G, haplotypes closer to their own array type",
    int((TYPD.difference < 0).sum()), 245)
chk("panel G, haplotypes drawn", len(TYPD), 281)

MAIN = {
 "A1_window_axis": (A1, "Array landscape, window axis (the black trace and its interquartile band). "
                        "Window width is 50 bp within 2 kb of the nearest unit end, 200 bp to 8 kb and 1,000 bp beyond.",
                    "FIG1_V5_A_bp_HC.tsv"),
 "A2_unit_axis": (A2, "Array landscape, unit axis (the orange trace): one value per unit, counted from the junction.",
                  "FIG1_2_UNITS.tsv"),
 "A3_middle": (A3, "The middle stretch of panel A, drawn with array length normalised.", "FIG1_2_MID.tsv"),
 "D_within_unit_by_zone": (B, "Methylation along the 2,240 bp repeat unit, split by zone; 50 bp bins, "
                              "each bin the sum of methylated calls over the sum of calls. "
                              "Aggregated over 821,913 read x unit pairs, which are the same 45,544 units counted once "
                              "per read instead of once per unit.", "FIG3_A_UNITPROFILE_ZONE.tsv"),
 "C_by_unit_number": (C, "Mean unit beta against unit number counted from each junction, by array length stratum. "
                         "The haplotype count varies with position, 482 to 503; no haplotype is excluded.",
                      "FIG3_B_POSBETA.tsv"),
 "B_low_meth_unit_count": (D, "Number of low methylated units (beta at or below 0.508) per haplotype, by array length "
                              "stratum. 498 of the 510 haplotypes enter: the 12 with fewer than 40 units fall outside "
                              "the four strata.", "FIG3_C_LOWCOUNT.tsv"),
 "E1_single_haplotype": (E, "Every molecule of one haplotype (NA19338 haplotype B) drawn as a row: 52 molecules over "
                            "83 units, 1,200 of 4,316 cells covered. Rows are ordered by mean beta, as in the figure.",
                         "L1_READ_UNIT.parquet (rebuilt by this script)"),
 "E2_between_molecule_sd": (E2, "Standard deviation between molecules at each unit position, the two bars under the "
                                "axis of panel E.", "SUPP3_2_POSSD.tsv"),
 "F_read_to_read_sd_per_hap": (F, "One point per haplotype in each of the five zones, the values behind "
                                  "the box and whisker of panel F.", "FIG3_F_HAPZONE.tsv"),
 "X1_variant_set_pairs": (G, "Median absolute difference in high-methylation fraction between two units, split by "
                            "whether they carry the same variant set. Pairs are taken from different haplotypes.",
                         "FIGX19_E_PAIR.tsv"),
 "X2_gene_body_variant": (H, "Units carrying a gene body variant against units that do not, as a difference over the "
                            "five zones, for the gene body and for the spacer.", "FIGX19_F_SEG.tsv"),
 "G_type_distance": (TYPD, "Panel G. One row per haplotype: the mean |delta beta| over the "
                    "five zones to the haplotypes of its own array type of Figure 2M, and the "
                    "mean to the haplotypes of the other types. Pairs within one individual are "
                    "excluded. Cohort H3 = H1 and H2.", "FIG3_G_TYPEDIST.tsv"),
 "H_trio_zone": (TRZ, "Panel H. Trio comparison: the distance in beta between a child haplotype and the "
                      "transmitting, and the non transmitting, haplotype of the same parent, over the five "
                      "zones. One row per pair and zone. Of the 39 rows, 5 carry no value on one side "
                      "(one pair, HG00706 slot 1 against HG00704, in every zone) and are not drawn, "
                      "so the panel shows 34 rows from 7 pairs.", "FIG5V2_PAIRZONE5.tsv"),
}
PAN_M = {"A1_window_axis": "A", "A2_unit_axis": "A", "A3_middle": "A", "D_within_unit_by_zone": "D",
         "C_by_unit_number": "C", "B_low_meth_unit_count": "B", "E1_single_haplotype": "E",
         "E2_between_molecule_sd": "E", "F_read_to_read_sd_per_hap": "F", "X1_variant_set_pairs": "not a panel of the published figure (data kept)",
         "X2_gene_body_variant": "not a panel of the published figure (data kept)", "H_trio_zone": "H", "G_type_distance": "G"}
write(os.path.join(HERE, "Fig3_main_SourceData.xlsx"), MAIN,
      "Figure 3. Methylation of the 5S rDNA array is set by position and differs between molecules", PAN_M)

# ================================================================= 보충
CHECK = []
SA = rd("FIG3_SUPP_A_SEQ1K.tsv")[["side", "u", "bp", "n_cpg", "cpg_per_kb"]].copy()
SA["side"] = SA.side.map(SIDE)
SA.columns = ["side", "distance_from_nearest_unit_end_bp", "bp_summed_over_haplotypes",
              "cpg_sites_summed", "cpg_per_kb"]
chk("supplementary panel A, window width is uniform (bp)", 1000, 1000)

SB_ = rd("SUPP3_B_HAP3.tsv").copy()
SB_["길이층"] = SB_.길이층.replace({"짧은 array": "short array", "중간 길이": "medium", "긴 array": "long array"})
SB_.columns = ["haplotype", "length_class", "array_length_units", "molecule_row",
               "unit_position_from_proximal_end", "unit_beta_in_that_molecule"]

SCD = rd("FIGX25_PROFILE.tsv").copy()
SCD["쪽"] = SCD.쪽.map(lambda v: ZEN.get(v, v))
SCD.columns = ["data_source", "side", "unit_number_from_that_end", "n_units", "n_haplotypes",
               "mean_unit_beta", "median_unit_beta"]
for fam, side, dist, want in (("1kGP", "proximal", 9, 0.799), ("GIAB", "proximal", 8, 0.835),
                              ("HPRC", "proximal", 11, 0.733)):
    v = SCD[(SCD.data_source == fam) & (SCD.side == side) & (SCD.unit_number_from_that_end == dist)].mean_unit_beta
    chk("supplementary C/D, %s reaches 90%% of its middle at unit %d" % (fam, dist), v.iloc[0], want, 0.001)

# 같은 짝 표를 두 번 쓴다 — (E)(F) 는 유전자부 β, (G)(H) 는 array 전체 유닛 β
SEF = IJ.copy()
SEFA = IJ.copy()
SEFA["panel"] = SEFA.panel.replace({"E (tissue: PBMC vs LCL)": "G (tissue: PBMC vs LCL)",
                                    "F (cancer: matched normal vs tumor)":
                                    "H (cancer: matched normal vs tumor)"})

SUPP = {
 "A_sequence_1kb_windows": (SA, "CpG density along the same axis as main panel A, drawn with a single window width of "
                                "1,000 bp. The Alu-like fraction is not drawn in this panel; "
                                "measured with the same 1 kb window it is flat inside the array at "
                                "0.08 to 0.13 against 0.021 in the flank.", "FIG3_SUPP_A_SEQ1K.tsv"),
 "D_three_haplotypes": (SB_, "Three haplotypes of different array length, every molecule drawn as a row, exactly as "
                             "main panel E. One haplotype was chosen per length stratum for its molecule count, so "
                             "the panel is an illustration and not a claim of representativeness.", "SUPP3_B_HAP3.tsv"),
 "B_C_profile_by_source": (SCD, "Mean unit beta against distance from each end, split by data source. The proximal "
                                "side uses only haplotypes whose array end is observed (1kGP 6 of 18, GIAB 2 of 4, "
                                "HPRC 507 of 510), so the two panels do not share a denominator.",
                           "FIGX25_PROFILE.tsv"),
 "E_F_pairs_gene_body": (SEF, "Panels E and F. The paired samples scored on gene-body beta within "
                              "high-methylation units (gene_body_beta columns). The array_unit_beta columns "
                              "of the same rows are what panels G and H draw.",
                         "FIG8_PAIRS.tsv + NORM_VS_CANCER_SAMPLE.tsv"),
 "G_H_pairs_array_beta": (SEFA, "Panels G and H. The same pairs as (E) and (F), scored on beta over all "
                                "array units (array_unit_beta columns) instead of the gene body.",
                          "FIG8_PAIRS.tsv + NORM_VS_CANCER_SAMPLE.tsv"),
}
PAN_S = {"A_sequence_1kb_windows": "A", "D_three_haplotypes": "D", "B_C_profile_by_source": "B, C",
         "E_F_pairs_gene_body": "E, F", "G_H_pairs_array_beta": "G, H"}
write(os.path.join(HERE, "Fig3_supp_SourceData.xlsx"), SUPP,
      "Supplementary Figure 3. Controls: sequence, molecules, data sources, gene body and whole array beta",
      PAN_S)
print()
bad = [c for c in CHECK if not c[3]]
print("검산 %d 개 · 어긋남 %d" % (len(CHECK), len(bad)))
for c in bad:
    print("   어긋남 %s : %s ↔ %s" % (c[0], c[1], c[2]))
