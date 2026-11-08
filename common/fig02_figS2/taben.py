# -*- coding: utf-8 -*-
r"""English column names and values for the Figure 2 and Supplementary Figure S2 result tables (2026-11-06).

The analysis code still works with its original (Korean) labels internally. This module is the single place
that maps those labels to English at the file boundary:

  to_en_file(path)   called by a producer right after it writes one of the tables below, rewrites the file in
                     English (header, values and comment lines), keeping every number and every byte of layout.
  open_ko(path)      called by code that reads one of these tables, returns a text stream with the original
                     labels, so the reading code needs no other change.
  read_csv_ko(path)  the same for pandas.

The mapping is exact and one to one inside each table, so open_ko(to_en_file(x)) gives back x byte for byte.
Tables that are not listed here pass through unchanged."""
import io, os, re

HG = re.compile(r"[가-힣]")

ZONE_COL = {"근위 1-3": "proximal_1-3", "근위 4-17": "proximal_4-17", "가운데": "central", "원위 1-10": "distal_1-10"}
ZONE_VAL = {"근위 1-3": "proximal 1-3", "근위 4-17": "proximal 4-17", "가운데": "central", "원위 1-10": "distal 1-10"}
FINE = {"C>T|비CpG": "C>T|nonCpG", "G>A|비CpG": "G>A|nonCpG", "C>T|CpG": "C>T|CpG", "G>A|CpG": "G>A|CpG",
        "결실 1bp": "del_1bp", "결실 2-9bp": "del_2-9bp", "결실 >=10bp": "del_>=10bp",
        "삽입 1bp": "ins_1bp", "삽입 2-9bp": "ins_2-9bp", "삽입 >=10bp": "ins_>=10bp",
        "블록 결실": "block_del", "블록 삽입": "block_ins", "기타": "other"}
FINE.update({s: s for s in ("C>A", "G>T", "C>G", "G>C", "T>C", "A>G", "T>G", "A>C", "T>A", "A>T")})
FINE_VAL = {"C>T|비CpG": "C>T|non-CpG", "G>A|비CpG": "G>A|non-CpG", "결실 1bp": "deletion 1 bp", "결실 2-9bp": "deletion 2-9 bp",
            "결실 >=10bp": "deletion >=10 bp", "삽입 1bp": "insertion 1 bp", "삽입 2-9bp": "insertion 2-9 bp",
            "삽입 >=10bp": "insertion >=10 bp", "블록 결실": "repeat-block deletion", "블록 삽입": "repeat-block insertion", "기타": "other"}
COARSE = {"1bp 결실": "1bp_deletion", "반복블록 결실": "repeat_block_deletion", "산화 C>A·G>T": "oxidation_C>A_or_G>T",
          "탈아미노 비CpG": "deamination_nonCpG", "탈아미노 CpG": "deamination_CpG", "그 밖 전이 T>C·A>G": "other_transitions_T>C_or_A>G",
          "그 밖 전환": "other_transversions", "긴 결실": "long_deletion", "1bp 삽입": "1bp_insertion", "긴 삽입": "long_insertion",
          "반복블록 삽입": "repeat_block_insertion"}
COARSE_VAL = {"1bp 결실": "1 bp deletion", "반복블록 결실": "repeat-block deletion", "산화 C>A·G>T": "oxidation C>A or G>T",
              "탈아미노 비CpG": "deamination non-CpG", "탈아미노 CpG": "deamination CpG", "그 밖 전이 T>C·A>G": "other transitions T>C or A>G",
              "그 밖 전환": "other transversions", "긴 결실": "long deletion", "1bp 삽입": "1 bp insertion", "긴 삽입": "long insertion",
              "반복블록 삽입": "repeat-block insertion"}
LAYER = {"전부": "all", "공유형 뺌": "shared removed", "단일유닛형만": "single-unit only"}
LAYER_COL = {"전부": "all", "공유형 뺌": "shared_removed", "단일유닛형만": "single_unit_only"}
ORIGIN = {"다개체 공유형": "shared", "배열내 다중형": "within_array_multiple", "단일유닛형": "single_unit"}
PROFILE = {"자리": "position", "쪽": "side", "접합부로부터_유닛": "units_from_junction", "하플로타입수": "n_haplotypes",
           "유닛당_변이수_합": "variants_per_unit_total", "출처": "source", "메틸베타": "meth_beta", "층": "layer"}


def _profile_col(h):
    if h in PROFILE: return PROFILE[h]
    for pre, en in (("절대_", "abs_"), ("몫_", "frac_")):
        if h.startswith(pre) and h[len(pre):] in FINE: return en + FINE[h[len(pre):]]
    if h.startswith("거친_") and h[3:] in COARSE: return "coarse_" + COARSE[h[3:]]
    for ko, en in LAYER_COL.items():
        if h.startswith(ko + "_") and h[len(ko) + 1:] in FINE: return en + "_" + FINE[h[len(ko) + 1:]]
    return None


T = {
 "X12_HAP_GROUPS.tsv": {"cols": {"무리": "group", "무리크기": "group_size", "보유토큰수": "n_tokens_carried"}},
 "X12_POP_VAR_DOWNSAMPLE.tsv": {"cols": {"반복": "replicate", "검정가능토큰": "testable_tokens"}},
 "X12_SIG_TOKENS.tsv": {"cols": {"무리": "group", "구역": "zone", "토큰": "token", "무리안_보유율": "carriage_in_group",
                                 "무리밖_보유율": "carriage_outside_group", "차": "difference"},
                        "vals": {"구역": {"가운데": "central", "근1_3": "prox1_3", "근4_17": "prox4_17", "원1_10": "dist1_10"}}},
 "X12_TYPE_GROUPS.tsv": {"cols": {"기준": "reference", "무리": "group", "하플로": "haplotypes", "개체": "individuals", "주집단": "main_population",
                                  "주집단몫": "main_population_fraction", "집단수": "n_populations", "CN중앙": "CN_median", "깊이중앙": "depth_median",
                                  "근위중앙": "proximal_median", "가운데중앙": "central_median", "원위중앙": "distal_median",
                                  "pc1중앙": "pc1_median", "pc2중앙": "pc2_median"},
                         "vals": {"기준": {"다수형": "majority"}}},
 "VTZONE_zones.tsv": {"cols": dict({"유형": "class", "최대_log2비": "max_log2_ratio", "최대_절대차": "max_abs_diff", "최대_귀무z": "max_null_z",
                                    "가장_먼_구역": "farthest_zone"}, **ZONE_COL),
                      "vals": {"유형": COARSE_VAL, "가장_먼_구역": ZONE_VAL}},
 "VTZONE4_byzone.tsv": {"cols": {"층": "layer", "구역": "zone", "무리": "group", "유닛당": "per_unit", "귀무z": "null_z", "구역안_순위": "rank_in_zone",
                                 "가운데대비_log2비": "log2_ratio_vs_central", "비_순위": "ratio_rank"},
                        "vals": {"층": LAYER, "구역": ZONE_VAL, "무리": FINE_VAL}},
 "VTSTRAT_cpg_profile.tsv": {"colfn": _profile_col},
 "fig_VTORG_single_region_type_composition_source_values.tsv": {"colfn": _profile_col,
                        "vals": {"출처": {"v11_io.load('primary') · 2단계 평균": "v11_io.load('primary'), two-stage mean"}}},
 "fig_VTPD_region_type_composition_source_values.tsv": {"colfn": _profile_col,
                        "vals": {"출처": {"v11_io.load('primary') · 2단계 평균": "v11_io.load('primary'), two-stage mean"}}},
 "S26b_boundary_drivers.tsv": {"cols": dict([("n_" + k, "n_" + v) for k, v in ORIGIN.items()] + [("dens_diff_" + k, "dens_diff_" + v) for k, v in ORIGIN.items()]),
                               "vals": {"set": {"귀무": "null", "실제": "observed"}, "location": dict(ZONE_VAL, **{"전체": "all"})}},
 "F8X2_pairs.tsv": {"cols": {"부류": "class", "간격": "gap", "짝수": "n_pairs", "하플로수": "n_haplotypes", "자카드평균": "jaccard_mean",
                             "공유자리평균": "shared_sites_mean", "합집합평균": "union_mean"},
                    "vals": {"부류": {"말단형-말단형": "terminal-terminal", "말단형-비말단형": "terminal-nonterminal", "비말단형-비말단형": "nonterminal-nonterminal"}}},
 "COSHARE_PERHAP.tsv": {"vals": {"B_stratum": {"B 희귀 (<1%)": "B rare (<1%)", "B 중간 (1-10%)": "B intermediate (1-10%)", "B 흔함 (>=10%)": "B common (>=10%)"}},
                        "comments": {"# 희귀 변이 A(코호트 유닛 빈도 < 1%) 공동 공유 · H1 557 가운데 창(0.2-0.8 · 창 유닛 23 이상) · 생산 COSHARE_BSTRAT.py · cls 0 = A 를 함께 가진 쌍 · 1 = 같은 거리에서 한쪽만 가진 쌍 · bin 1-10·11-100 은 같은 사건 혼입이라 그림에서 뺀다":
                                     "# co-sharing of a rare variant A (cohort unit frequency below 1%), central window 0.2-0.8 with at least 23 window units, made by COSHARE_BSTRAT.py, cls 0 = pairs that both carry A, cls 1 = pairs at the same distance where only one unit carries A, bins 1-10 and 11-100 mix copies of one event and are not drawn"}},
 "REGCHAR_axis.tsv": {"cols": {"특징": "feature", "길이층": "length_stratum", "자리": "position", "쪽": "side", "깊이": "depth", "유닛수": "n_units",
                               "평균점수": "mean_score", "문턱넘는몫": "frac_above_threshold", "초과몫": "excess_frac"},
                      "vals": {"특징": {"근위다움": "proximalness", "원위다움": "distalness"}}},
 "REGCHAR_relcount.tsv": {"cols": {"항목": "item", "N칸": "N_bin", "N평균": "N_mean", "N상대": "N_relative", "개수평균": "count_mean",
                                   "개수상대": "count_relative", "하플로수": "n_haplotypes"},
                          "vals": {"항목": {"고정 개별변이 2235TG": "fixed single variant 2235TG",
                                            "고정 개별변이 D1096:1(말단형 진단)": "fixed single variant D1096:1 (terminal-type diagnostic)",
                                            "비례 개별변이 2205TC(U84C)": "proportional single variant 2205TC (U84C)",
                                            "가운데 유닛(양쪽 깊이>10)": "central units (both depths >10)",
                                            "가운데 유닛(양쪽 깊이>30)": "central units (both depths >30)",
                                            "가운데 유닛(양쪽 깊이>5)": "central units (both depths >5)",
                                            "근위다움 영역(깊이 10 안)": "proximal-character region (within depth 10)",
                                            "원위다움 영역(깊이 10 안)": "distal-character region (within depth 10)"}}},
}

T["POP_VAR_DOWNSAMPLE.tsv"] = T["X12_POP_VAR_DOWNSAMPLE.tsv"]   # the producer writes this name, the figure reads the X12_ copy

# eighteen more tables (the tables are written in English by the producing script,)
_UNREL = {"상동 (같은 개인)": "homologous (same individual)", "무관 · 같은 세부집단": "unrelated, same subpopulation",
          "무관 · 다른 세부집단": "unrelated, different subpopulation", "무관 · 같은 대집단": "unrelated, same superpopulation",
          "무관 · 다른 대집단": "unrelated, different superpopulation", "무관 · 전체": "unrelated, all"}
T.update({
 "METH_TERM_units.tsv": {"comments": {"# 말단형 vs 비말단형 유닛 메틸 · H1 557 · 하드콜 β(beta_call · dedup pool) · read>=5 · 유닛 20 개 이상 하플로타입 · beta_norm = β - 그 하플로타입 유닛 β 중앙값 · 말단형 = D1096:1 보유":
                                      "# methylation of terminal-type vs non-terminal units, H1 557, hard-call beta (beta_call, dedup pool), read>=5, haplotypes with at least 20 units, beta_norm = beta minus the median unit beta of that haplotype, terminal type = carries D1096:1"}},
 "X2_zone_region.tsv": {"comments": {"# 방향 : revcomp Y -> HIGH_IS_PROX · N -> HIGH_IS_PROX (U84C 실측)": "# orientation: revcomp Y -> HIGH_IS_PROX, N -> HIGH_IS_PROX (measured with U84C)"},
                        "cells": {"근위단": "proximal end", "중앙": "middle", "원위단": "distal end"}},
 "Q1_results.tsv": {"comments": {"# interior 슬롯 531 · 유닛 44969 · 상동 개체 236 · 서명 12150": "# interior slots 531, units 44969, homologous individuals 236, signatures 12150"},
                    "cells": {"부분": "partial", "T1 array내": "T1 within array", "T2 상동": "T2 homologous", "T3 무관(맞춤)": "T3 unrelated (matched)"}},
 "Q2b_results.tsv": {"comments": {"# P08 Q2b · ruler majref 6.868174 · seed 20261010 · nsim 2499(구간 999) · 2026-10-10": "# P08 Q2b, ruler majref 6.868174, seed 20261010, nsim 2499 (interval 999), 2026-10-10"}},
 "S14_length.tsv": {"cells": {"가운데구획": "central compartment", "근위구획": "proximal compartment", "원위구획": "distal compartment"}},
 "S26b_position_profile.tsv": {"cells": {"다개체 공유형": "multi-individual shared", "배열내 다중형": "within-array multiple", "단일유닛형": "single-unit"}},
 "T1_pair_summary.tsv": {"vals": {"relation": {"무관": "unrelated", "부모": "parent", "부모1": "parent1", "부모2": "parent2", "부(추론정본·제외)": "father (inferred reference, excluded)", "모": "mother"}}},
 "W15_minority_group_gaps.tsv": {"cells": {"원위 말단형": "distal terminal-type", "원위 그 밖의 소수형": "other distal minority types", "근위 소수형": "proximal minority type"}},
 "X12_POP_VAR_PERM.tsv": {"cells": {"통계량": "statistic", "관측": "observed", "귀무중앙": "null_median", "귀무최소": "null_min", "귀무최대": "null_max", "치환수": "n_permutations",
                                    "내려뽑기수": "n_subsamples", "가장작은P": "smallest_P", "하플로타입": "haplotypes", "개인": "individuals", "난수씨": "seed",
                                    "AFR 고유 토큰 중앙값": "AFR private-token median", "AFR − 다른 집단 최대 (중앙값)": "AFR minus maximum of the other populations (median)"}},
 "X12_WINDOWS.tsv": {"cells": {"창": "window", "창수": "n_windows"}},
 "X1_population_control.tsv": {"comments": {"# P08 X1 집단 통제 · ruler majref 6.868174 · seed 20261010 · B 2000 · N_perm 999 · 2026-10-10": "# P08 X1 population control, ruler majref 6.868174, seed 20261010, B 2000, N_perm 999, 2026-10-10"},
                               "cells": dict(_UNREL, **{"순열-세부집단": "permutation, subpopulation", "순열-대집단": "permutation, superpopulation"})},
 "X5_lag2_noncircular.tsv": {"comments": {"# P08 X5 · 순환하지 않는 정의 (medoid 거리 상위 사분위) · seed 20261010 · N_perm 999": "# P08 X5, non-circular definition (top quartile of medoid distance), seed 20261010, N_perm 999"}},
 "X7_homogenisation_evidence.tsv": {"cells": {"사유 (그 haplotype 에만)": "private (that haplotype only)", "매우 희귀 (<0.1%)": "very rare (<0.1%)", "희귀 (0.1-1%)": "rare (0.1-1%)",
                                              "중간 (1-10%)": "intermediate (1-10%)", "흔함 (>=10%)": "common (>=10%)"}},
 "X9_patch_floor.tsv": {"comments": {"# P08 X9 · 패치 길이의 '교환 없는 바닥' · X4 와 같은 설정 (interior · 쌍 4000) · seed 20261010": "# P08 X9, no-exchange floor of the patch length, same settings as X4 (interior, 4000 pairs), seed 20261010"},
                        "cells": dict(_UNREL)},
 "XA_decay_constant_robustness.tsv": {"comments": {"# P08 XA · 감쇠 상수 a 의 민감도 · seed 20261010 · B 2000 · 2026-10-10": "# P08 XA, sensitivity of the decay constant a, seed 20261010, B 2000, 2026-10-10"}},
 "Y4_cellline_summary.tsv": {"cells": {"없음": "none"}},
 "ZB_neighbour_divergence.tsv": {"cells": {"말단형–보통(이웃)": "terminal-ordinary (neighbor)", "보통–보통(원위 이웃)": "ordinary-ordinary (distal neighbor)", "가운데 이웃": "central neighbor", "근위 이웃": "proximal neighbor"}},
 "ZM_terminal_similarity.tsv": {"cells": {"말단형": "terminal-type", "포함": "included", "같은 사람": "same person", "다른 사람": "different person", "제외": "excluded", "보통": "ordinary"}},
})
T["P08_X2_zone_region.tsv"] = T["X2_zone_region.tsv"]   # the figure-table copy of the P08 output


def _maps(name, header):
    t = T[name]
    col = {}
    for h in header:
        if h in t.get("cols", {}): col[h] = t["cols"][h]
        elif "colfn" in t and t["colfn"](h): col[h] = t["colfn"](h)
    return col


_NUM = re.compile(r"\d+(?:\.\d+)?")
def _num_comment(ln, cm):
    """2026-11-07 : a '#' line that differs from a mapped line only in its numbers (same numbers in the same order on both sides)"""
    for k, v in cm.items():
        nk, nv = _NUM.findall(k), _NUM.findall(v)
        if not nk or nk != nv: continue
        pat = "^" + "".join(re.escape(p) + ("(\\d+(?:\\.\\d+)?)" if i < len(nk) else "") for i, p in enumerate(_NUM.split(k))) + "$"
        m = re.match(pat, ln)
        if not m: continue
        parts = _NUM.split(v)
        return "".join(p + (m.group(i + 1) if i < len(nk) else "") for i, p in enumerate(parts))
    return ln


def _convert(text, name, to_en):
    t = T[name]
    nl = "\r\n" if "\r\n" in text else "\n"
    lines = text.split(nl)
    end = lines[-1] == ""
    if end: lines = lines[:-1]
    cm = t.get("comments", {})
    if not to_en: cm = {v: k for k, v in cm.items()}
    cells = t.get("cells", {})          # exact whole-cell map applied to every cell (any line, any column)
    if not to_en: cells = {v: k for k, v in cells.items()}
    out, header, col, vals = [], None, None, None
    for ln in lines:
        if ln.startswith("#"):
            out.append(cm[ln] if ln in cm else _num_comment(ln, cm)); continue
        c = ln.split("\t")
        if header is None:
            header = c
            if to_en:
                col = _maps(name, c); hdr_ko = c
            else:
                inv = {}
                for h_en in c:
                    for k, v in _maps_all(name).items():
                        if v == h_en: inv[h_en] = k
                col = inv; hdr_ko = [inv.get(h, h) for h in c]
            vm = t.get("vals", {})
            vals = []
            for h in hdr_ko:
                m = vm.get(h)
                vals.append((m if to_en else {v: k for k, v in m.items()}) if m else None)
            out.append("\t".join(cells.get(y, y) for y in (col.get(h, h) for h in c))); continue
        out.append("\t".join(cells.get(y, y) for y in ((vals[j].get(x, x) if j < len(vals) and vals[j] else x) for j, x in enumerate(c))))
    return nl.join(out) + (nl if end else "")


_ALL = {}
def _maps_all(name):
    if name not in _ALL:
        t = T[name]; m = dict(t.get("cols", {}))
        if "colfn" in t:
            for k in list(PROFILE) + ["%s%s" % (p, f) for p in ("절대_", "몫_") for f in FINE] + ["거친_" + c for c in COARSE] + \
                     ["%s_%s" % (l, f) for l in LAYER_COL for f in FINE]:
                v = t["colfn"](k)
                if v: m[k] = v
        _ALL[name] = m
    return _ALL[name]


def to_en_text(text, name):
    if name not in T: return text
    out = _convert(text, name, True)
    if HG.search(out): raise ValueError("Korean left in %s: %r" % (name, HG.findall(out)[:5]))
    if _convert(out, name, False) != text: raise ValueError("round trip failed for %s" % name)
    return out


def ko_text(text, name):
    if name not in T or not text: return text
    return _convert(text, name, False)


def to_en_file(path):
    name = os.path.basename(path)
    if name not in T: return path
    with io.open(path, encoding="utf-8", newline="") as f: s = f.read()
    e = to_en_text(s, name)
    with io.open(path, "w", encoding="utf-8", newline="") as f: f.write(e)
    return path


def open_ko(path, mode="r", encoding="utf8", **kw):
    name = os.path.basename(str(path))
    if name not in T: return io.open(path, mode, encoding=encoding, **kw)
    with io.open(path, encoding="utf-8", newline="") as f: s = f.read()
    k = ko_text(s, name)
    if kw.get("newline", None) != "": k = k.replace("\r\n", "\n")
    return io.StringIO(k)


def read_csv_ko(path, *a, **kw):
    import pandas as pd
    return pd.read_csv(open_ko(path, newline=""), *a, **kw)
