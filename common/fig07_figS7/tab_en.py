# -*- coding: utf-8 -*-
"""English column names and values for the Figure 7 family result tables.

The analysis code works with its original (Korean) labels internally. This module maps those labels to
English at the file boundary:

  to_en_file(path)   called by a producer on a table it has just written; rewrites the file in English
                     (header and values), keeping every number, and returns the path.
  open_ko(path)      called by code that reads one of these tables; returns a text stream with the
                     original labels, so the reading code needs no other change.

The mapping is one to one inside each table, so open_ko(to_en_file(x)) gives back the text of x.
A file that is not listed in T passes through unchanged."""
import gzip, io, os, re

HG = re.compile(r"[ㄱ-ㆎ가-힣]")

SPECIES = {"제브라피시": "zebrafish", "애기장대": "Arabidopsis thaliana", "금화조": "zebra finch",
           "수마트라오랑우탄": "Sumatran orangutan", "보르네오오랑우탄": "Bornean orangutan",
           "고릴라": "gorilla", "보노보": "bonobo", "침팬지": "chimpanzee", "시아망": "siamang"}
SIDE = {"낮은 좌표": "lower coordinate", "높은 좌표": "higher coordinate"}
YESNO = {"예": "yes", "아니오": "no"}
MODAL = {"단봉": "unimodal", "이봉": "bimodal"}
CEN_METHOD = {
    "위성 밀도 봉우리 (UCSC GenArk repeatMasker Satellite >=50%/100kb 연속 · 유전체 중앙의 40~259 배)":
        "satellite-density peak (UCSC GenArk repeatMasker Satellite >=50%/100kb run, 40~259 times the genome median)",
    "주석된 중심절 영역 (Eichler T2T centromeric_regions · Satellite 범위 · 두 하플로 합침)":
        "annotated centromeric region (Eichler T2T centromeric_regions, Satellite range, the two haplotypes merged)",
    "서열에서 측정 (CEN180 주기 178 bp 탠덤 덩어리 · 일치율>=0.80)":
        "measured from sequence (tandem block with the CEN180 period of 178 bp, autocorrelation>=0.80)",
    "서열에서 측정 (Tgut191A 191 bp · Tgut716A 716 bp 의 배수 주기 · 일치율>=0.55 · 10 kb 창)":
        "measured from sequence (period a multiple of Tgut191A 191 bp or Tgut716A 716 bp, autocorrelation>=0.55, 10 kb windows)"}
# a pair holding {} is a template: each {} stands for a part that is carried over, itself mapped if it matches a pair
DETAIL = [("검정불가(귀무 적합 폭주)", "not testable (null fit diverged)"),
          ("단봉(2성분 유의하지 않음)", "unimodal (two components not significant)"),
          ("판정불가(분해능 부족 · 간격이 1/CpG 의 {} 배)", "not callable (insufficient resolution, component gap {} times 1/CpG)"),
          ("단봉(봉우리 하나 · 어깨 골깊이비 {})", "unimodal (one peak, shoulder trough-depth ratio {})"),
          ("단봉(봉우리 하나)", "unimodal (one peak)"),
          ("단봉+소수꼬리(소수 성분 {})", "unimodal with a minor tail (minor component {})"),
          ("이봉", "bimodal"),
          ("보류(재현 {}/{})", "undecided (reproduced {}/{})")]
METHOD = [("우도비(이봉 · p {}/{} · 컷 {} · LRT p<={} · LR {})", "likelihood ratio (bimodal, p {}/{}, cut {}, LRT p<={}, LR {})"),
          ("사람 5S 비율 목표 f={} · 가장 가까운 컷 {} (달성 {} · {})", "human 5S fraction target f={}, nearest cut {} (achieved {}, {})")] + DETAIL
VERDICT = {"검정불가": "not testable", "단봉": "unimodal", "판정불가": "not callable", "이봉": "bimodal",
           "단봉+소수꼬리": "unimodal with a minor tail", "보류": "undecided"}

T = {
 "CENTROMERE.tsv": {
     "cols": {"종": "species", "접근번호": "accession", "염색체": "chromosome", "시작": "start", "끝": "end", "방법": "method",
              "시작Mb": "start_Mb", "끝Mb": "end_Mb", "길이Mb": "length_Mb"},
     "vals": {"종": SPECIES, "방법": CEN_METHOD}},
 "CHROM_NAMES.tsv": {
     "cols": {"접근번호": "accession", "종": "species", "염색체번호": "chromosome_number", "출처": "source",
              "우리표_길이": "chrom_map_length", "길이일치": "length_match"},
     "vals": {"종": SPECIES, "길이일치": {"확인": "confirmed", "우리표에 없음": "not in the chromosome map", "조회값 없음": "no NCBI value"}}},
 "FIGDATA_A1_chrom_map.tsv": {
     "cols": {"종": "species", "종키": "species_key", "염색체": "chromosome", "염색체이름": "chromosome_name",
              "염색체길이": "chromosome_length", "시작": "start", "끝": "end", "사본수": "n_copies", "유닛길이": "unit_length",
              "유전자길이": "gene_length", "분석에씀": "used_in_analysis"},
     "vals": {"종": SPECIES}},
 "FIGDATA_SUPP4_ATH_CONTEXT.tsv": {
     "cols": {"문맥": "context", "read유닛": "read_units", "자리수중앙": "median_n_sites", "평균β": "mean_beta"}},
 "ORIENTATION.tsv": {
     "cols": {"종": "species", "염색체길이": "chromosome_length", "array시작": "array_start", "array끝": "array_end",
              "낮은끝까지": "to_low_end", "높은끝까지": "to_high_end", "가까운끝": "nearer_end", "원위쪽": "distal_side",
              "근위쪽": "proximal_side", "가까운끝까지_Mb": "to_nearer_end_Mb", "먼끝까지_Mb": "to_farther_end_Mb",
              "염색체내_상대위치": "relative_position_in_chromosome", "중심절시작": "centromere_start", "중심절끝": "centromere_end",
              "중심절까지_Mb": "to_centromere_Mb", "중심절겹침": "overlaps_centromere", "근위근거": "proximal_basis",
              "끝자와같나": "vs_nearest_end_rule"},
     "vals": {"종": SPECIES, "가까운끝": SIDE, "원위쪽": SIDE, "근위쪽": SIDE, "중심절겹침": YESNO,
              "근위근거": {"중심절 (염색체이름)": "centromere (chromosome name)", "중심절 (접근번호)": "centromere (accession)",
                           "중심절 (염색체이름(하플로 떼고))": "centromere (chromosome name, haplotype suffix removed)"},
              "끝자와같나": {"같음": "same", "다름": "different"}}},
 "PERARRAY_LANDSCAPE.tsv.gz": {
     "cols": {"종": "species", "이름표": "label", "쪽": "side", "위치": "location", "호출수": "n_calls", "메틸호출수": "n_meth_calls"},
     "vals": {"종": SPECIES, "쪽": SIDE, "위치": {"array 안": "inside the array"}}},
 "PERARRAY_META.tsv": {
     "cols": {"종": "species", "이름표": "label", "array시작": "array_start", "array끝": "array_end", "유닛길이": "unit_length",
              "유닛수": "n_units", "유전자시작": "gene_start", "유전자끝": "gene_end", "부분유닛_낮은쪽": "partial_unit_low_side",
              "부분유닛_높은쪽": "partial_unit_high_side", "유닛1쪽": "unit1_side", "read유닛_ALL": "read_units_ALL",
              "read유닛_HIGH": "read_units_HIGH", "read유닛_LOW": "read_units_LOW", "read유닛_MID": "read_units_MID"},
     "vals": {"종": SPECIES, "유닛1쪽": SIDE}},
 "PROFILE_BP_v13.tsv.gz": {
     "cols": {"종": "species", "묶음": "group", "무리": "class", "유닛길이": "unit_length", "위치": "position",
              "호출수": "n_calls", "메틸수": "n_meth"},
     "vals": {"종": SPECIES}},
 "PROFILE_META_v13.tsv": {
     "cols": {"종": "species", "묶음": "group", "유닛길이": "unit_length", "원유닛": "raw_units", "완전유닛": "complete_units",
              "틀맞춤유닛": "frame_matched_units", "지도유닛": "mapped_units", "버림_CpG적음": "discarded_few_CpG",
              "버림_덮임부족": "discarded_low_coverage", "흔한자리": "common_positions", "최대자리관측": "max_position_observations",
              "CpG중앙": "CpG_median", "참조CpG중앙": "reference_CpG_median", "q1컷": "q1_cut", "q4컷": "q4_cut",
              "무리방법": "class_method", "판정": "verdict", "판정상세": "verdict_detail", "LRT_귀무최대": "LRT_null_max",
              "성분_저": "component_low", "성분_고": "component_high", "성분_저몫": "component_low_weight",
              "소수성분몫": "minor_component_weight", "간격배_분해능": "gap_over_resolution",
              "반반_저차최대": "halves_low_diff_max", "반반_고차최대": "halves_high_diff_max", "반반_몫차최대": "halves_weight_diff_max",
              "반반허용": "halves_tolerance", "재현벌": "replicates_passed", "적합봉수": "fitted_peaks", "골깊이비": "trough_depth_ratio",
              "한성분_알파": "one_component_alpha", "한성분_베타": "one_component_beta", "옛판정_v12": "old_verdict_v12",
              "three조건": "three_conditions", "귀무통과율": "null_pass_rate", "과분산귀무통과율": "overdispersed_null_pass_rate",
              "과분산귀무_알파": "overdispersed_null_alpha", "과분산귀무_베타": "overdispersed_null_beta",
              "과분산귀무_쓸수있나": "overdispersed_null_usable", "분할방식": "split_method", "달성LOW비율": "achieved_LOW_fraction",
              "소수성분": "minor_component", "컷": "cut", "유전자시작": "gene_start", "유전자끝": "gene_end"},
     "vals": {"종": SPECIES, "무리방법": METHOD, "판정": VERDICT, "판정상세": DETAIL, "옛판정_v12": MODAL,
              "three조건": {"통과": "pass", "탈락": "fail"},
              "과분산귀무_쓸수있나": {"예": "yes", "아니오(α·β<1 이라 귀무가 U자)": "no (α, β < 1, so the null is U-shaped)"},
              "분할방식": {"우도비": "likelihood ratio", "가장 가까운 컷": "nearest cut"}}},
 "UNITBETA_PERARRAY_v13.tsv.gz": {
     "cols": {"종": "species", "묶음": "group", "유닛번호": "unit_number", "시작좌표": "start_coordinate", "무리": "class"},
     "vals": {"종": SPECIES}},
}

_RX = {}


def _rx(t):
    if t not in _RX:
        _RX[t] = re.compile("^" + "(.+?)".join(re.escape(p) for p in t.split("{}")) + "$")
    return _RX[t]


def _cell(x, exact, tmpl):
    if x in exact:
        return exact[x]
    for a, b in tmpl:
        m = _rx(a).match(x)
        if m:
            parts = b.split("{}")
            return "".join(p + (_cell(m.group(k + 1), exact, tmpl) if k < len(parts) - 1 else "") for k, p in enumerate(parts))
    return x


def _pairs(m, to_en):
    p = list(m.items()) if isinstance(m, dict) else list(m)
    if not to_en:
        p = [(b, a) for a, b in p]
    return {a: b for a, b in p if "{}" not in a}, [(a, b) for a, b in p if "{}" in a]


def _convert(text, name, to_en):
    t = T[name]
    cols = t["cols"] if to_en else {v: k for k, v in t["cols"].items()}
    lines = text.split("\n")
    head = lines[0].split("\t")
    ko = head if to_en else [cols.get(h, h) for h in head]
    maps = [_pairs(t["vals"][h], to_en) if h in t.get("vals", {}) else None for h in ko]
    out = ["\t".join(cols.get(h, h) for h in head)]
    for ln in lines[1:]:
        c = ln.split("\t")
        out.append("\t".join(_cell(x, *maps[k]) if k < len(maps) and maps[k] else x for k, x in enumerate(c)))
    return "\n".join(out)


def _read(path):
    if path.endswith(".gz"):
        with gzip.open(path, "rt", encoding="utf-8", newline="") as f:
            return f.read()
    with io.open(path, encoding="utf-8", newline="") as f:
        return f.read()


def to_en_text(text, name):
    out = _convert(text, name, True)
    if HG.search(out):
        raise ValueError("Korean left in %s: %r" % (name, HG.findall(out)[:5]))
    if _convert(out, name, False) != text:
        raise ValueError("round trip failed for %s" % name)
    return out


def to_en_file(path):
    name = os.path.basename(path)
    if name not in T:
        return path
    e = to_en_text(_read(path), name).encode("utf-8")
    if path.endswith(".gz"):
        with io.open(path, "wb") as f, gzip.GzipFile(filename="", mode="wb", fileobj=f, mtime=0) as g:
            g.write(e)
    else:
        with io.open(path, "wb") as f:
            f.write(e)
    return path


def open_ko(path):
    name = os.path.basename(str(path))
    if name not in T:
        return path
    return io.StringIO(_convert(_read(path), name, False))
