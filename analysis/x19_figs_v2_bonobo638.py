# -*- coding: utf-8 -*-
# 뒤 그림 v2 — 보노보를 guppy 6.3.8 read 만으로 다시 그림 (원 그림·원 코드 무수정)
#   원 그림 코드(x2 · x10 · x12)를 그대로 exec 하고 입력과 저장 이름만 바꾼다 (분석 논리·시드 그대로) :
#     R1  x2  : AVF 합산에서 온 array 행(meth 가 AVF frac 합산과 같은 행 — 애기장대 CG 재계산 행·붉은털·생쥐는 제외)
#               의 meth/flank/delta 를 정수 n_5mC 판(x16)으로 · 보노보는 6.3.8 만 값(x17) · 호출기 딱지 '6.3.8 만'
#     R4  x10 : read_unit 에서 보노보 6.3.7 read 제거(x18 표) · 유인원 페이지만
#     R5  x12 : coding 표에서 보노보 6.3.7 read 제거 · 유닛전체 판정은 v2 표에서 · 유인원 페이지만
#   거름 순서는 x17 과 같게(중앙값 문턱은 원래대로 전체 read 로 정한 뒤 6.3.7 을 뺌) → 수치가 x17 과 맞아야 함
# Saves figures/<name>_v2.png, .svg (svg.fonttype=none, text kept as text) and .pdf (fonttype 42)
#   표 : tables/OTHERSP_READUNIT_PERARRAY_v2_bonobo638.tsv · tables/OTHERSP_READUNIT_CODING_PERARRAY_v2_bonobo638.tsv
#   로그 : X19_FIGS_V2.txt (+ 각 원 코드 로그는 X19_x2.txt · X19_x10.txt · X19_x12.txt 로 바꿔 씀)
import io, os, shutil, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import os as _os_need, sys as _sys_need


def _need(k):
    """No defaults are set: if an environment variable is missing the script stops instead of silently reading the wrong location."""
    v = _os_need.environ.get(k)
    if not v:
        _sys_need.exit("environment variable %s is not set (the deposited copy has no default; "
                       "see README.md, section 'How to run it')" % k)
    return v


OUT = _need("RDNA_PREV_OUT_DIR")
HO = _os_need.path.dirname(_need("RDNA_SVG_DIR"))
SVGDIR = HO + "/30_figures_SVG"
LOG = io.open(OUT + "/X19_FIGS_V2.txt", "w", encoding="utf8", newline="\n")


def P(s=""):
    print(s)
    LOG.write(s + "\n")
    LOG.flush()


def _save(fig, name):
    base = OUT + "/figures/" + name
    fig.savefig(base + ".png")
    with matplotlib.rc_context({"svg.fonttype": "none", "pdf.fonttype": 42}):
        fig.savefig(base + ".svg", format="svg")
        fig.savefig(base + ".pdf", format="pdf")
    for ext in (".svg", ".pdf"):
        shutil.copy2(base + ext, SVGDIR + "/" + name + ext)
    P("  저장 : figures/%s .png .svg .pdf (+ 30_그림_SVG 사본)" % name)


def run(fname, reps, tag):
    src = io.open(OUT + "/code/" + fname, encoding="utf8").read()
    for a, b in reps:
        assert src.count(a) == 1, "%s : 바꿀 문장을 정확히 한 번 못 찾음 — %s" % (fname, a[:60])
        src = src.replace(a, b)
    ns = {"__name__": "x19_" + tag, "_save": _save}
    exec(compile(src, fname + "(v2)", "exec"), ns)
    plt.close("all")
    return ns


FILT = ("_BC = pd.read_csv(OUT + \"/tables/BONOBO_READ_BASECALLER.tsv\", sep=\"\\t\")\n"
        "_637 = set(_BC[_BC.basecaller == \"guppy-6.3.7\"].read_id)\n"
        "_n0 = len(D)\n"
        "D = D[~((D.species_key == \"bonobo_mPanPan1\") & D.read_id.isin(_637))].copy()\n"
        "print(\"  [v2] 보노보 6.3.7 read-유닛 %d 개 뺌\" % (_n0 - len(D)))\n")

# ---------------------------------------------------------------- R1 (x2)
OVR = ("T = pd.read_csv(OUT + \"/tables/OTHERSP_ARRAY_METH.tsv\", sep=\"\\t\")\n"
       "_N = pd.read_csv(OUT + \"/tables/OTHERSP_ARRAY_METH_n5mC.tsv\", sep=\"\\t\").set_index(\"array_id\")\n"
       "_m = (T.array_id.isin(_N.index)) & ((T.meth - T.array_id.map(_N.meth_frac)).abs() < 1e-4)\n"
       "T.loc[_m, \"meth\"] = T.loc[_m, \"array_id\"].map(_N.meth_int)\n"
       "T.loc[_m, \"flank\"] = T.loc[_m, \"array_id\"].map(_N.flank_int)\n"
       "T.loc[_m, \"delta\"] = T.loc[_m, \"meth\"] - T.loc[_m, \"flank\"]\n"
       "_BB = pd.read_csv(OUT + \"/tables/BONOBO_BY_BASECALLER.tsv\", sep=\"\\t\")\n"
       "_b = _BB[_BB.grp == \"6.3.8\"].set_index(\"metric\")\n"
       "_bo = T.species_key == \"bonobo_mPanPan1\"\n"
       "T.loc[_bo, \"meth\"] = _b.loc[\"A1 array β (AVF 정의)\", \"value\"]\n"
       "T.loc[_bo, \"flank\"] = _b.loc[\"A1 flank β (L+R)\", \"value\"]\n"
       "T.loc[_bo, \"delta\"] = T.loc[_bo, \"meth\"] - T.loc[_bo, \"flank\"]\n"
       "T.loc[_bo, \"n_call\"] = _b.loc[\"A1 array β (AVF 정의)\", \"n\"]\n"
       "print(\"  [v2] 정수 n_5mC 로 바꾼 행 %d · 보노보 6.3.8 만 meth %.4f flank %.4f delta %+.4f\" % "
       "(_m.sum(), T.loc[_bo, 'meth'].iloc[0], T.loc[_bo, 'flank'].iloc[0], T.loc[_bo, 'delta'].iloc[0]))\n")
ns = run("x2_fig_meth.py", [
    ('LOG = io.open(OUT + "/X2_FIG.txt", "w", encoding="utf8", newline="\\n")',
     'LOG = io.open(OUT + "/X19_x2.txt", "w", encoding="utf8", newline="\\n")'),
    ('T = pd.read_csv(OUT + "/tables/OTHERSP_ARRAY_METH.tsv", sep="\\t")\n', OVR),
    ('"보노보": "6.3.7+6.3.8 혼합"', '"보노보": "guppy6.3.8 만"'),
    # 표시 고침(수치 무관) : 원 그림은 chrom_name 이 비면 이름표가 "nan" — 새 pandas 의 astype(str) 가 NaN 을 "nan" 글자로
    #   바꾸지 않아 원 코드의 != "nan" 비교가 늘 참이 되던 것. fillna 로 chrom(접근번호/염색체)을 채움.
    ('labT["lab"] = np.where(labT.chrom_name.astype(str) != "nan",\n                       labT.chrom_name.astype(str), labT.chrom.astype(str))',
     'labT["lab"] = labT.chrom_name.fillna(labT.chrom).astype(str)'),
    ('fig.savefig(OUT + "/figures/Fig_othersp_array_meth.png")', '_save(fig, "Fig_othersp_array_meth_v2")'),
], "x2")
T = ns["T"]
g = T[T.species_key.isin(["bonobo_mPanPan1", "chimp_mPanTro3", "gorilla_mGorGor1", "orang_bornean_mPonPyg2", "orang_sumatran_mPonAbe1"])]
P("R1 v2 : 대형 유인원 5 array − flank %+.3f ~ %+.3f (보노보 %+.4f)" % (g.delta.min(), g.delta.max(), T[T.species_key == "bonobo_mPanPan1"].delta.iloc[0]))

# ---------------------------------------------------------------- R4 (x10) 유인원 페이지
ns = run("x10_readunit_per_array.py", [
    ('LOG = io.open(OUT + "/X10_READUNIT_PERARRAY.txt", "w", encoding="utf8", newline="\\n")',
     'LOG = io.open(OUT + "/X19_x10.txt", "w", encoding="utf8", newline="\\n")'),
    ('D = D[D.n_calls >= 0.5 * med].copy()\n', 'D = D[D.n_calls >= 0.5 * med].copy()\n' + FILT),
    ('T.to_csv(OUT + "/tables/OTHERSP_READUNIT_PERARRAY.tsv"', 'T.to_csv(OUT + "/tables/OTHERSP_READUNIT_PERARRAY_v2_bonobo638.tsv"'),
    ('for page, grps in PAGES:', 'PAGES = [p for p in PAGES if p[0] == "유인원"]\nfor page, grps in PAGES:'),
    ('fig.savefig(OUT + "/figures/Fig_readunit_perarray_%s.png" % page)', '_save(fig, "Fig_readunit_perarray_%s_v2" % page)'),
], "x10")
b = ns["T"][ns["T"].array_id.str.contains("bonobo")].iloc[0]
P("R4 v2 : 보노보 read-유닛 %d · OFF %.4f · %s · 봉우리 %.2f/%.2f (x17 : OFF 0.0875 · 이봉 0.89/0.13)" % (b.n, b.off, b.call, b.peak1, b.peak2))

# ---------------------------------------------------------------- R5 (x12) 유인원 페이지
ns = run("x12_coding_perarray.py", [
    ('LOG = io.open(OUT + "/X12_CODING_PERARRAY.txt", "w", encoding="utf8", newline="\\n")',
     'LOG = io.open(OUT + "/X19_x12.txt", "w", encoding="utf8", newline="\\n")'),
    ('D = D[(D.n_calls_gene >= np.maximum(3, 0.5 * med))].copy()\n', 'D = D[(D.n_calls_gene >= np.maximum(3, 0.5 * med))].copy()\n' + FILT),
    ('X10 = pd.read_csv(OUT + "/tables/OTHERSP_READUNIT_PERARRAY.tsv"', 'X10 = pd.read_csv(OUT + "/tables/OTHERSP_READUNIT_PERARRAY_v2_bonobo638.tsv"'),
    ('T.to_csv(OUT + "/tables/OTHERSP_READUNIT_CODING_PERARRAY.tsv"', 'T.to_csv(OUT + "/tables/OTHERSP_READUNIT_CODING_PERARRAY_v2_bonobo638.tsv"'),
    ('for page, grps in PAGES:', 'PAGES = [p for p in PAGES if p[0] == "유인원"]\nfor page, grps in PAGES:'),
    ('fig.savefig(OUT + "/figures/Fig_coding_readunit_perarray_%s.png" % page)', '_save(fig, "Fig_coding_readunit_perarray_%s_v2" % page)'),
], "x12")
b = ns["T"][ns["T"].array_id.str.contains("bonobo")].iloc[0]
P("R5 v2 : 보노보 coding read-유닛 %d · coding OFF %.4f · %s (x17 : 0.1217 · 이봉)" % (b.n, b.off, b.call))
LOG.close()
