# -*- coding: utf-8 -*-
# The later figures, version 2: the bonobo is redrawn from the guppy 6.3.8 reads alone, leaving the original figures and the original code untouched.
#   The original drawing code (x2, x10 and x12) is exec'd unchanged, with only the input and the output name altered, so the analysis and the seeds are the same:
#     R1, from x2: the array rows coming from the AVF sum, those whose methylation equals the AVF fraction sum, excluding the recomputed Arabidopsis CG row, the rhesus and the mouse,
#               have their methylation, flank and delta taken from the integer n_5mC version in x16, with the bonobo on its 6.3.8-only values from x17 and labelled accordingly.
#     R4, from x10: the bonobo 6.3.7 reads are removed from the read-unit table using the x18 table, for the ape pages only.
#     R5, from x12: the bonobo 6.3.7 reads are removed from the coding table, the whole-unit call is taken from the v2 table, for the ape pages only.
#   The order of filtering matches x17: the median threshold is first set on all the reads, and the 6.3.7 reads are then removed, so the numbers must agree with x17.
# Saves figures/<name>_v2.png, .svg (svg.fonttype=none, text kept as text) and .pdf (fonttype 42)
#   Tables: tables/OTHERSP_READUNIT_PERARRAY_v2_bonobo638.tsv and tables/OTHERSP_READUNIT_CODING_PERARRAY_v2_bonobo638.tsv
#   Logs: X19_FIGS_V2.txt, with each original script's log redirected to X19_x2.txt, X19_x10.txt and X19_x12.txt.
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
                       "see README.md, section 'Environment variables')" % k)
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
    P("  saved: figures/%s .png .svg .pdf" % name)


def run(fname, reps, tag):
    src = io.open(OUT + "/code/" + fname, encoding="utf8").read()
    for a, b in reps:
        assert src.count(a) == 1, "%s: the sentence to replace was not found exactly once: %s" % (fname, a[:60])
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
    # Display only, no number changes: an empty chrom_name is filled in with the chromosome, either an accession or a chromosome name,
    #   so that the label does not read "nan".
    ('labT["lab"] = np.where(labT.chrom_name.astype(str) != "nan",\n                       labT.chrom_name.astype(str), labT.chrom.astype(str))',
     'labT["lab"] = labT.chrom_name.fillna(labT.chrom).astype(str)'),
    ('fig.savefig(OUT + "/figures/Fig_othersp_array_meth.png")', '_save(fig, "Fig_othersp_array_meth_v2")'),
], "x2")
T = ns["T"]
g = T[T.species_key.isin(["bonobo_mPanPan1", "chimp_mPanTro3", "gorilla_mGorGor1", "orang_bornean_mPonPyg2", "orang_sumatran_mPonAbe1"])]
P("R1 v2: great apes 5 array − flank %+.3f ~ %+.3f (bonobo %+.4f)" % (g.delta.min(), g.delta.max(), T[T.species_key == "bonobo_mPanPan1"].delta.iloc[0]))

# ---------------------------------------------------------------- R4 (x10), the ape pages
ns = run("x10_readunit_per_array.py", [
    ('LOG = io.open(OUT + "/X10_READUNIT_PERARRAY.txt", "w", encoding="utf8", newline="\\n")',
     'LOG = io.open(OUT + "/X19_x10.txt", "w", encoding="utf8", newline="\\n")'),
    ('D = D[D.n_calls >= 0.5 * med].copy()\n', 'D = D[D.n_calls >= 0.5 * med].copy()\n' + FILT),
    ('T.to_csv(OUT + "/tables/OTHERSP_READUNIT_PERARRAY.tsv"', 'T.to_csv(OUT + "/tables/OTHERSP_READUNIT_PERARRAY_v2_bonobo638.tsv"'),
    ('for page, grps in PAGES:', 'PAGES = [p for p in PAGES if p[0] == "유인원"]\nfor page, grps in PAGES:'),
    ('fig.savefig(OUT + "/figures/Fig_readunit_perarray_%s.png" % page)', '_save(fig, "Fig_readunit_perarray_%s_v2" % page)'),
], "x10")
b = ns["T"][ns["T"].array_id.str.contains("bonobo")].iloc[0]
P("R4 v2: bonobo read-units %d, OFF %.4f, %s, peaks %.2f/%.2f (x17: OFF 0.0875, bimodal 0.89/0.13)" % (b.n, b.off, b.call, b.peak1, b.peak2))

# ---------------------------------------------------------------- R5 (x12), the ape pages
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
P("R5 v2: bonobo coding read-units %d, coding OFF %.4f, %s (x17: 0.1217, bimodal)" % (b.n, b.off, b.call))
LOG.close()
