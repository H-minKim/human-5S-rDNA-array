# -*- coding: utf-8 -*-
# The read x unit methylation distribution from the coding CpG sites alone, that is the 5S gene, separately per species and array.
#   Input: tables/OTHERSP_READUNIT_CODING.tsv.gz, built from the BAM MM/ML directly, with its whole-unit values checked against the existing read_unit table.
#   For comparison with human: b_gene and n_gene from L1, past the gate.
#   Conditions: a unit with a gene hit, coding CpG calls at least half that array's coding median, and at least 3.
#   The call follows the same rule as x10: a second peak at least 0.30 away, a trough no higher than half that peak, and at least 2% of molecules on its side.
#          The width is h = max(0.035, 0.5 / median CpG count), which widens where there are few coding CpG sites.
#   Also reported: the overdispersion ratio against a binomial null in which the molecules at a position share that position's mean p at the same number of calls.
#          OFF is beta < 0.35; the whole-unit and spacer beta of the same read-unit are recorded alongside, although the figure shows the coding part only.
import io, importlib.util
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams["font.family"] = "Malgun Gothic"
rcParams["axes.unicode_minus"] = False
rcParams["figure.dpi"] = 110
import os as _os_need, sys as _sys_need


def _need(k):
    """No defaults are set: if an environment variable is missing the script stops instead of silently reading the wrong location."""
    v = _os_need.environ.get(k)
    if not v:
        _sys_need.exit("environment variable %s is not set (the deposited copy has no default; "
                       "see README.md, section 'Environment variables')" % k)
    return v


OUT = _need("RDNA_PREV_OUT_DIR")
src = io.open(OUT + "/code/x10_readunit_per_array.py", encoding="utf8").read()
LOGLINE = 'LOG = io.open(OUT + "/X10_READUNIT_PERARRAY.txt", "w", encoding="utf8", newline="\\n")'
assert LOGLINE in src, "x10 log line not found: risk of overwriting the x10 log"
head = src.split("P(\"=\" * 104)")[0].replace(LOGLINE, "LOG = io.StringIO()")
G = {}
exec(head, G)
call = G["call"]
LOG = io.open(OUT + "/X12_CODING_PERARRAY.txt", "w", encoding="utf8", newline="\n")
NM = G["NM"]
rng = np.random.RandomState(11)


def P(s=""):
    try:
        print(s)
    except UnicodeEncodeError:
        print(s.encode("ascii", "replace").decode())
    LOG.write(s + "\n")
    LOG.flush()


def overdisp(b, n, key):
    p = pd.Series(b).groupby(key).transform("mean").values
    sim = rng.binomial(n, np.clip(p, 0, 1)) / n
    return np.var(b - p) / max(np.var(sim - p), 1e-12)


D = pd.read_csv(OUT + "/tables/OTHERSP_READUNIT_CODING.tsv.gz", sep="\t")
D = D.drop_duplicates(["species_key", "sample", "read_id", "array_id", "array_pos"])
D = D[D.n_gene_hits > 0].copy()
med = D.groupby("array_id").n_calls_gene.transform("median")
D = D[(D.n_calls_gene >= np.maximum(3, 0.5 * med))].copy()
D["bg"] = D.n_5mC_gene / D.n_calls_gene
D["ba"] = D.n_5mC / D.n_calls
nsp = D.n_calls - D.n_calls_gene
D["bs"] = np.where(nsp > 0, (D.n_5mC - D.n_5mC_gene) / nsp.clip(lower=1), np.nan)
D["종"] = D.species_key.map(NM)
ST = pd.read_csv(OUT + "/tables/OTHERSP_ARRAY_STRUCT.tsv", sep="\t")
UL = dict(zip(ST.array_id, ST.unit_len))
NU = dict(zip(ST.array_id, ST.n_unit))
X10 = pd.read_csv(OUT + "/tables/OTHERSP_READUNIT_PERARRAY.tsv", sep="\t")
X10c = dict(zip(X10.array_id, X10["call"]))

ROWS, PANEL = [], {}
H = pd.read_parquet(_need("RDNA_HUMAN_L1_PARQUET"),
                    columns=["gate", "b_gene", "n_gene", "b_call", "b_sp", "sample", "slot", "array_pos"])
H = H[(H.gate == 1) & (H.n_gene >= 5) & H.b_gene.notna()]
key = H["sample"].astype(str) + "_" + H.slot.astype(str) + "#" + H.array_pos.astype(str)
r, _, _ = call(H.b_gene.values, float(H.n_gene.median()))
r.update(group="사람", 종="사람 v11 (참고)", array_id="사람 510 하플로타입 합산", unit_len=2240, n_unit=np.nan,
         n=len(H), calls=float(H.n_gene.median()), off=(H.b_gene < 0.35).mean(), mean=H.b_gene.mean(),
         off_all=(H.b_call < 0.35).mean(), mean_sp=H.b_sp.mean(),
         od=overdisp(H.b_gene.values, H.n_gene.values.astype(int), key.values), call_all="이봉")
ROWS.append(r); PANEL[("사람", r["array_id"])] = (H.b_gene.values, r)

for (sp, aid), g in D.groupby(["종", "array_id"]):
    x = g.bg.values
    cm = float(g.n_calls_gene.median())
    r, _, _ = call(x, cm)
    grp = "유인원" if sp in ("보노보", "침팬지", "고릴라", "보르네오오랑", "수마트라오랑", "시아망") else sp
    k = (g.array_id + "#" + g.array_pos.astype(str)).values
    r.update(group=grp, 종=sp, array_id=aid, unit_len=UL.get(aid, np.nan), n_unit=NU.get(aid, np.nan),
             n=len(x), calls=cm, off=(x < 0.35).mean(), mean=x.mean(), off_all=(g.ba < 0.35).mean(),
             mean_sp=np.nanmean(g.bs), od=overdisp(x, g.n_calls_gene.values.astype(int), k),
             call_all=X10c.get(aid, "-"))
    if len(x) < 100:
        r["call"] = "표본 부족 (" + r["call"] + ")"
    ROWS.append(r); PANEL[(grp, aid)] = (x, r)

T = pd.DataFrame(ROWS)
T.to_csv(OUT + "/tables/OTHERSP_READUNIT_CODING_PERARRAY.tsv", sep="\t", index=False, float_format="%.5f")
P("=" * 100)
P("coding (5S gene) CpG only: verdict on the read x unit methylation distribution, per array")
P("=" * 100)
P("   species          array                      units  read-units codingCpG  coding_mean OFF(coding) OFF(whole_unit)  overdispersion  "
  "peak1  peak2  trough_ratio  minor_side  coding verdict      (whole-unit verdict)")
for r in T.sort_values(["group", "종", "mean"]).itertuples(index=False):
    P("   %-10s %-26s %5s %8s   %4.0f       %.3f     %.3f       %.3f      %5.1f  %.2f %s %s %.3f  %-14s (%s)"
      % (r.종, r.array_id.split("__")[-1][:26], ("%.0f" % r.unit_len) if r.unit_len == r.unit_len else "-",
         format(r.n, ","), r.calls, getattr(r, "mean"), r.off, r.off_all, r.od, r.peak1,
         ("%.2f" % r.peak2) if r.peak2 == r.peak2 else "  - ",
         ("%.2f" % r.trough_ratio) if r.trough_ratio == r.trough_ratio else "  - ", r.minor_mass, r.call, r.call_all))
P()
P("  verdict tally (coding)")
for grp, gg in T.groupby("group"):
    P("   %-10s array %2d : %s" % (grp, len(gg), " · ".join("%s %d" % (k, v) for k, v in gg.call.value_counts().items())))

COL = {"사람": "#c0392b", "유인원": "#1f5fa9", "금화조": "#e67e22", "애기장대": "#7e57c2", "제브라피시": "#2e7d32"}
PAGES = [("유인원", ["사람", "유인원"]), ("금화조", ["금화조"]), ("애기장대", ["애기장대"]), ("제브라피시", ["제브라피시"])]
for page, grps in PAGES:
    keys = [k for k in PANEL if k[0] in grps]
    keys.sort(key=lambda k: (k[0] != "사람", PANEL[k][1]["unit_len"] if PANEL[k][1]["unit_len"] == PANEL[k][1]["unit_len"] else 0,
                             PANEL[k][1]["종"], -PANEL[k][1]["n"]))
    ncol = 4 if len(keys) <= 12 else 6
    nrow = int(np.ceil(len(keys) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 2.9 * nrow + 1.0), squeeze=False)
    for i, k in enumerate(keys):
        ax = axs[i // ncol][i % ncol]
        x, r = PANEL[k]
        nb = int(round(r["calls"]))
        edges = np.clip((np.arange(nb + 2) - 0.5) / nb, 0, 1)
        ax.hist(x, bins=edges, weights=np.ones(len(x)) / len(x), color=COL[k[0]], alpha=0.85)
        ax.set_xlim(0, 1)
        acc = r["array_id"].split("__")[-1]
        acc = acc.split("_")[0] + ":" + (acc.split("_")[1] if "_" in acc else "")
        ttl = "%s  %s" % (r["종"], acc) if k[0] != "사람" else r["array_id"]
        ax.set_title("%s\n유닛 %s bp · read-유닛 %s · coding CpG %.0f\ncoding %s (OFF %.1f%%) | 유닛전체 %s (OFF %.1f%%)"
                     % (ttl, ("%.0f" % r["unit_len"]) if r["unit_len"] == r["unit_len"] else "-", format(r["n"], ","),
                        r["calls"], r["call"], 100 * r["off"], r["call_all"], 100 * r["off_all"]), fontsize=7.8)
        ax.tick_params(labelsize=7)
        if i % ncol == 0:
            ax.set_ylabel("read-유닛 비율", fontsize=8)
        if i // ncol == nrow - 1:
            ax.set_xlabel("read 안 유닛의 coding β", fontsize=8)
    for j in range(len(keys), nrow * ncol):
        axs[j // ncol][j % ncol].axis("off")
    fig.suptitle("%s — coding(5S 유전자) CpG 만 : 각 read 의 각 유닛 메틸 분포 · array 마다 한 칸" % page, fontsize=13)
    fig.tight_layout(rect=[0, 0.02, 1, 0.97])
    fig.savefig(OUT + "/figures/Fig_coding_readunit_perarray_%s.png" % page)
    plt.close(fig)
    P("  saved: figures/Fig_coding_readunit_perarray_%s.png" % page)
LOG.close()
