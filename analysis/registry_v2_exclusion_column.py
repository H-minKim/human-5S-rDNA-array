# -*- coding: utf-8 -*-
# The last stage of the reference registry: two columns, an exclusion flag and its reason, are added to produce the _v2 table.
#
#   The paths come from environment variables.
#
#   The exclusion list is empty: the deposited _v2 table has the flag set to N in all 616 rows,
#     and section 8 of that folder's README gives the reason.
#     ('NA18565', '2') is listed under `_SUPERSEDED` and is not excluded
#     (all five spanning reads give 37-38 units against 35 canonical, so the assembly has collapsed a unit).
#
#   Check: the table this writes is byte-identical to the deposited one (md5 92cadda22dae).
import io, os, sys, collections
sys.stdout.reconfigure(encoding="utf-8")

WORK_ROOT = os.environ["WORK_ROOT"]
REG_DIR = os.environ.get("REGISTRY_V11_DIR",
                         os.path.join(WORK_ROOT, "00_handover",
                                      "handover_v11_data_ledger_2026-10-04"))
IN = os.path.join(REG_DIR, "V11F4_DATASET_hap_flags.tsv")
OUT = os.environ.get("REGISTRY_V2_OUT",
                     os.path.join(REG_DIR, "V11F4_DATASET_hap_flags_v2.tsv"))

# Items that are not excluded and are not put in `EXCL`
_SUPERSEDED = {
    ("NA18565", "2"): "관통 read 5개가 전부 37~38 유닛(정본 35) — 조립본 유닛 접힘"
                      " · 2026-10-05 사용자 결정으로 되돌림(제외하지 않음)",
}
EXCL = {}

f = io.open(IN, encoding="utf8")
h = f.readline().rstrip("\n").split("\t")
ix = {c: i for i, c in enumerate(h)}
R = [L.rstrip("\n").split("\t") for L in f if L.strip()]
out = []
for r in R:
    k = (r[0], r[1])
    out.append(r + [("Y" if k in EXCL else "N"), EXCL.get(k, "-")])
with io.open(OUT, "w", encoding="utf8", newline="\n") as g:
    g.write("\t".join(h + ["excluded", "exclusion_reason"]) + "\n")
    g.write("\n".join("\t".join(r) for r in out) + "\n")
print("-> %s · 열 %d · 행 %d · 제외 %d"
      % (os.path.basename(OUT), len(h) + 2, len(out),
         sum(1 for r in out if r[-2] == "Y")))

# -- whether the denominator changes, which is the point of the exclusion column --
jx = {c: i for i, c in enumerate(h + ["excluded", "exclusion_reason"])}
fl = lambda r: set(r[ix["flags"]].split(";")) if r[ix["flags"]] != "-" else set()


def cnt(sel):
    s = collections.Counter(r[0] for r in sel)
    return len(sel), len(s), sum(1 for v in s.values() if v == 2)


H = [r for r in out if r[ix["cohort"]] == "HPRC"]
E = [r for r in out if r[ix["qc_eligible"]] == "1"]
C = [r for r in out if r[ix["ont_data_source"]] == "HPRC_ONT_submissions"]
F2 = [r for r in H if not ({"T1_single_spanning_read", "C2_read_consensus"} & fl(r))]
D = [r for r in out if int(r[ix["meth_rows"]]) > 0]
ok = lambda L2: [r for r in L2 if r[jx["excluded"]] == "N"]
print()
print("| 분모 | 전 (슬롯/개체/짝) | 후 |")
print("|---|---|---|")
for nm, L2 in (("B_HPRC정본", H), ("C_ONT출처", C), ("D_메틸있음", D),
               ("E_QC적격", E), ("F2_조립본만", F2)):
    a, b = cnt(L2), cnt(ok(L2))
    print("| %s | %d / %d / %d | %d / %d / %d |" % ((nm,) + a + b))
