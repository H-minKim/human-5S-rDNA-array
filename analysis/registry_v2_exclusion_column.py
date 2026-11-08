# -*- coding: utf-8 -*-
# 참조서열 대장의 마지막 단 — `제외` · `제외사유` 두 열을 붙여 `_v2` 를 만든다
#
#   이 글은 2026-10-04 에 셸 heredoc 으로 돌았고 파일로 저장되지 않았다.
#   셈은 바꾸지 않았다. 바꾼 것 : heredoc 껍데기를 벗기고 자리를 환경변수로 뺐다.
#
#   `EXCL` 에 대하여 — 여기에 판이 둘 있었다. 헷갈리지 않게 적는다.
#     기록에서 꺼낸 판은 `EXCL` 에 `('NA18565','2')` 가 들어 있다
#     (사유 : "관통 read 5 개가 전부 37~38 유닛(정본 35) — 조립본 유닛 접힘").
#     그런데 기탁된 `_v2` 는 616 행이 전부 `제외=N` 이다.
#     그 폴더 `README` §8 이 까닭을 적는다 :
#     ⇒ 기탁본에 맞는 것은 **`EXCL` 가 빈 판**이다. 그래서 여기서는 비워 둔다.
#        뺐던 항목은 지우지 않고 `_SUPERSEDED` 에 남긴다 — 있었던 일이다.
#
#   검산 : 이 글이 내는 표는 기탁본과 바이트까지 같다 (md5 92cadda22dae · 2026-10-22 확인).
import io, os, sys, collections
sys.stdout.reconfigure(encoding="utf-8")

WORK_ROOT = os.environ["WORK_ROOT"]
REG_DIR = os.environ.get("REGISTRY_V11_DIR",
                         os.path.join(WORK_ROOT, "00_handover",
                                      "handover_v11_data_ledger_2026-10-04"))
IN = os.path.join(REG_DIR, "V11F4_DATASET_hap_flags.tsv")
OUT = os.environ.get("REGISTRY_V2_OUT",
                     os.path.join(REG_DIR, "V11F4_DATASET_hap_flags_v2.tsv"))

# Items that are not excluded; kept here for the record only, and not put in `EXCL`
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

# ── 분모가 바뀌는지 — `제외` 를 넣은 뜻이 여기 있다 ──────────────────
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
