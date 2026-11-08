# -*- coding: utf-8 -*-
# Fig S1 B 의 입력 표 `FIG2_data_hifi.tsv` 를 만드는 글
#
#   이 글도 2026-09-30 01:47:56 에 셸 heredoc 으로 돌았고 저장되지 않았다
#   (`FIG2_data_unitcmp` 보다 19 초 뒤 · 같은 작업의 둘째 토막이다).
#   바꾼 것 : heredoc 껍데기를 벗기고 입력·출력 자리를 환경변수로 뺐다.
#
#   입력 하나 : `03_result_tables\2026-09-07_map_test\JXHIFI_ALL.tsv.gz`
import io, os, sys, gzip, collections
sys.stdout.reconfigure(encoding="utf-8")

WORK_ROOT = os.environ["WORK_ROOT"]
D = os.environ.get("JXHIFI_DIR", os.path.join(WORK_ROOT, "03_result_tables", "2026-09-07_map_test"))
TABLE_DIR = os.environ.get("TABLE_DIR",
                           os.path.join(WORK_ROOT, "04_figures", "_tidy_2026-09-30"))
OUT = os.environ.get("HIFI_OUT",
                     os.path.join(TABLE_DIR, "FIG6_canonical_test_table", "rawdata", "FIG2_data_hifi.tsv"))

# ──────────────────────── 아래는 원문 그대로 ────────────────────────
LAB = {"asmunits3": "HPRC", "arr": "hifiasm", "fill60": "hifiasm+fill", "tier1ref": "span",
       "resub2": "resub2", "refix": "refix", "sweep31": "sweep31", "nopurge": "nopurge",
       "resub": "resub"}
SUB = {75: ("G", "A"), 369: ("A", "G"), 563: ("C", "G"), 568: ("G", "A"), 1096: ("T", "G"),
       1344: ("A", "C"), 1387: ("C", "T"), 1398: ("G", "C"), 1683: ("G", "T"), 1756: ("A", "C")}
rows = [l.rstrip("\n").split("\t") for l in
        gzip.open(os.path.join(D, "JXHIFI_ALL.tsv.gz"), "rt", encoding="utf8")]
rows = [r for r in rows if len(r) >= 8]


def toks(s):
    d = {}
    if s in ("-", ""): return d
    for x in s.split(";"):
        p = x.rsplit(":", 1)
        if len(p) == 2: d[p[0]] = int(p[1])
    return d


def flip(S):
    o = set(S)
    for p, (cb, mb) in SUB.items():
        cur = [t for t in S if len(t) > 1 and t[:-1] == str(p) and t[-1] in "ACGT"]
        al = cur[0][-1] if cur else cb
        for t in cur: o.discard(t)
        if al != mb: o.add("%d%s" % (p, al))
    if "437D1" in S: o.discard("437D1")
    elif not any(t.startswith("437D") or t.startswith("437I") for t in S): o.add("437REF")
    return o


res = {}
for ref, F in (("CHM13", lambda x: x), ("MAJ", flip)):
    BY = collections.defaultdict(lambda: [0, 0, 0]); SL = collections.defaultdict(set)
    for r in rows:
        if r[7] == "NA": continue
        d = int(r[4])
        if d < 5: continue
        T = toks(r[6]); J = F(set(k for k, v in T.items() if v / float(d) >= 0.5))
        M = F(set() if r[7] == "-" else set(r[7].split(",")))
        k = LAB.get(r[5], r[5]); SL[k].add((r[0], r[1]))
        BY[k][0] += len(J & M); BY[k][1] += len(J - M); BY[k][2] += len(M - J)
    res[ref] = (BY, SL)
with io.open(OUT, "w", encoding="utf8", newline="\n") as g:
    g.write("source\tslots\tref\trecall\tprecision\tn_map_var\n")
    B0, S0 = res["CHM13"]; B1, _ = res["MAJ"]
    for k in sorted(B0, key=lambda x: -(B0[x][0] + B0[x][2])):
        c, fp, fn = B0[k]
        if c + fn < 200: continue
        for ref, BB in (("CHM13", B0), ("MAJ", B1)):
            c2, f2, n2 = BB[k]
            g.write("%s\t%d\t%s\t%.4f\t%.4f\t%d\n"
                    % (k, len(S0[k]), ref, 100.0 * c2 / max(1, c2 + n2),
                       100.0 * c2 / max(1, c2 + f2), c2 + n2))
print(io.open(OUT, encoding="utf8").read())
