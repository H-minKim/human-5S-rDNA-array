# -*- coding: utf-8 -*-
# Fig S1 A 의 입력 표 `FIG2_data_unitcmp.tsv` 를 만드는 글
#
#   이 글은 2026-09-30 01:47:35 에 셸 heredoc 으로 한 번 돌았고 파일로 저장되지 않았다.
#   셈은 한 줄도 바꾸지 않았다 — 바꾼 것은 셋뿐이다 :
#     ⑴ heredoc 껍데기(`cd $S && python - <<'PY'`)를 벗겼다
#     ⑵ 입력·출력 자리를 환경변수로 뺐다 (원래는 스크래치 폴더의 상대 이름이었다)
#     ⑶ `# -*- coding -*-` 와 이 머리글을 더했다
#   원문 md5(기록에서 꺼낸 그대로) : ccd1b588014d
#
#   입력 둘은 `FIG6_canonical_test_table\rawdata\` 에 보존돼 있다 :
#     `genepos.tsv` · `pubspan2_pairs.tsv`
import io, os, sys, collections, statistics
sys.stdout.reconfigure(encoding="utf-8")

WORK_ROOT = os.environ["WORK_ROOT"]
TABLE_DIR = os.environ.get("TABLE_DIR",
                           os.path.join(WORK_ROOT, "04_figures", "_tidy_2026-09-30"))
S = os.environ.get("UNITCMP_DIR", os.path.join(TABLE_DIR, "FIG6_canonical_test_table", "rawdata"))
OUT = os.environ.get("UNITCMP_OUT", os.path.join(S, "FIG2_data_unitcmp.tsv"))

# ──────────────────────────────────────────────────────────────────────
# 검정 ① : 짝별 유닛 수 차 · 길이 차        아래는 원문 그대로
# ──────────────────────────────────────────────────────────────────────
G = {}
f = io.open(os.path.join(S, "genepos.tsv"), encoding="utf8")
h = f.readline().rstrip("\n").split("\t"); ix = {k: i for i, k in enumerate(h)}
for L in f:
    r = L.rstrip("\n").split("\t"); G[r[0]] = {k: r[ix[k]] for k in ix}
rows = [L.rstrip("\n").split("\t") for L in io.open(os.path.join(S, "pubspan2_pairs.tsv"), encoding="utf8")]
hh = rows[0]; jx = {k: i for i, k in enumerate(hh)}; rows = rows[1:]
U = 2240.0
out = []
for r in rows:
    A = G[r[jx["pubfile"]]]; B = G[r[jx["spanfile"]]]
    if int(A["n_gene"]) <= 0 or int(B["n_gene"]) <= 0: continue
    ua = int(int(A["total_bp"]) // U); ub = int(int(B["total_bp"]) // U)
    out.append((r[0], r[1], r[jx["src"]], int(r[jx["units"]]), ua, ub, ub - ua,
                round((int(B["total_bp"]) - int(A["total_bp"])) / U, 4)))
with io.open(OUT, "w", encoding="utf8", newline="\n") as g:
    g.write("sample\tslot\tsrc\tasm_units_geneN\tasm_units_len\tspan_units_len\td_unit\td_len_unit\n")
    for x in out: g.write("\t".join(map(str, x)) + "\n")
c = collections.Counter((x[2], x[6]) for x in out)
print("검정① 짝 %d (HPRC %d · hifiasm %d)" % (len(out), sum(1 for x in out if x[2] == "arr"), sum(1 for x in out if x[2] == "hiarr")))
for s in ("arr", "hiarr"):
    n = sum(v for k, v in c.items() if k[0] == s)
    print("   %s : 일치 %d/%d (%.1f%%)" % (s, c[(s, 0)], n, 100.0 * c[(s, 0)] / n))
d = [x[7] for x in out]
print("   길이차 중앙 %+.4f 유닛 (%+.0f bp) · |차|<0.5 %d/%d (%.1f%%)"
      % (statistics.median(d), statistics.median(d) * U, sum(1 for x in d if abs(x) < 0.5), len(d),
         100.0 * sum(1 for x in d if abs(x) < 0.5) / len(d)))
