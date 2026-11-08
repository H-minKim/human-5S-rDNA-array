# -*- coding: utf-8 -*-
"""
EXPORT_origin_finalize, finalize the inherited/recurrent calls.
  (1) Add an effect size threshold: a call needs not only p<0.05 but also (top-lineage share - null) >= 0.10.
      Measured: the median effect among significant calls is +0.500 (1.4% below 0.10), among non-significant +0.013,
      so most calls survive the threshold, but tiny differences in variants with very many carriers are filtered out.
  (2) Compare with the lineage test values already produced by S3 to confirm the implementation.
  (3) Representative examples are drawn near the median effect size (drawing the most common variants would mislead).
Input : export\V11_ORIGIN_CALL.tsv (output of EXPORT_origin_call.py), primary\S3_lineage_test_results.tsv
Output : updates the same two files (origin_call recomputed)
"""
import io, os, collections
import numpy as np
import v11_io as G

OUTD = os.path.join(G.OUT, "export")
P = os.path.join(OUTD, "V11_ORIGIN_CALL.tsv")
EFF_MIN = 0.10

hd = None; R = []
for L in io.open(P, encoding="utf8"):
    f = L.rstrip("\n").split("\t")
    if hd is None: hd = f; H = {c: i for i, c in enumerate(hd)}; continue
    R.append(f)
def F(r, c):
    v = r[H[c]]
    try: return float(v)
    except: return float("nan")

# ---------- (2) S3 comparison ----------
S3 = {}
p3 = os.path.join(G.OUT, "primary", "S3_lineage_test_results.tsv")
if os.path.exists(p3):
    h3 = None
    for L in io.open(p3, encoding="utf8"):
        f = L.rstrip("\n").split("\t")
        if h3 is None: h3 = f; H3 = {c: i for i, c in enumerate(h3)}; continue
        try: S3[f[0]] = (float(f[H3["top_lineage_share"]]), float(f[H3["p_lineage"]]))
        except: pass
mine = {r[0]: (F(r, "top_lineage_share"), F(r, "p_lineage")) for r in R}
com = [x for x in S3 if x in mine and S3[x][0] == S3[x][0] and mine[x][0] == mine[x][0]]
print("=== S3 계통 검정값과 대조 ===")
if com:
    a = np.array([S3[x][0] for x in com]); b = np.array([mine[x][0] for x in com])
    d = np.abs(a - b)
    print("  공통 변이 %d · 최다계통 점유 차이 : 중앙 %.4f · 최대 %.4f · 0.05 이내 %.1f%%"
          % (len(com), np.median(d), d.max(), 100 * (d <= 0.05).mean()))
    sa = np.array([S3[x][1] < 0.05 for x in com]); sb = np.array([mine[x][1] < 0.05 for x in com])
    print("  유의 판정 일치 %.1f%% (S3 유의 %d · 이번 %d)" % (100 * (sa == sb).mean(), sa.sum(), sb.sum()))
    print("  ※ S3 는 순열 300 · 여기는 500 이라 p 가 완전히 같지는 않다(점유율은 같아야 한다)")
else:
    print("  공통 변이 없음")

# ---------- (1) recompute calls ----------
cnt_old = collections.Counter(r[H["origin_call"]] for r in R)
for r in R:
    top = F(r, "top_lineage_share"); nul = F(r, "top_lineage_null"); p = F(r, "p_lineage")
    per = F(r, "units_per_carrier_array"); hf = F(r, "hap_freq")
    rec = r[H["recurrence_context"]] == "Y"
    if hf >= 0.90:
        call = "조상형"
    elif top != top:
        call = "판정불가"
    else:
        eff = top - nul
        if p < 0.05 and eff >= EFF_MIN: call = "계승형"
        elif p < 0.05: call = "판정불가(효과작음)"
        elif p >= 0.20 and per <= 1.5 and rec: call = "재발시사형"
        elif p >= 0.20: call = "분산형"
        else: call = "판정불가"
    r[H["origin_call"]] = call
g = io.open(P, "w", encoding="utf8", newline="\n")
g.write("\t".join(hd) + "\n")
for r in R: g.write("\t".join(r) + "\n")
g.close()

cnt = collections.Counter(r[H["origin_call"]] for r in R)
print("")
print("=== 판정 (효과크기 문턱 %.2f 적용 후) ===" % EFF_MIN)
print("  판정            변이수   효과 중앙   배열당 사본 중앙   보유 하플로 중앙   재발문맥")
for c, n in cnt.most_common():
    sub = [r for r in R if r[H["origin_call"]] == c]
    eff = [F(r, "top_lineage_share") - F(r, "top_lineage_null") for r in sub]
    eff = [e for e in eff if e == e]
    print("  %-16s %5d   %s   %6.2f            %5.0f          %5.1f%%" % (
        c, n, ("%+.3f" % np.median(eff)) if eff else "   -  ",
        np.median([F(r, "units_per_carrier_array") for r in sub]),
        np.median([F(r, "n_haplotypes") for r in sub]),
        100 * np.mean([r[H["recurrence_context"]] == "Y" for r in sub])))
print("  (문턱 적용 전: %s)" % dict(cnt_old))

# ---------- (3) representative examples: near the median effect size ----------
print("")
print("=== 판정별 대표 예시 (효과크기 중앙 근처에서 뽑음) ===")
for c in ("계승형", "재발시사형", "분산형", "조상형"):
    sub = [r for r in R if r[H["origin_call"]] == c]
    if not sub: continue
    if c == "조상형":
        sub = sorted(sub, key=lambda r: -F(r, "n_units"))[:3]
    else:
        eff = np.array([F(r, "top_lineage_share") - F(r, "top_lineage_null") for r in sub])
        med = np.nanmedian(eff)
        sub = [sub[i] for i in np.argsort(np.abs(eff - med))[:3]]
    for r in sub:
        print("  %-12s %-12s 유닛 %6s · 배열 %4s · 배열당 %5s · 최다계통 %s (귀무 %s) p=%s · 재발문맥 %s" % (
            c, r[0], r[H["n_units"]], r[H["n_haplotypes"]], r[H["units_per_carrier_array"]],
            r[H["top_lineage_share"]] or "-", r[H["top_lineage_null"]] or "-",
            r[H["p_lineage"]] or "-", r[H["recurrence_context"]]))

# ---------- update group table columns ----------
src = os.path.join(OUTD, "V11_VARIANT_GROUPS.tsv")
bak = src + ".noorigin"
D = {r[0]: r for r in R}
hd2 = None; out = []
for L in io.open(bak, encoding="utf8"):
    f = L.rstrip("\n").split("\t")
    if hd2 is None:
        hd2 = f
        out.append(f + ["origin_call", "p_lineage", "top_lineage_share", "top_lineage_null",
                        "units_per_carrier_array", "recurrence_context"])
        continue
    r = D.get(f[0])
    if r is None:
        nu_ = float(f[hd2.index("n_units")]); nh_ = max(float(f[hd2.index("n_haplotypes")]), 1)
        out.append(f + ["판정불가(보유 하플로 부족)", "", "", "", "%.2f" % (nu_ / nh_), ""])
    else:
        out.append(f + [r[H["origin_call"]], r[H["p_lineage"]], r[H["top_lineage_share"]],
                        r[H["top_lineage_null"]], r[H["units_per_carrier_array"]],
                        r[H["recurrence_context"]]])
g = io.open(src, "w", encoding="utf8", newline="\n")
for f in out: g.write("\t".join(f) + "\n")
g.close()
print("")
print("그룹표 갱신:", src)
