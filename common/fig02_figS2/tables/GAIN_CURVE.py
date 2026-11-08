# -*- coding: utf-8 -*-
r"""Put "why cut there" into numbers, output the split gain curve and the penalty threshold as a table.

   Decision rule of S24 (as in the code):
     - the gain of cutting at position t is g(t) = SSE(whole) - SSE(left) - SSE(right)
     - pick the largest g and cut only when it exceeds the penalty pen = PEN x number of features x ln(DMAX)
     - after cutting, repeat the same within each piece (binary segmentation)
   => with this table the figure can show that "3|4 and 9|10 are the only positions whose gain exceeds the threshold".

   beta is excluded (unified 2026-10-09), outputs go only to my folder."""
import io, os, sys, subprocess, time, shutil
sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
TMP = os.path.join(HERE, "_s24nb")
MIR = os.path.join(HERE, "REPRO", "mirror_h1f")
OUT = os.path.join(HERE, "_gain")
os.makedirs(OUT, exist_ok=True)

# Run only the head of the S24 copy (feature construction), then compute the gain curve directly
CODE = r'''
import io, os, sys, numpy as np
sys.argv = ["S24", "primary"]
src = io.open(r"%s", encoding="utf8").read()
head = src.split('g = io.open(os.path.join(OUT, "S24_optimal_zoning.tsv")')[0]
exec(compile(head, "s24head", "exec"))
import math
rows = []
for side in ("prox", "dist"):
    ds, X = profile(side)
    Z = zscore(X)
    nf = Z.shape[1]
    for PEN in (0.5, 1.0, 2.0):
        pen = PEN * nf * math.log(DMAX)
        # 되풀이하며 각 단계의 이득 곡선을 모두 적는다
        segs = [(0, len(Z))]
        step = 0
        while True:
            cand = None
            curves = []
            for si, (a, b) in enumerate(segs):
                S_ = Z[a:b]; m = len(S_)
                if m < 2 * MINSEG: continue
                cs1 = np.cumsum(S_, 0); cs2 = np.cumsum(S_ * S_, 0); t1 = cs1[-1]; t2 = cs2[-1]
                mu = t1 / m; c0 = float((t2 - 2 * mu * t1 + m * mu * mu).sum())
                for t in range(MINSEG, m - MINSEG + 1):
                    l1 = cs1[t - 1]; l2 = cs2[t - 1]
                    cl = float((l2 - l1 * l1 / t).sum())
                    r1 = t1 - l1; r2 = t2 - l2; rr = m - t
                    cr = float((r2 - r1 * r1 / rr).sum())
                    gph = c0 - (cl + cr)
                    cut_depth = int(ds[a + t - 1])     # 이 깊이와 다음 깊이 사이를 자른다
                    curves.append((si, cut_depth, gph))
                    if cand is None or gph > cand[0]: cand = (gph, si, t, cut_depth)
            for si, cd, gph in curves:
                rows.append((side, "%%.1f" %% PEN, step, si, cd, "%%.4f" %% gph, "%%.4f" %% pen,
                             "1" if (cand and cd == cand[3] and si == cand[1]) else "0"))
            if cand is None or cand[0] <= pen: break
            g_, si, t, cd = cand
            a, b = segs.pop(si); segs.append((a, a + t)); segs.append((a + t, b)); segs.sort()
            step += 1
            if step > 8: break
g = io.open(r"%s", "w", encoding="utf8", newline="\n")
g.write("side\tpen\tstep\tseg\tcut_after_depth\tgain\tthreshold\tchosen\n")
for r in rows: g.write("\t".join(str(x) for x in r) + "\n")
g.close()
print("이득 곡선 줄 %%d · 특징 수 %%d" %% (len(rows), nf))
''' % (os.path.join(TMP, "S24_optimal_zoning.py").replace("\\", "/"),
       os.path.join(OUT, "S24_gain_curve.tsv").replace("\\", "/"))

run = os.path.join(TMP, "_gain_run.py")
io.open(run, "w", encoding="utf8", newline="\n").write(CODE)
e = {**os.environ, "PYTHONIOENCODING": "utf-8", "PYTHONHASHSEED": "0", "S24_NOBETA": "1",
     "REPRO_MIRROR": MIR, "V11_ROOT": MIR, "V11_TIER": "h1",
     "PYTHONPATH": os.path.join(HERE, "hook_tier") + os.pathsep + os.environ.get("PYTHONPATH", "")}
t0 = time.time()
r = subprocess.run([sys.executable, run], cwd=TMP, env=e, capture_output=True, text=True,
                   encoding="utf8", errors="replace", timeout=5400)
print("이득 곡선 %s %.1fs" % ("ok" if r.returncode == 0 else "실패", time.time() - t0))
print((r.stdout or "").strip()[-400:])
if r.returncode: print((r.stderr or "")[-1500:]); sys.exit(1)

f = os.path.join(OUT, "S24_gain_curve.tsv")
L = [l.rstrip("\n").split("\t") for l in io.open(f, encoding="utf8") if l.strip()]
H = {c: i for i, c in enumerate(L[0])}
print("\n■ 벌점 1.0 에서 채택된 자리와 그 이득")
for side in ("prox", "dist"):
    rs = [r for r in L[1:] if r[H["side"]] == side and r[H["pen"]] == "1.0"]
    th = float(rs[0][H["threshold"]])
    print("   %-5s 문턱 %.1f" % (side, th))
    for r in rs:
        if r[H["chosen"]] == "1":
            print("      단계 %s · 깊이 %s 뒤에서 자름 · 이득 %.1f %s"
                  % (r[H["step"]], r[H["cut_after_depth"]], float(r[H["gain"]]),
                     "(문턱 넘음 ⇒ 채택)" if float(r[H["gain"]]) > th else "(문턱 못 넘음 ⇒ 멈춤)"))
