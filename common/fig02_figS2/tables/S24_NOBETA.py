# -*- coding: utf-8 -*-
r"""Can i and j be reconciled? Reruns j (S24) without methylation beta.

   Why: i (S18) uses only variant features (20 characters / 160 individual) and has no beta.
        j (S24) uses 20 variant features + beta + total variants per unit = 22 features,
        and the top feature of distal zone 1 is `beta_meth-1.4`, so beta may be pulling that boundary.
   => Run the same procedure without beta and see whether the distal boundary moves toward the mode of i (7).

   The original code is not modified. A copy is made and run, and outputs stay only in my folder."""
import io, os, sys, re, shutil, subprocess, time
sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = __import__("os").environ.get("V11_0913_CODE", __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
TMP = os.path.join(HERE, "_s24nb")
REPRO = os.path.join(HERE, "REPRO")
if os.path.exists(TMP): shutil.rmtree(TMP)
os.makedirs(TMP)
for f in os.listdir(SRC):
    if f.endswith(".py"): shutil.copyfile(os.path.join(SRC, f), os.path.join(TMP, f))

p = os.path.join(TMP, "S24_optimal_zoning.py")
t = io.open(p, encoding="utf8").read()
a = '''        ds.append(d); X.append(np.concatenate([v, [b]]))'''
assert t.count(a) == 1
t = t.replace(a, '''        ds.append(d)
        X.append(v if NOBETA else np.concatenate([v, [b]]))''')
a2 = '''NF = len(FB) + 1                      # + beta'''
assert t.count(a2) == 1
t = t.replace(a2, '''NOBETA = os.environ.get("S24_NOBETA", "") == "1"   # β 를 빼고 돌릴 때
NF = len(FB) + (0 if NOBETA else 1)   # + beta''')
t = t.replace('names = FB + ["beta_meth"]', 'names = FB + ([] if NOBETA else ["beta_meth"])')
io.open(p, "w", encoding="utf8", newline="\n").write(t)
print("사본 고침 :", p)

for tag, nb in (("withbeta", "0"), ("nobeta", "1")):
    out = os.path.join(HERE, "_s24_" + tag)
    if os.path.exists(out): shutil.rmtree(out)
    os.makedirs(out)
    e = {**os.environ, "PYTHONIOENCODING": "utf-8", "PYTHONHASHSEED": "0",
         "S24_NOBETA": nb,
         "REPRO_MIRROR": out, "V11_ROOT": os.path.join(REPRO, "mirror_h1f"), "V11_TIER": "h1",
         "PYTHONPATH": os.path.join(HERE, "hook_tier") + os.pathsep + os.environ.get("PYTHONPATH", "")}
    t0 = time.time()
    r = subprocess.run([sys.executable, p, "primary"], cwd=TMP, env=e,
                       capture_output=True, text=True, encoding="utf8", errors="replace", timeout=5400)
    print("\n■ %-9s %s %.1fs" % (tag, "ok" if r.returncode == 0 else "실패", time.time() - t0))
    if r.returncode:
        print((r.stderr or "")[-1500:]); continue
    for ln in (r.stdout or "").split("\n"):
        if ln.strip(): print("   ", ln.strip()[:150])
    f = os.path.join(out, "primary", "S24_optimal_zoning.tsv")
    if os.path.exists(f):
        L = [l.rstrip("\n").split("\t") for l in io.open(f, encoding="utf8") if l.strip() and not l.startswith("#")]
        H = {c: i for i, c in enumerate(L[0])}
        for side in ("prox", "dist"):
            z = sorted([x for x in L[1:] if x[H["side"]] == side and x[H["pen"]] == "1.0"],
                       key=lambda x: int(x[H["zone"]]))
            print("   %-5s 벌점 1.0 : %s" % (side, " | ".join("%s-%s" % (x[H["depth_from"]], x[H["depth_to"]]) for x in z)))
            for x in z: print("        %s" % x[H["top_features"]][:80])
