# -*- coding: utf-8 -*-
r"""Runs S24 separately for the version with beta and the version without beta.
   Existing outputs in the mirror are backed up and restored (no overwriting)."""
import io, os, sys, shutil, subprocess, time, hashlib
sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
TMP = os.path.join(HERE, "_s24nb")
MIR = os.path.join(HERE, "REPRO", "mirror_h1f")
OUTF = os.path.join(MIR, "primary", "S24_optimal_zoning.tsv")
md5 = lambda q: hashlib.md5(open(q, "rb").read()).hexdigest()
bak = None
if os.path.exists(OUTF):
    bak = OUTF + ".bak_s24run"
    shutil.copy2(OUTF, bak)
    print("■ 거울 기존 산출 백업 ·", md5(OUTF)[:10])

def zones(f, pen="1.0"):
    L = [l.rstrip("\n").split("\t") for l in io.open(f, encoding="utf8") if l.strip() and not l.startswith("#")]
    H = {c: i for i, c in enumerate(L[0])}
    out = {}
    for side in ("prox", "dist"):
        z = sorted([x for x in L[1:] if x[H["side"]] == side and x[H["pen"]] == pen], key=lambda x: int(x[H["zone"]]))
        out[side] = [("%s-%s" % (x[H["depth_from"]], x[H["depth_to"]]), x[H["top_features"]]) for x in z]
    return out

res = {}
for tag, nb in (("withbeta", "0"), ("nobeta", "1")):
    e = {**os.environ, "PYTHONIOENCODING": "utf-8", "PYTHONHASHSEED": "0", "S24_NOBETA": nb,
         "REPRO_MIRROR": MIR, "V11_ROOT": MIR, "V11_TIER": "h1",
         "PYTHONPATH": os.path.join(HERE, "hook_tier") + os.pathsep + os.environ.get("PYTHONPATH", "")}
    t0 = time.time()
    r = subprocess.run([sys.executable, os.path.join(TMP, "S24_optimal_zoning.py"), "primary"],
                       cwd=TMP, env=e, capture_output=True, text=True, encoding="utf8", errors="replace", timeout=5400)
    print("\n■ %-9s %s %.1fs" % (tag, "ok" if r.returncode == 0 else "실패", time.time() - t0))
    if r.returncode:
        print((r.stderr or "")[-1200:]); continue
    for ln in (r.stdout or "").split("\n"):
        if "가운데 기준" in ln or "구역" in ln: print("   ", ln.strip()[:140])
    d = os.path.join(HERE, "_s24_" + tag)
    os.makedirs(d, exist_ok=True)
    shutil.copy2(OUTF, os.path.join(d, "S24_optimal_zoning.tsv"))
    res[tag] = zones(os.path.join(d, "S24_optimal_zoning.tsv"))

if bak:
    shutil.copy2(bak, OUTF); os.remove(bak)
    print("\n■ 거울 산출 되돌림 ·", md5(OUTF)[:10])

print("\n■ 결과 (벌점 1.0)")
for tag in ("withbeta", "nobeta"):
    if tag not in res: continue
    for side in ("prox", "dist"):
        print("   %-9s %-5s : %s" % (tag, side, " | ".join(z for z, _ in res[tag][side])))
        for z, f in res[tag][side]: print("        %-9s %s" % (z, f[:76]))
