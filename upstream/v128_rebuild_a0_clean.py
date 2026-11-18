# -*- coding: utf-8 -*-
# Rebuild FLANK_METH_HAP_A0_CLEAN.tsv, recomputing two samples (HG02717, HG02886) from the SAM copy whose MD5 matches the recorded one.
#   The copy of {HG02717,HG02886}.placed50k.final.sam under $V_ROOT/13_junction50k/final_SAM/node3/rj50k_final2/ does not match the MD5 recorded in MD5_13.txt;
#   the copy in 04_placement/junction50k_v2/node3/rj50k_final2/ matches it, so only those two samples are recomputed from it.
#   The method: the same computation as the a0 = ref_front version is run on those two files alone, and for the rest the per-file intermediate results of that version, in 04_logs/v118_partial, are used as they are,
#          then the final table is rebuilt the same way, by a groupby sum and a mean of beta.
#   Outputs: 02_tables/FLANK_METH_HAP_A0_CLEAN.tsv and 04_logs/V128_REBUILD_CLEAN.txt
import io, os, glob, hashlib
import numpy as np, pandas as pd
OUT = __import__("os").environ["RDNA_DATA_ROOT"] + "/v11_readvar"; T = OUT + "/02_tables"
OFF = __import__("os").environ["RDNA_DATA_ROOT"] + "/for_paper/08_versions/v11/13_junction50k/mapping_table/FLANK_OFFSET.tsv"
CLEAN = OUT + "/04_placement/junction50k_v2/node3/rj50k_final2"
PART = OUT + "/04_logs/v118_partial"
LOG = io.open(OUT + "/04_logs/V128_REBUILD_CLEAN.txt", "w", encoding="utf8", newline="\n")
def P(s=""): print(s, flush=True); LOG.write(str(s) + "\n"); LOG.flush()
md5 = lambda f: hashlib.md5(open(f, "rb").read()).hexdigest()
WINS = [(0, 2000), (2000, 5000), (5000, 10000), (10000, 20000), (20000, 50000)]
COMP = str.maketrans("ACGTacgtN", "TGCAtgcaN")
BAD = {"HG02717": "92ce16e6df0d9c44f3f8bd8a9cf3ba4c", "HG02886": "caf4665861b72e2c911c31a369fb403e"}

O = pd.read_csv(OFF, sep="\t", dtype={"slot": str})
REF = {r.ref_name: (r.sample, r.slot, int(r.ref_front), int(r.ref_front) + int(r.canon_len)) for r in O.itertuples(index=False)}

def ref_span(cig):
    n = 0; num = ""
    for ch in cig:
        if ch.isdigit(): num += ch
        elif ch in "MDN=X": n += int(num); num = ""
        else: num = ""
    return n

def qr_pairs(cig, pos):
    out = []; q = 0; r = pos; num = ""
    for ch in cig:
        if ch.isdigit(): num += ch; continue
        n = int(num); num = ""
        if ch in "M=X": out.append((q, r, n)); q += n; r += n
        elif ch in "IS": q += n
        elif ch in "DN": r += n
        elif ch == "H": pass
    return out

def read_calls(flag, seq, cig, pos, mm, ml):
    rev = bool(flag & 16)
    orig = seq.translate(COMP)[::-1] if rev else seq
    probs_off = 0; target = None
    for part in mm.rstrip(";").split(";"):
        f = part.split(",")
        nskip = len(f) - 1
        if f[0].startswith("C+m"):
            target = [int(x) for x in f[1:]]; break
        probs_off += nskip
    if target is None: return []
    cpos = [i for i, ch in enumerate(orig) if ch == "C"]
    idx = []; k = -1
    for s in target:
        k += s + 1
        if k >= len(cpos): break
        idx.append(cpos[k])
    pr = ml[probs_off:probs_off + len(idx)]
    if len(pr) < len(idx): idx = idx[:len(pr)]
    L = len(seq)
    segs = qr_pairs(cig, pos)
    qi = [(L - 1 - i) if rev else i for i in idx]
    order = np.argsort(qi); qi = np.array(qi)[order]; pv = np.array(pr, float)[order]
    res = []; j = 0
    for q0, r0, n in segs:
        while j < len(qi) and qi[j] < q0: j += 1
        while j < len(qi) and qi[j] < q0 + n:
            res.append((r0 + (qi[j] - q0), pv[j] / 255.0)); j += 1
    return res

def one(f):
    acc = {}
    with io.open(f, encoding="utf8", errors="replace") as fh:
        for line in fh:
            if line[0] == "@": continue
            t = line.rstrip("\n").split("\t")
            if len(t) < 11: continue
            rn = t[2]
            if rn not in REF: continue
            flag = int(t[1])
            if flag & 0x900: continue
            pos = int(t[3]) - 1; cig = t[5]; seq = t[9]
            if cig == "*" or seq == "*": continue
            a0, a1 = REF[rn][2], REF[rn][3]
            r0 = pos; r1 = pos + ref_span(cig)
            if not any((a0 - hi <= r1 and r0 <= a0 - lo) or (a1 + lo <= r1 and r0 <= a1 + hi) for lo, hi in WINS): continue
            mm = ml = None
            for x in t[11:]:
                if x.startswith("MM:Z:"): mm = x[5:]
                elif x.startswith("ML:B:C,"):
                    ml = [int(v) for v in x[7:].split(",") if v != ""]
                    if not ml: ml = None
            if mm is None or ml is None: continue
            for rp, p in read_calls(flag, seq, cig, pos, mm, ml):
                if rp < a0: d = a0 - rp; side = "원위"
                elif rp >= a1: d = rp - a1 + 1; side = "근위"
                else: continue
                for lo, hi in WINS:
                    if lo <= d < hi:
                        k = (REF[rn][0], REF[rn][1], side, "%d-%dk" % (lo // 1000, hi // 1000))
                        a = acc.setdefault(k, [0, 0.0, set()])
                        a[0] += 1; a[1] += p; a[2].add(t[0]); break
    return pd.DataFrame([dict(sample=k[0], slot=k[1], side=k[2], win=k[3], n_call=v[0], beta=v[1] / max(v[0], 1), n_read=len(v[2]))
                         for k, v in acc.items()])

P("[0] verify the clean copy")
rows = []
for smp, want in BAD.items():
    f = os.path.join(CLEAN, "%s.placed50k.final.sam" % smp)
    got = md5(f)
    P("   %s : %s (recorded %s) -> %s" % (smp, got[:8], want[:8], "일치" if got == want else "불일치"))
    assert got == want, smp
    d = one(f); d["src"] = smp
    rows.append(d)
    P("      recomputed : rows %d, distal 0-2k n_call %s" % (len(d), d[(d.side == "원위") & (d.win == "0-2k")].n_call.tolist()))
NEW = pd.concat(rows, ignore_index=True)

P()
P("[1] rebuild the table from the v118 intermediate results, with the corrupted part replaced by clean values")
parts = sorted(glob.glob(PART + "/*.tsv"))
P("   intermediate files %d" % len(parts))
keep = []
for p in parts:
    b = os.path.basename(p)
    if b.split("_")[0] in BAD:
        P("   excluded (based on the corrupted SAM) : %s" % b); continue
    keep.append(pd.read_csv(p, sep="\t", dtype={"slot": str}))
F = pd.concat(keep + [NEW.drop(columns=["src"])], ignore_index=True)
F = F.groupby(["sample", "slot", "side", "win"], as_index=False).agg(
    n_call=("n_call", "sum"), n_read=("n_read", "sum"), beta=("beta", "mean"))
F.to_csv(T + "/FLANK_METH_HAP_A0_CLEAN.tsv", sep="\t", index=False, float_format="%.4f")
P("   saved 02_tables/FLANK_METH_HAP_A0_CLEAN.tsv, rows %d, haplotypes %d" % (len(F), F.groupby(["sample", "slot"]).ngroups))

P()
P("[2] compare with the reference version (mirror_fix)")
M = pd.read_csv(__import__("os").environ["REGISTRY_DIR"] + "/verify/flank_orientation/mirror_fix/02_tables/FLANK_METH_HAP.tsv", sep="\t", dtype={"slot": str})
K = ["sample", "slot", "side", "win"]
J = F.merge(M, on=K, how="outer", suffixes=("_a", "_b"), indicator=True)
P("   rows: this version %d, reference version %d, both %d, one side only %d" % (len(F), len(M), (J._merge == "both").sum(), (J._merge != "both").sum()))
P("   β max difference %.8f, n_call max difference %s, n_read max difference %s"
  % ((J.beta_a - J.beta_b).abs().max(), (J.n_call_a - J.n_call_b).abs().max(), (J.n_read_a - J.n_read_b).abs().max()))
d = J[(J.beta_a - J.beta_b).abs() > 1e-9]
P("   differing rows %d %s" % (len(d), sorted(set(d["sample"])) if len(d) else ""))
P()
P("[3] difference from the v118 table built from the corrupted SAM (for reference)")
A = pd.read_csv(T + "/FLANK_METH_HAP_A0.tsv", sep="\t", dtype={"slot": str})
J2 = F.merge(A, on=K, suffixes=("_c", "_d"))
d2 = J2[(J2.beta_c - J2.beta_d).abs() > 1e-9]
P("   differing rows %d, samples %s" % (len(d2), sorted(set(d2["sample"]))))
LOG.close()
