# -*- coding: utf-8 -*-
# Builds FIG_data_table.tsv and FIG_data_chain.tsv, the input tables for Supplementary Figure 1E.
#
#
#   A further stage builds T1_final.tsv and T2_final.tsv.
#     The computation below reads only column 0 of those two files, the sample name.
#     The fourth column, done, is filled from node logs, and
#     done does not affect which rows are included: for x in sorted(T) takes all of them.
#     Without those logs the done column is N throughout.
#
#   Inputs:
#     v11  02_canonical/hap_seq_node_v8.tsv, 02_canonical/CANON_V10_PATCH.tsv
#          03_maps/g8canonv10/CONSENSUS_MAP_G8_full.tsv
#          09_tests/2026-09-29_tier_test/log\*goctrl*.log, *go45*.log
#     03_result_tables/tier1_span.tsv
#     FIG5_cut_test/rawdata/chainpick2.out, chainpick_nodes.out
import io, os, sys, glob, collections
sys.stdout.reconfigure(encoding="utf-8")

WORK_ROOT = os.environ["WORK_ROOT"]
V = os.environ.get("V_ROOT", os.path.join(WORK_ROOT, "08_versions", "v11"))
L = os.environ.get("TIER_LOG_DIR",
                   os.path.join(V, "09_tests", "2026-09-29_tier_test", "log"))
TABLE_DIR = os.environ.get("TABLE_DIR",
                           os.path.join(WORK_ROOT, "04_figures", "_tidy_2026-09-30"))
CHAIN_DIR = os.environ.get("CHAIN_DIR", os.path.join(TABLE_DIR, "FIG5_cut_test", "rawdata"))
RESULT_DIR = os.environ.get("RESULT_DIR", os.path.join(WORK_ROOT, "03_result_tables"))
OUTD = os.environ.get("TRIM_OUT_DIR", CHAIN_DIR)

# -- 1) the T1 and T2 lists; for the done column see the header above
info = collections.defaultdict(dict)
for _L in io.open(os.path.join(V, "02_canonical", "hap_seq_node_v8.tsv"), encoding="utf8"):
    r = _L.rstrip("\n").split("\t")
    if len(r) >= 3:
        p = r[2].split("__"); info[r[0]][r[1]] = (p[-1].split("_")[0], p[1], p[0])
rep = set()
for i, _L in enumerate(io.open(os.path.join(V, "02_canonical", "CANON_V10_PATCH.tsv"), encoding="utf8")):
    if i:
        r = _L.rstrip("\n").split("\t"); rep.add((r[0], r[1]))
span = collections.defaultdict(set)
for _L in io.open(os.path.join(RESULT_DIR, "tier1_span.tsv"), encoding="utf8"):
    r = _L.rstrip("\n").split("\t")
    if len(r) >= 2: span[r[0]].add(r[1])
T1 = []; T2 = []
for s, v in info.items():
    if len(v) != 2 or any((s, sl) in rep for sl in v) or len(span.get(s, set())) < 2: continue
    b = [v[sl][0] for sl in v]; nd = sorted(set(v[sl][1] for sl in v))
    st = "+".join(sorted(v[sl][2] for sl in v))
    if all(x == "hiarr" for x in b): T1.append((s, "+".join(nd), st))
    elif all(x == "arr" for x in b): T2.append((s, "+".join(nd), st))
print("T1 candidates %d, T2 candidates %d" % (len(T1), len(T2)))

# -- 2) --
c = collections.Counter()
f = io.open(os.path.join(V, "03_maps", "g8canonv10", "CONSENSUS_MAP_G8_full.tsv"), encoding="utf8")
h = f.readline().rstrip("\n").split("\t"); ix = {k: i for i, k in enumerate(h)}
for x in f:
    r = x.rstrip("\n").split("\t")
    if r[ix["seg_type"]] == "UNIT": c[(r[0], r[1])] += 1
sc = collections.defaultdict(dict)
for fn in ("chainpick2.out", "chainpick_nodes.out"):
    for x in io.open(os.path.join(CHAIN_DIR, fn), encoding="utf8"):
        p = [y.strip() for y in x.strip().strip("|").split("|")]
        if len(p) < 6 or not p[1].isdigit(): continue
        smp, a, b, v = p[0], int(p[1]), int(p[2]), float(p[3].replace("*", ""))
        sc[smp][(a, b)] = v; sc[smp][(b, a)] = v


def pg(pat):
    o = collections.defaultdict(set)
    for p in glob.glob(os.path.join(L, pat)):
        for x in io.open(p, encoding="utf8"):
            r = x.rstrip("\n").split("\t")
            t = next((y for y in r if y.endswith("_trim")), None)
            if t and r[-1].isdigit(): o[t.replace("_trim", "")].add(int(r[-1]))
    return o


CT = pg("*goctrl*.log"); T4 = pg("*go45*.log")
T = {"T1": set(x[0] for x in T1) | {"sample"}, "T2": set(x[0] for x in T2) | {"sample"}}
EX = {"HG04153", "HG00146", "HG01136", "HG02392", "HG00344"}
tier = {s: n for n in ("T1", "T2") for s in (T[n] - EX)}
elig = [s for s in sorted(set(CT) & set(T4)) if s in tier]
TH = 0.80


def assign(ou, g):
    av = list(g); real = {}
    for tol in (0, 1):
        for u in ou:
            if u in real: continue
            m = [x for x in av if abs(x - u) <= tol]
            if m: b = min(m, key=lambda z: abs(z - u)); av.remove(b); real[u] = b
    return real, av


def match(ou, pick):
    av = list(pick); done = set(); n = 0
    for tol in (0, 1):
        for u in ou:
            if u in done: continue
            m = [x for x in av if abs(x - u) <= tol]
            if m: av.remove(min(m, key=lambda z: abs(z - u))); done.add(u); n += 1
    return n


rows = []; pairs = []
for s in elig:
    ou = sorted([c[k] for k in c if k[0] == s], reverse=True)
    g = sorted(T4[s], reverse=True)
    pick = []
    for u in g:
        if len(pick) >= len(ou): break
        if any(sc.get(s, {}).get((u, p2), 0) >= TH for p2 in pick): continue
        pick.append(u)
    real, fake = assign(ou, g); rv = list(real.values())
    rows.append(dict(sample=s, tier=tier[s], n_hap=len(ou), canon=",".join(map(str, ou)),
                     cands=",".join(map(str, g)), n_cand=len(g), n_sel=len(pick),
                     n_ok=match(ou, pick), n_avail=match(ou, g), chain=1 if s in sc else 0))
    for i in range(len(rv)):
        for j in range(i + 1, len(rv)):
            v = sc.get(s, {}).get((rv[i], rv[j]))
            if v is not None: pairs.append((s, rv[i], rv[j], v, "must_separate"))
    for u in fake:
        vs = [(sc.get(s, {}).get((u, r2), -1), r2) for r2 in rv]; vs = [x for x in vs if x[0] >= 0]
        if vs: v, r2 = max(vs); pairs.append((s, r2, u, v, "must_collapse"))
with io.open(os.path.join(OUTD, "FIG_data_table.tsv"), "w", encoding="utf8", newline="\n") as g:
    ks = ["sample", "tier", "n_hap", "canon", "cands", "n_cand", "n_sel", "n_ok", "n_avail", "chain"]
    g.write("\t".join(ks) + "\n")
    for r in rows: g.write("\t".join(str(r[k]) for k in ks) + "\n")
with io.open(os.path.join(OUTD, "FIG_data_chain.tsv"), "w", encoding="utf8", newline="\n") as g:
    g.write("sample\tA\tB\tjaccard\tkind\n")
    for p in sorted(pairs, key=lambda z: z[3]): g.write("\t".join(map(str, p)) + "\n")
print("individuals %d, haplotypes %d, assemblies %d, selected %d, pairs %d, individuals with a chain %d"
      % (len(rows), sum(r["n_hap"] for r in rows), sum(r["n_avail"] for r in rows),
         sum(r["n_ok"] for r in rows), len(pairs), sum(r["chain"] for r in rows)))
