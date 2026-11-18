# -*- coding: utf-8 -*-
# Locating a read fragment, by gene ends, then the unit chain, then the variant fingerprint, then dynamic programming over units. Aligning a 6-unit block with minimap2 fails and is not used.
#   Step 3b: canonical-unit repair (revision 1, 2 sites).

# -*- coding: utf-8 -*-
# Revision 1, stage 1: cut from each consensus read the fragment running from the start of the left anchor unit to the end of the right anchor unit, locating it by dynamic programming over units.
#
#   Aligning a 6-unit block to the read with minimap2 fails: the primary and secondary scores differ by 0.7%,
#                            and the left and right blocks take the same primary position on the read. The ONT error exceeds the between-unit variation, so base-level alignment cannot locate a position in a tandem array.
#   Instead the functions of uindel2.py, the same call that chose the consensus reads, are used as they are:
#         gene ends in the read by blastn, then the unit chain, then the per-unit variant fingerprint with asm20 --cs, then dynamic programming against the map, giving read unit against map array_pos.
#         The fragment runs from the start of the read unit paired with map unit aL to the end of the read unit paired with aR, in the read orientation given by the majority gene strand.
#   Taken unchanged from uindel2: gene_hits, gene_ends, chain_from, call_vars, score, dp and load_map, extracted as functions through ast.
#   Arguments: path to uindel2.py, individual, slot, left anchor, right anchor, read table, fastq.gz, output folder.
import ast, io, os, sys, gzip, random, socket, statistics, collections
src_path, smp, sl, aL, aR, rtab, fq, OUT = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]), sys.argv[6], sys.argv[7], sys.argv[8]
tree = ast.parse(io.open(src_path, encoding="utf8").read())
keep_nodes = []
for node in tree.body:
    if isinstance(node, (ast.Import, ast.ImportFrom, ast.FunctionDef)):
        keep_nodes.append(node)
    elif isinstance(node, ast.Assign):
        names = [t.id for t in node.targets if isinstance(t, ast.Name)] + [e.id for t in node.targets if isinstance(t, ast.Tuple) for e in t.elts if isinstance(e, ast.Name)]
        if any(n in ("TGT", "FQ", "NP", "ONLY", "OUTD") for n in names):
            continue
        keep_nodes.append(node)
mod = ast.Module(body=keep_nodes, type_ignores=[])
G = {"__name__": "uindel2_funcs"}
exec(compile(mod, src_path, "exec"), G)
NL, TAB, H = G["NL"], G["TAB"], G["H"]
# the two tables uindel2 fills at module level, FN and SUMF, are filled here by the same rule
for L in io.open(os.path.join(H, "hap_seq_node_v8.tsv"), encoding="utf8"):
    r = L.rstrip(NL).split(TAB)
    if len(r) >= 3 and r[0] != "sample": G["FN"][(r[0], r[1])] = r[2]
f = io.open(os.path.join(H, "g8canonv10", "CONSENSUS_SUM_G8_full.tsv"), encoding="utf8"); f.readline()
for L in f:
    r = L.rstrip(NL).split(TAB); G["SUMF"][(r[0], r[1])] = dict(N=int(r[3]), flip=(r[4] == "Y"), op=int(r[5]), od=int(r[6]))
os.makedirs(OUT, exist_ok=True)
W = os.path.join(OUT, "tmp_%s_%s" % (smp, sl)); os.makedirs(W, exist_ok=True)
rc, rdfa = G["rc"], G["rdfa"]
# the reads used for the consensus
use = {}
f = io.open(rtab, encoding="utf8"); h = f.readline().rstrip(NL).split(TAB)
ci = [i for i, n in enumerate(h) if n.startswith("use_for_consensus")][0]
for L in f:
    r = L.rstrip(NL).split(TAB)
    if len(r) > ci and r[ci] == "Y": use[r[h.index("read16")]] = r[h.index("status")]
SEQ = {}
with gzip.open(fq, "rt") as fh:
    while True:
        hd = fh.readline()
        if not hd: break
        s = fh.readline().strip(); fh.readline(); fh.readline()
        rid = hd[1:].split()[0]
        if rid[:16] in use and rid[:16] not in SEQ: SEQ[rid[:16]] = (rid, s.upper())
seqs = dict((k, v[1]) for k, v in SEQ.items())
HT = G["gene_hits"](seqs, W, "a")
flip = {k for k, hs in HT.items() if sum(1 for x in hs if x[0] == "-") > sum(1 for x in hs if x[0] == "+")}
if flip:
    HT2 = G["gene_hits"](dict((k, rc(seqs[k])) for k in flip), W, "b")
    for k in flip: seqs[k] = rc(seqs[k]); HT[k] = HT2.get(k, [])
MAPS, FIX = G["load_map"](smp)
mv, aps = MAPS[sl]
CH = {k: G["chain_from"](G["gene_ends"](HT.get(k, []), -10 ** 9, 10 ** 9)) for k in seqs}
useqs, uidx = [], []
for k, ch in CH.items():
    for i, it in enumerate(ch):
        if it[0] == "U": useqs.append(seqs[k][max(0, it[1] - 1 - G["PAD"]):it[2] + G["PAD"]]); uidx.append((k, i))
VV = G["call_vars"](useqs, W, "u") if useqs else {}
FP = collections.defaultdict(dict)
for n_, (k, i) in enumerate(uidx): FP[k][i] = VV.get(n_, set()) - FIX
random.seed(20260916)
pool_r = [v for d in FP.values() for v in d.values() if v]; pool_m = [v for v in mv if isinstance(v, set) and v]
base = statistics.median([G["jac"](random.choice(pool_r), random.choice(pool_m)) for _ in range(2000)]) if pool_r and pool_m else 0.1
# the units in between: the length of the canonical window from aL to aR in the map
U = {}
f = io.open(os.path.join(H, "g8canonv10", "CONSENSUS_MAP_G8_full.tsv"), encoding="utf8"); ix = dict((n, i) for i, n in enumerate(f.readline().rstrip(NL).split(TAB)))
for L in f:
    r = L.rstrip(NL).split(TAB)
    if (r[0], r[1]) == (smp, sl) and r[ix["seg_type"]] == "UNIT": U[int(r[ix["array_pos"]])] = (int(r[ix["unit_start"]]), int(r[ix["unit_end"]]))
asm_win = U[aR][1] - U[aL][0] + 1
host = socket.gethostname()
T = collections.Counter(); LOG = []
outf = os.path.join(OUT, "seg_%s_%s_%d.%s.fa" % (smp, sl, aL, host))
with io.open(outf, "w", newline=NL) as out:
    for k in sorted(seqs):
        T["read"] += 1
        ch = CH[k]
        # dynamic programming per fragment as divided by the breaks, the same division as uindel2
        segs, cur = [], []
        for i, it in enumerate(ch):
            if it[0] == "B":
                if cur: segs.append(cur)
                cur = []
            else:
                cur.append(i)
        if cur: segs.append(cur)
        hitL = hitR = None; nseg = 0; evs = []
        for sg in segs:
            R = [FP[k].get(i, set()) if ch[i][0] == "U" else set() for i in sg]
            if sum(1 for i in sg if ch[i][0] == "U") < G["MINSEG"]: continue
            nseg += 1
            sc, ops = G["dp"](R, mv, base)
            evs += [e for e in G["events"](ops, aps, R)]
            for o in ops:
                if o[0] != "M": continue
                ci_ = sg[o[1]]; ap = aps[o[2]]
                if ch[ci_][0] != "U": continue
                if ap == aL: hitL = ch[ci_]
                if ap == aR: hitR = ch[ci_]
        if hitL is None or hitR is None:
            T["앵커 유닛이 DP 짝에 없음"] += 1
            LOG.append((k, "앵커 짝 없음 L=%s R=%s · 사건 %s" % (hitL is not None, hitR is not None, evs))); continue
        s0, e1 = hitL[1], hitR[2]                     # 1-based: the start of the left anchor unit and the end of the right anchor unit
        if e1 - s0 < 500:
            T["조각 짧음/순서"] += 1; LOG.append((k, "짧음")); continue
        seg = seqs[k][s0 - 1:e1]
        out.write(">%s|%s|0|0%s%s%s" % (SEQ[k][0], use.get(k, "?"), NL, seg, NL))
        T["조각"] += 1
        LOG.append((k, "조각 %d bp · 조립본 창 %d 대비 %+.2f 유닛 · 조각 안 read 유닛 %d · 사건 %s" % (
            len(seg), asm_win, (len(seg) - asm_win) / 2240.0,
            sum(1 for it in ch if it[0] == "U" and it[1] >= s0 and it[2] <= e1),
            ",".join("%s%d@%s" % (e[0], e[1], e[2]) for e in evs) or "-")))
io.open(os.path.join(OUT, "SEG_LOG_%s_%s.tsv" % (smp, sl)), "w", encoding="utf8", newline=NL).write(
    "read16\tresult\n" + "".join("%s\t%s\n" % r for r in LOG))
for fn in os.listdir(W):
    try: os.remove(os.path.join(W, fn))
    except OSError: pass
os.rmdir(W)
print("%s h%s anchors %d/%d, assembly window %d bp, for consensus %d, fastq found %d, fingerprint threshold %.3f" % (smp, sl, aL, aR, asm_win, len(use), len(seqs), base))
for kk in sorted(T): print("   %-26s %d" % (kk, T[kk]))
for r in LOG: print("   %s  %s" % r)
