# -*- coding: utf-8 -*-
# GM12878(HG001) hifiasm hap1, hap2 array -> unit split, orientation, flank numbering, variants vs canon_v2
#   hiarr method: first check whether both 50kb flanks are on the same contig (if not, the number of that end is NA)
#   unit boundary = 5S gene anchor (gene at the unit end), unit 1 = distal flank side (same axis as K562 hapA)
#   output ~/p11_encode/gm12878/units/GM_UNITS.tsv : hap, unit(distal 1), prox_idx, variants
import subprocess, re, os
W = os.path.expanduser("~/p11_encode"); O = W + "/gm12878/units"; os.makedirs(O, exist_ok=True)
C = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
rc = lambda s: "".join(C.get(c, "N") for c in reversed(s))
G = "GTCTACGGCCATACCACCCTGAACGCGCCCGATCTCGTCTGATCTCGGAAGCTAAGCAGGGTCGGGCCTGGTTAGTACTTGGATGGGAGACCGCCTGGGAATACCGGGTGCTGTAGGCTTT"[:119]
def scan(s, g):
    cand = set()
    for s0 in range(0, 108, 12):
        sd = g[s0:s0 + 12]; j = s.find(sd)
        while j >= 0: cand.add(j - s0); j = s.find(sd, j + 1)
    return sorted(p for p in cand if 0 <= p and p + 119 <= len(s) and sum(1 for x, y in zip(s[p:p + 119], g) if x != y) <= 12)
def rdfa(p):
    d = {}; n = None
    for L in open(p):
        if L.startswith(">"): n = L[1:].split()[0]; d[n] = []
        else: d[n].append(L.strip())
    return {k: "".join(v) for k, v in d.items()}
out = open(O + "/GM_UNITS.tsv", "w"); out.write("hap\tcontig\tunit\tprox_idx\tn_units\tvariants\n")
for hap in ("hap1", "hap2"):
    ctg = rdfa(W + "/gm12878/hifi/asm/HG001.bp.%s.p_ctg.fa" % hap)
    # contig with the most 5S gene anchors
    best = max(ctg, key=lambda k: len(scan(ctg[k], G)) + len(scan(ctg[k], rc(G))))
    s = ctg[best]
    if len(scan(s, rc(G))) > len(scan(s, G)): s = rc(s)          # so that the gene is on the + strand
    open(O + "/%s.fa" % hap, "w").write(">%s\n%s\n" % (hap, s))
    fl = {}
    for side in ("distal", "proximal"):
        for L in subprocess.run(["minimap2", "-c", "-x", "asm10", W + "/ref/FLANK_%s_50000.fa" % side, O + "/%s.fa" % hap], capture_output=True, text=True).stdout.splitlines():
            f = L.split("\t")
            if int(f[10]) >= 20000 and (side not in fl or int(f[10]) > fl[side][2]): fl[side] = ((int(f[2]) + int(f[3])) / 2, f[4], int(f[10]))
    P = scan(s, G); us = [(P[k - 1] + 119 if k > 0 else max(0, P[0] + 119 - 2240), P[k] + 119) for k in range(len(P))]   # same rule as K562 asmunits: include the end unit that ends with the first gene
    # unit 1 = distal flank side
    if "distal" in fl: dfront = fl["distal"][0] < us[0][0]
    elif "proximal" in fl: dfront = not (fl["proximal"][0] < us[0][0])
    else: dfront = True
    order = list(range(len(us))) if dfront else list(range(len(us)))[::-1]
    open(O + "/u.fa", "w").write("".join(">u%d\n%s\n" % (i, s[a:b]) for i, (a, b) in enumerate(us)))
    V = {}
    for L in subprocess.run(["minimap2", "-c", "--cs=short", "-x", "asm20", W + "/ref/canon_v2.fa", O + "/u.fa"], capture_output=True, text=True).stdout.splitlines():
        f = L.split("\t"); i = int(f[0][1:])
        if i in V: continue
        cs = [x for x in f[12:] if x.startswith("cs:Z:")][0][5:]; t = int(f[7]); vs = []
        for op, val in re.findall(r"([:*+\-])([0-9a-z]+)", cs):
            if op == ":": t += int(val)
            elif op == "*": vs.append("%d%s%s" % (t + 1, val[0].upper(), val[1].upper())); t += 1
            elif op == "+": vs.append("I%d:%d" % (t, len(val)))
            elif op == "-": vs.append("D%d:%d" % (t + 1, len(val))); t += len(val)
        V[i] = vs
    n = len(us)
    for k, i in enumerate(order):
        out.write("%s\t%s\t%d\t%d\t%d\t%s\n" % (hap, best, k + 1, n - k, n, ",".join(V.get(i, [])) or "-"))
    print(hap, best, "길이", len(s), "유닛", n, "flank", {k: (int(v[0]), v[1], v[2]) for k, v in fl.items()}, "unit1=원위쪽" if dfront else "unit1=서열 뒤")
out.close()
