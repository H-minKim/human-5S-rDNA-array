# -*- coding: utf-8 -*-
# Extract new map rows from sequence consensus 2 (so that the placement reference also switches to consensus 2)
#
#  Input : argv[1] = list of "sample slot",  ~/c2seq/<sample>_<slot>/cons2.fa
#         cons2.fa is built in the same frame as the aoplace reference (after the RCF flip), so it is not flipped again.
#  Output : ~/c2seq/MAPROWS_<host>.tsv   same 20 columns and same rules as CONSENSUS_MAP_G6
#           PARTIAL_PROX array_pos 0 = 1..first gene start-1
#           UNIT         array_pos 1..n = gene start -> next gene start
#           UNIT_LAST    array_pos n+1 = last gene start -> end
#           N = last array_pos (same as G6: NA19090 1 is UNIT 67 + UNIT_LAST 1 = 68)
#           if a spacing falls outside 2100~2400, write chain_ok 0 but keep the numbering continuous
#         ~/c2seq/SEQOVR_<host>.tsv    sample slot sequence path (sequence with N removed, the same sequence as the map coordinates)
#  Unit cutting matches the current pipeline (G6), at the gene HSP start. (G8 only on)
#  If several HSPs fall in the same 1 kb bin, keep only the longest.
import io
import os
import re
import sys
import socket
import subprocess

NL = chr(10)
TAB = chr(9)
TH = "8"
R = os.path.expanduser("~/s5s_tuc/ref")
CANON = os.path.join(R, "unit_canon.fa")
GENE = os.path.join(R, "5S_gene_query.fa")
CS = re.compile("([:*+-])([A-Za-z0-9]+)")
HOST = socket.gethostname()


def say(s):
    sys.stderr.write(s + NL)


def rd(p):
    return "".join(x.strip() for x in io.open(p, encoding="utf8").read().split(NL)[1:])


def run(cmd):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if r.returncode != 0:
        say("FAIL: %s %s" % (cmd, r.stderr[-400:]))
        raise SystemExit(1)
    return r.stdout


def starts(seq, W):
    fa = os.path.join(W, "mk.fa")
    io.open(fa, "w", newline=NL).write(">s" + NL + seq + NL)
    run("makeblastdb -in %s -dbtype nucl -out %s.db >/dev/null" % (fa, fa))
    out = run("blastn -task blastn -query %s -db %s.db -outfmt '6 sstart send length' "
              "-dust no -word_size 9 -evalue 1e-5 -max_target_seqs 100000000" % (GENE, fa))
    best = {}
    for L in out.splitlines():
        f = L.split(TAB)
        if len(f) < 3:
            continue
        s = min(int(f[0]), int(f[1])) - 1
        ln = int(f[2])
        b = s // 1000
        if b not in best or ln > best[b][1]:
            best[b] = (s, ln)
    return sorted(v[0] for v in best.values())


def varof(seqs, W):
    q = os.path.join(W, "mk.v.fa")
    with io.open(q, "w", newline=NL) as g:
        for i, s in enumerate(seqs):
            g.write(">U%d%s%s%s" % (i, NL, s, NL))
    out = run("minimap2 -x asm5 -c --cs -t %s --secondary=no %s %s" % (TH, CANON, q))
    V = {}
    for L in out.splitlines():
        f = L.split(TAB)
        if len(f) < 12:
            continue
        i = int(f[0][1:])
        if i in V:
            continue
        cs = ""
        for x in f[12:]:
            if x.startswith("cs:Z:"):
                cs = x[5:]
        tp = int(f[7])
        o = set()
        for m in CS.finditer(cs):
            op, v = m.group(1), m.group(2)
            if op == ":":
                tp += int(v)
            elif op == "*":
                o.add("%d%s" % (tp + 1, v[1].upper()))
                tp += 1
            elif op == "-":
                o.add("%dD%d" % (tp + 1, len(v)))
                tp += len(v)
            elif op == "+":
                o.add("%dI%d" % (tp + 1, len(v)))
        V[i] = o
    return V


MR = io.open(os.path.expanduser("~/c2seq/MAPROWS_%s.tsv" % HOST), "w", encoding="utf8", newline=NL)
SO = io.open(os.path.expanduser("~/c2seq/SEQOVR_%s.tsv" % HOST), "w", encoding="utf8", newline=NL)
for L in io.open(sys.argv[1], encoding="utf8"):
    p = L.split()
    if len(p) < 2:
        continue
    smp, sl = p[0], p[1]
    W = os.path.expanduser("~/c2seq/%s_%s" % (smp, sl))
    fa = os.path.join(W, "cons2.fa")
    if not os.path.exists(fa):
        say("  %s_%s cons2.fa 없음 - 건너뜀" % (smp, sl))
        continue
    raw = rd(fa)
    seq = raw.replace("N", "")
    st = starts(seq, W)
    if len(st) < 2:
        say("  %s_%s 유전자 HSP %d - 건너뜀" % (smp, sl, len(st)))
        continue
    segs = [("PARTIAL_PROX", 0, st[0])]
    for i in range(len(st) - 1):
        segs.append(("UNIT", st[i], st[i + 1]))
    segs.append(("UNIT_LAST", st[-1], len(seq)))
    V = varof([seq[a:b] for _t, a, b in segs], W)
    N = len(segs) - 1
    nbad = nv = 0
    for i, (t, a, b) in enumerate(segs):
        ln = b - a
        ok = 1
        if t == "UNIT" and not (2100 <= ln <= 2400):
            ok = 0
            nbad += 1
        gap = ln if t == "UNIT" else 0
        if t == "PARTIAL_PROX":
            row = [smp, sl, str(N), "0", t, "1", str(b), str(ln), "0", "1", "-", "-", "-", "-", "-"]
        else:
            vs = V.get(i)
            vl = ",".join(sorted(vs)) if vs else "-"
            nv += len(vs) if vs else 0
            ni = sum(1 for x in (vs or ()) if "I" in x[1:])
            nd = sum(1 for x in (vs or ()) if "D" in x[1:])
            row = [smp, sl, str(N), str(i), t, str(a + 1), str(b), str(ln), str(gap), str(ok),
                   str(len(vs)) if vs is not None else "-", str(ni), str(nd), "-", vl]
        MR.write(TAB.join(row + ["-"] * 5) + NL)
    arr = os.path.join(W, "cons2_arr.fa")
    io.open(arr, "w", newline=NL).write(">%s_%s_cons2%s%s%s" % (smp, sl, NL, seq, NL))
    SO.write(TAB.join([smp, sl, arr]) + NL)
    say("  %s_%s : 길이 %d (N %d 제거) · UNIT %d · UNIT_LAST 1 · N=%d · 간격밖 %d · 변이 %d"
        % (smp, sl, len(seq), raw.count("N"), len(st) - 1, N, nbad, nv))
MR.close()
SO.close()
print("done")
