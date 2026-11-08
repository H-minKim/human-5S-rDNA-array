# -*- coding: utf-8 -*-
# Recall read units: cut with a margin (±PAD) and align to the context-padded canonical unit (2026-09-12 user design)
#   reference = [last CTX bp of the canonical unit][canonical unit 2,240bp][first CTX bp of the canonical unit]
#     In a tandem array, what precedes and follows a unit is the end and start of the same canonical unit. Adding context puts the margin on both sides,
#     and the central unit is free of alignment edge problems. Joining two canonical units would split which copy it attaches to, so that is not used.
#   why: aligning exactly cut units trims the tails short (524 of 2,172 in the pilot), and filling them (character matching)
#        lets one indel shift everything after it, creating false substitutions. G8 units have the gene (2122-2240) at the end, so the damage concentrates in the gene part.
#   variants: only the alignment result (cs) is used. No filling. Of target coordinates t (0-based), only those in [CTX, CTX+2240) are written as canonical unit coordinates (t-CTX+1).
#   status: OK (alignment covers the whole unit segment), NO_FULL, SEQ_MISMATCH, NOALN
#   args: JX|AO, input unit table, output table   (env var PAD default 50, CTX default 100)
import io
import os
import re
import sys
import subprocess
import collections

NL = chr(10)
TAB = chr(9)
H = os.path.expanduser
KIND, INP, OUTP = sys.argv[1], sys.argv[2], sys.argv[3]
PAD = int(os.environ.get("PAD", "50"))
CTX = int(os.environ.get("CTX", "100"))
CANON = H("~/s5s_tuc/ref/unit_canon.fa")
CS = re.compile(r"([-:*+=])([A-Za-z0-9]+)")
CANSEQ = "".join(x.strip() for x in io.open(CANON, encoding="utf8").read().split(NL)[1:]).upper()
UL = len(CANSEQ)
REFSEQ = CANSEQ[-CTX:] + CANSEQ + CANSEQ[:CTX]
REF = "/tmp/canon_ctx%d.fa" % CTX
if not os.path.exists(REF):
    io.open(REF, "w", newline=NL).write(">canon_ctx" + NL + REFSEQ + NL)


def call(seqs, tag):
    q = "/tmp/rp2_%s.fa" % tag
    with io.open(q, "w", newline=NL) as g:
        for i, s in enumerate(seqs):
            g.write(">U%d%s%s%s" % (i, NL, s, NL))
    r = subprocess.run(["minimap2", "-x", "asm20", "-c", "--cs=long", "-t", "4", "--secondary=no", REF, q],
                       capture_output=True, text=True)
    OUT = {}
    for L in r.stdout.splitlines():
        f2 = L.split(TAB)
        if len(f2) < 12:
            continue
        i = int(f2[0][1:])
        if i in OUT:
            continue
        cs = ""
        for x in f2[12:]:
            if x.startswith("cs:Z:"):
                cs = x[5:]
        OUT[i] = (int(f2[2]), int(f2[3]), int(f2[7]), int(f2[8]), cs)
    os.remove(q)
    return OUT


def tokens(qseq, g):
    """Extract only the variants in the unit segment [CTX, CTX+UL) of the alignment (cs), in canonical unit coordinates"""
    qs, qe, ts, te, cs = g
    tp = ts
    toks, reb = [], []
    for m in CS.finditer(cs):
        op, v = m.group(1), m.group(2)
        if op == "=":
            reb.append(v.upper())
            tp += len(v)
        elif op == ":":
            n = int(v)
            reb.append(REFSEQ[tp:tp + n])
            tp += n
        elif op == "*":
            if CTX <= tp < CTX + UL:
                toks.append("%d%s" % (tp - CTX + 1, v[1].upper()))
            reb.append(v[1].upper())
            tp += 1
        elif op == "-":
            if CTX <= tp < CTX + UL:
                toks.append("%dD%d" % (tp - CTX + 1, len(v)))
            tp += len(v)
        elif op == "+":
            if CTX <= tp < CTX + UL:
                toks.append("%dI%d" % (tp - CTX + 1, len(v)))
            reb.append(v.upper())
    full = (ts <= CTX) and (te >= CTX + UL)
    return sorted(set(toks)), "".join(reb), full


if KIND == "JX":
    C_SMP, C_RID, C_JOB, C_F0, C_F1, C_VD, C_VL = 0, 1, 3, 7, 8, 11, 12
else:
    C_SMP, C_RID, C_JOB, C_F0, C_F1, C_VD, C_VL = 0, 1, 4, 9, 10, 13, 14
rows = [L.rstrip(NL).split(TAB) for L in io.open(INP, encoding="utf8") if L.strip()]
by = collections.defaultdict(list)
for i, r in enumerate(rows):
    by[(r[C_SMP], r[C_JOB])].append(i)


def samfile(smp, job):
    if KIND == "JX":
        return H("~/rjrepro/%s.placed.sam" % smp) if job == "JX" else H("~/rjrepro_c2/%s.placed.sam" % smp)
    return H("~/fix1/aog8f1/sam/%s.%s.sam" % (smp, job))


T = collections.Counter()
for (smp, job), idx in sorted(by.items()):
    sp = samfile(smp, job)
    if not os.path.exists(sp):
        T["SAM 없음"] += len(idx)
        continue
    want = set(rows[i][C_RID] for i in idx)
    SEQ = {}
    for L in io.open(sp, encoding="utf8"):
        if L.startswith("@"):
            continue
        f = L.split(TAB, 11)
        if f[0] in want and f[0] not in SEQ:
            SEQ[f[0]] = (f[9], int(f[1]) & 16)
    seqs, meta = [], []
    for i in idx:
        r = rows[i]
        if r[C_RID] not in SEQ:
            T["read 없음"] += 1
            continue
        s, rev = SEQ[r[C_RID]]
        L0 = len(s)
        f0, f1 = int(r[C_F0]), int(r[C_F1])
        a, b = (L0 - f1, L0 - f0) if rev else (f0, f1)
        a2, b2 = max(0, a - PAD), min(L0, b + PAD)
        seqs.append(s[a2:b2])
        meta.append((i, a - a2, b2 - b))
    if not seqs:
        continue
    V = call(seqs, "%s_%s_%s" % (KIND, smp, job))
    for k, (i, padl, padr) in enumerate(meta):
        r = rows[i]
        g = V.get(k)
        note = "pad%d;%d" % (padl, padr)
        if g is None:
            r[C_VD], r[C_VL] = "NOALN", "-"
            T["NOALN"] += 1
        else:
            toks, reb, full = tokens(seqs[k], g)
            if reb != seqs[k][g[0]:g[1]]:
                r[C_VD], r[C_VL] = "SEQ_MISMATCH", "-"
                T["SEQ_MISMATCH"] += 1
            elif not full:
                r[C_VD], r[C_VL] = "NO_FULL", "-"
                T["NO_FULL"] += 1
                T["  NO_FULL 중 read 끝이라 여유 부족"] += (padl < PAD or padr < PAD)
            else:
                r[C_VD], r[C_VL] = "OK", (",".join(toks) or "-")
                T["OK"] += 1
        if len(r) > (17 if KIND == "AO" else 13):
            r[-1] = note
        else:
            r.append(note)
with io.open(OUTP, "w", encoding="utf8", newline=NL) as o:
    for r in rows:
        o.write(TAB.join(r) + NL)
print("%s %s · %s" % (KIND, os.path.basename(INP), " · ".join("%s %d" % kv for kv in sorted(T.items()))))
