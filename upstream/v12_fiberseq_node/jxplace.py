# -*- coding: utf-8 -*-
"""Junction placement: the position (array_pos) of a read whose slot is fixed is set by the flank anchor

  Before this step : the slot is fixed by slotvote.py from a vote on flank variants (accuracy 99.3% at margin >=0.90).
                The chain gate (K=12) is therefore not used here: that gate exists to
                fix the slot, and with the slot fixed a read with a single unit also receives a position.

  Steps
    (1) align the read to the canonical sequence of its slot with flanks added (map-ont, --secondary=no)
    (2) alignment coordinate -> canonical coordinate (- flank5) -> find the unit interval in the G6 map
    (3) per unit, cut the read sequence, align it to the canonical unit (unit_canon.fa) and call variants
    (4) six-column AOJV output : individual, read (16 characters), slot, array_pos, variant list, scored

  Coordinate convention : the coordinate in the flank-added file = the canonical coordinate + flank5
"""
import io
import os
import sys
import subprocess
import collections
import re

CIG = re.compile(r"(\d+)([MIDNSHP=X])")
CS = re.compile(r"([:*+\-])([A-Za-z0-9]+)")
H = os.path.expanduser
HSN = H(os.environ["HSN"])
CMPD = H(os.environ["CMPD"])
MAPF = H(os.environ["MAPF"])
READS = H(os.environ["READS"])
CANON = H("~/s5s_tuc/ref/unit_canon.fa")
SUMF = H(os.environ["SUMF"])
C_TR = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def rc(x):
    return x.translate(C_TR)[::-1]

TH = os.environ.get("TH", "6")
MARGIN = float(os.environ.get("MARGIN", "0.90"))
# the read-name length: 16 for ONT and 0, the whole name, for PacBio, whose first 16 characters are the movie and would collide in every case
RIDLEN = int(os.environ.get("RIDLEN", "16"))
VOTE, OUT = sys.argv[1], sys.argv[2]
NL = "\n"

SEQF = {}
for x in io.open(HSN, encoding="utf8"):
    f = x.rstrip(NL).split("\t")
    if len(f) >= 3:
        SEQF[(f[0], f[1])] = f[2]

# the revcomp table, located by column name: the 6th in G6 and the 5th in the G8 family
RCF = {}
for i, x in enumerate(io.open(SUMF, encoding="utf8")):
    fz = x.rstrip(NL).split("	")
    if i == 0:
        if "revcomp" not in fz:
            raise SystemExit("SUMF has no revcomp column: %s" % SUMF)
        RCI = fz.index("revcomp")
        continue
    if len(fz) > RCI:
        RCF[(fz[0], fz[1])] = (fz[RCI] == "Y")

# the map: (individual, slot) -> [(array_pos, start0, end, N)], in canonical coordinates with no flank
UB = collections.defaultdict(list)
f = io.open(MAPF, encoding="utf8")
h = f.readline().rstrip(NL).split("\t")
ix = dict((n, i) for i, n in enumerate(h))
for L in f:
    r = L.rstrip(NL).split("\t")
    if len(r) < 15 or r[ix["seg_type"]].startswith("PARTIAL"):
        continue
    UB[(r[0], r[1])].append((int(r[ix["array_pos"]]), int(r[ix["unit_start"]]) - 1,
                             int(r[ix["unit_end"]]), int(r[ix["N"]])))
for k in UB:
    UB[k].sort()
CANLEN = dict((k, max(e for (_u, _s0, e, _N) in v)) for k, v in UB.items())


def load(p):
    hdr = ""
    buf = []
    for L in io.open(p, encoding="utf8"):
        if L.startswith(">"):
            hdr = L.strip()
        else:
            buf.append(L.strip())
    d = dict(t.split("=", 1) for t in hdr.split("|")[1:] if "=" in t)
    return "".join(buf).upper(), int(d.get("flank5", 0))


# reading the slot assignment
ASSIGN = collections.defaultdict(dict)      # individual -> read -> slot
for i, L in enumerate(io.open(VOTE, encoding="utf8")):
    if i == 0:
        continue
    r = L.rstrip(NL).split("\t")
    if len(r) < 8:
        continue
    if float(r[7]) < MARGIN:
        continue
    ASSIGN[r[0]][r[1]] = "1" if int(r[5]) > int(r[6]) else "2"

o = io.open(OUT, "w", encoding="utf8", newline=NL)
T = collections.Counter()

for smp in sorted(ASSIGN):
    rd2sl = ASSIGN[smp]
    fq = os.path.join(READS, "%s.hit.fq" % smp)
    if not os.path.exists(fq):
        continue
    SEQ = {}
    it = io.open(fq, encoding="utf8")
    while True:
        hh = it.readline()
        if not hh:
            break
        sq = it.readline().strip()
        it.readline()
        it.readline()
        rid = hh[1:].split()[0].split("\t")[0]
        if rid in rd2sl:
            SEQ[rid] = sq
    if not SEQ:
        continue

    for sl in ("1", "2"):
        ids = [r for r in SEQ if rd2sl[r] == sl]
        if not ids:
            continue
        key = (smp, sl)
        if key not in UB or (smp, sl) not in SEQF:
            continue
        p = os.path.join(CMPD, SEQF[(smp, sl)])
        if not os.path.exists(p):
            continue
        s, f5 = load(p)
        # Polarity: the unit_start and unit_end of the map are in the frame after revcomp has been applied, as at fsplace lines 206-210.
        #   The flanked canonical sequence is stored in its original orientation, so a slot with revcomp = Y must be flipped here.
        #   Without the flip the unit indices come out reversed throughout, as N + 1 - p.
        if RCF.get((smp, sl)):
            f3 = len(s) - f5 - CANLEN.get((smp, sl), 0)
            s = rc(s)
            f5 = f3                      # flipping it puts the 3' flank where the 5' flank was
        ref = "/tmp/jx_%s_%s.fa" % (smp, sl)
        io.open(ref, "w", newline=NL).write(">%s_h%s%s%s%s" % (smp, sl, NL, s, NL))
        q = "/tmp/jx_%s_%s.q.fa" % (smp, sl)
        with io.open(q, "w", newline=NL) as g:
            for r in ids:
                g.write(">%s%s%s%s" % (r, NL, SEQ[r], NL))

        # 1) the read against its own slot's reference
        pr = subprocess.Popen(["minimap2", "-a", "-x", "map-ont", "-t", TH,
                               "--secondary=no", ref, q],
                              stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
        cut = []        # (read, array_pos, fragment sequence)
        for L in pr.stdout:
            if L.startswith("@"):
                continue
            fl = L.rstrip(NL).split("\t")
            if len(fl) < 10 or fl[2] == "*":
                continue
            if int(fl[1]) & 0x900:
                continue
            rid = fl[0]
            seq = fl[9]
            rpos = int(fl[3]) - 1        # the reference coordinate, that is the canonical sequence plus flank5
            qp = 0
            m2 = {}
            for mm in CIG.finditer(fl[5]):
                n_, op = int(mm.group(1)), mm.group(2)
                if op == "S":
                    qp += n_
                elif op == "H":
                    pass
                elif op in "M=X":
                    for t in range(n_):
                        m2[rpos + t] = qp + t
                    rpos += n_
                    qp += n_
                elif op == "I":
                    qp += n_
                elif op in "DN":
                    rpos += n_
            if not m2:
                continue
            # 2) the map's unit intervals are carried into reference coordinates, by adding f5, and only the units the read covers entirely are kept
            for (u, s0, e0, N) in UB[key]:
                A, B = s0 + f5, e0 + f5 - 1
                if A in m2 and B in m2:
                    cut.append((rid, u, seq[m2[A]:m2[B] + 1]))
        pr.stdout.close()
        pr.wait()
        for x in (ref, q):
            try:
                os.remove(x)
            except OSError:
                pass
        if not cut:
            continue

        # 3) the variants are called per unit
        uq = "/tmp/jx_%s_%s.u.fa" % (smp, sl)
        with io.open(uq, "w", newline=NL) as g:
            for i, (_r, _u, sq) in enumerate(cut):
                g.write(">U%d%s%s%s" % (i, NL, sq, NL))
        rr = subprocess.run(["minimap2", "-x", "asm5", "-c", "--cs", "-t", TH,
                             "--secondary=no", CANON, uq],
                            capture_output=True, text=True)
        V = {}
        for L in rr.stdout.splitlines():
            fl = L.split("\t")
            if len(fl) < 12:
                continue
            i = int(fl[0][1:])
            if i in V:
                continue
            cs = ""
            for x in fl[12:]:
                if x.startswith("cs:Z:"):
                    cs = x[5:]
            tp = int(fl[7])
            out = set()
            for mm in CS.finditer(cs):
                op, v = mm.group(1), mm.group(2)
                if op == ":":
                    tp += int(v)
                elif op == "*":
                    out.add("%d%s" % (tp + 1, v[1].upper()))
                    tp += 1
                elif op == "-":
                    out.add("%dD%d" % (tp + 1, len(v)))
                    tp += len(v)
                elif op == "+":
                    out.add("%dI%d" % (tp + 1, len(v)))
            V[i] = out
        try:
            os.remove(uq)
        except OSError:
            pass

        for i, (rid, u, _sq) in enumerate(cut):
            vv = V.get(i)
            if vv is None:
                T["변이못부름"] += 1
                continue
            o.write("%s\t%s\t%s\t%d\t%s\t1%s"
                    % (smp, (rid[:RIDLEN] if RIDLEN else rid), sl, u, (",".join(sorted(vv)) or "-"), NL))
            T["유닛"] += 1
        T["read"] += len(set(r for r, _u, _s in cut))
    sys.stderr.write("  %-10s read %d · 유닛 %d%s" % (smp, T["read"], T["유닛"], NL))
o.close()
sys.stderr.write("  합계 " + " · ".join("%s %d" % kv for kv in sorted(T.items())) + NL)
