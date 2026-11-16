# -*- coding: utf-8 -*-
# 5S in other species: the G8 unit map for an array with several reference units, choosing the best-matching reference for each unit.
#
# This follows species_g8map.py (single reference), changed only where having several references makes a difference.
#
# Only one thing becomes more complicated: the length window that decides whether a gap is accepted as a unit.
#   The unit boundary itself comes from the gene ends and does not depend on the reference (gene end + 1 .. next gene end).
#   Whether the gap is accepted as a unit is decided by the length window LO..HI, and that window comes from the reference's period.
#   Two references mean two windows, so this script:
#     1) builds a separate window [LO_k, HI_k] per group from its own period (the median unit length of its members)
#     2) accepts a gap as a unit if it falls in any one window (the union of the windows), recording which window accepted it
#     3) checks for multiples: if one window's HI is at least twice another window's LO, a missed gene can
#        let two merged units pass as a normal unit of the other type, so a warning is recorded
#     4) records where the windows overlap: in that range, length alone cannot tell which type a gap is
#        (which is acceptable, because the type is decided by sequence alignment and not by length)
#     5) counts the units recovered relative to using the single array reference.
#   Because the windows are unioned, no unit can be lost: a single window is a subset of the union.
#
# Choosing a reference per unit: align to the context version of every reference and take the top score. The confidence is normalised by Dmax, the distance between the references.
#   margin_frac >= 0.50 is OK, 0.20-0.50 is WEAK, below 0.20 is AMBIG and falls back to the common frame of the array reference.
# No filling: variants come from the alignment (cs / BTOP) alone. If the chosen reference's unit interval is not fully covered, the result is NO_FULL.
#
# Arguments: array_id, fasta, gene.fa, reference folder, group files (comma separated), PAD, CTX, GL, output folder, array start, array end, chrom, strand.
# MIT License · Copyright (c) 2026 Hungmin Kim
import io
import os
import sys
import csv
import subprocess
import collections

NL = chr(10)
TAB = chr(9)
(AID, FA, GENEF, RDIR, GFILES, PADS, CTXS, GLS, OUT, S0, E1, CHROM, STRAND) = sys.argv[1:14]
PAD, CTX, GL = int(PADS), int(CTXS), int(GLS)
s0, e1 = int(S0), int(E1)
os.makedirs(os.path.join(OUT, "g8map"), exist_ok=True)
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
rc = lambda s: "".join(COMP.get(x, "N") for x in reversed(s))


def rdfa(p):
    return "".join(x.strip() for x in io.open(p, encoding="utf8").read().split(NL)[1:]).upper()


def faidx(reg):
    r = subprocess.run(["samtools", "faidx", FA, reg], capture_output=True, text=True)
    return "".join(r.stdout.splitlines()[1:]).upper()


def cs_vars(cs, rs):
    out, p, i = [], rs, 0
    while i < len(cs):
        ch = cs[i]
        if ch == ":":
            j = i + 1
            while j < len(cs) and cs[j].isdigit():
                j += 1
            p += int(cs[i + 1:j])
            i = j
        elif ch == "*":
            out.append((p + 1, "S", cs[i + 2].upper()))
            p += 1
            i += 3
        elif ch == "-":
            j = i + 1
            while j < len(cs) and cs[j].isalpha():
                j += 1
            out.append((p + 1, "D", j - i - 1))
            p += j - i - 1
            i = j
        elif ch == "+":
            j = i + 1
            while j < len(cs) and cs[j].isalpha():
                j += 1
            out.append((p + 1, "I", j - i - 1))
            i = j
        else:
            i += 1
    return out


def fmt(tok):
    out, i = [], 0
    while i < len(tok):
        p, k, v = tok[i]
        if k == "S":
            out.append("%d%s" % (p, v))
            i += 1
        elif k == "D":
            n = v
            while i + 1 < len(tok) and tok[i + 1][1] == "D" and tok[i + 1][0] == p + n:
                n += tok[i + 1][2]
                i += 1
            out.append("%dD%d" % (p, n))
            i += 1
        else:
            n = v
            while i + 1 < len(tok) and tok[i + 1][1] == "I" and tok[i + 1][0] == p:
                n += tok[i + 1][2]
                i += 1
            out.append("%dI%d" % (p, n))
            i += 1
    return out


def window(period):
    # The human G8v4 rule carried across to other species: whichever is the more permissive of -6% / +7% of the period and +-30 bp.
    return (int(min(period * 0.94, period - 30)), int(max(period * 1.07, period + 30)))


LEN = dict((L.split(TAB)[0], int(L.split(TAB)[1])) for L in io.open(FA + ".fai", encoding="utf8"))
refs = collections.OrderedDict()
refs["arrayref"] = rdfa(os.path.join(RDIR, AID + ".unit.fa"))
PERIOD = {"arrayref": len(refs["arrayref"])}
for f in [x for x in GFILES.split(",") if x]:
    k = os.path.basename(f).split(".")[-3]
    refs[k] = rdfa(f)
    h = io.open(f, encoding="utf8").read().split(NL)[0]
    PERIOD[k] = len(refs[k])
    for tok in h.split():
        if tok.startswith("period="):
            PERIOD[k] = int(tok[7:])
WIN = dict((k, window(PERIOD[k])) for k in refs)
GK = [k for k in refs if k != "arrayref"]
ULO = min(w[0] for w in WIN.values())
UHI = max(w[1] for w in WIN.values())
# check whether one window's HI is at least twice another window's LO
warn = []
for a in WIN:
    for b in WIN:
        if a != b and WIN[a][1] >= 2 * WIN[b][0]:
            warn.append("%s 창 상한 %d 이 %s 창 하한 %d 의 2 배 이상" % (a, WIN[a][1], b, WIN[b][0]))
ov = (max(w[0] for w in WIN.values()), min(w[1] for w in WIN.values()))
ovtxt = "%d-%d" % ov if ov[0] <= ov[1] else "없음"

ctxp = os.path.join(OUT, "g8map", AID + ".ctx.fa")
with io.open(ctxp, "w", encoding="utf8", newline=NL) as fh:
    for k, s in refs.items():
        fh.write(">%s%s%s%s" % (k, NL, s[-CTX:] + s + s[:CTX], NL))


def aln_score(qseq, tseq):
    qp = os.path.join(OUT, "g8map", AID + ".q.fa")
    tp = os.path.join(OUT, "g8map", AID + ".t.fa")
    io.open(qp, "w", encoding="utf8", newline=NL).write(">q%s%s%s" % (NL, qseq, NL))
    io.open(tp, "w", encoding="utf8", newline=NL).write(
        ">t%s%s%s" % (NL, tseq[-CTX:] + tseq + tseq[:CTX], NL))
    r = subprocess.run("minimap2 -cx asm20 -t 2 %s %s 2>/dev/null" % (tp, qp),
                       shell=True, capture_output=True, text=True)
    b = 0
    for L in r.stdout.splitlines():
        f = L.split(TAB)
        for x in (f[12:] if len(f) > 12 else []):
            if x.startswith("AS:i:"):
                b = max(b, int(x[5:]))
    os.remove(qp)
    os.remove(tp)
    return b


DM = 0
if len(GK) >= 2:
    DM = min((aln_score(refs[k], refs[k]) - max(aln_score(refs[k], refs[j]) for j in GK if j != k))
             for k in GK)

# ---- build the boundaries from the gene positions, independently of any reference ----
mar = 3 * max(PERIOD.values()) + 3000
L0, R1 = max(1, s0 - mar), min(LEN[CHROM], e1 + mar)
seq = faidx("%s:%d-%d" % (CHROM, L0, R1))
rp = os.path.join(OUT, "g8map", AID + ".region.fa")
io.open(rp, "w", encoding="utf8", newline=NL).write(">region%s%s%s" % (NL, seq, NL))
rb = subprocess.run("blastn -task blastn -word_size 11 -dust no -evalue 1e-5 -query %s -subject %s "
                    "-outfmt '6 qstart qend sstart send bitscore' 2>/dev/null" % (GENEF, rp),
                    shell=True, capture_output=True, text=True)
G = []
for L in rb.stdout.splitlines():
    q1, q2, x1, x2, bs = L.split(TAB)
    q1, q2, x1, x2 = int(q1), int(q2), int(x1), int(x2)
    if x1 <= x2:
        G.append((x1 - (q1 - 1), x2 + (GL - q2), "+"))
    else:
        G.append((x2 - (GL - q2), x1 + (q1 - 1), "-"))
G.sort()
keep, last = [], -1
for (a, b, st) in G:
    if a > last:
        keep.append((a, b, st))
        last = b
major = collections.Counter(x[2] for x in keep).most_common(1)[0][0] if keep else STRAND
if major == "-":
    n = len(seq)
    wseq = rc(seq)
    gg = sorted([(n - b + 1, n - a + 1) for (a, b, st) in keep if st == major])
else:
    wseq = seq
    gg = sorted([(a, b) for (a, b, st) in keep if st == major])

qs, meta, brk = [], [], 0
pos = 0
n_rescued = 0
for i in range(len(gg) - 1):
    a, b = gg[i][1] + 1, gg[i + 1][1]
    L = b - a + 1
    acc = [k for k in refs if WIN[k][0] <= L <= WIN[k][1]]
    if not acc:
        brk += 1
        continue
    if "arrayref" not in acc:
        n_rescued += 1
    pos += 1
    la, lb = max(1, a - PAD), min(len(wseq), b + PAD)
    meta.append((pos, a, b, la, lb, L, ",".join(acc)))
    qs.append(">u%d%s%s%s" % (pos, NL, wseq[la - 1:lb], NL))
qp = os.path.join(OUT, "g8map", AID + ".units.fa")
io.open(qp, "w", encoding="utf8", newline=NL).write("".join(qs))

# ---- align the unit to every reference ----
best = collections.defaultdict(dict)
rm = subprocess.run("minimap2 -cx asm20 --cs -N 10 -p 0.5 -t 4 %s %s 2>/dev/null" % (ctxp, qp),
                    shell=True, capture_output=True, text=True)
for L in rm.stdout.splitlines():
    f = L.split(TAB)
    if len(f) < 12:
        continue
    q, tgt, rs, re_ = f[0], f[5], int(f[7]), int(f[8])
    bs, cs = 0, ""
    for x in f[12:]:
        if x.startswith("AS:i:"):
            bs = int(x[5:])
        elif x.startswith("cs:Z:"):
            cs = x[5:]
    if tgt not in best[q] or bs > best[q][tgt][0]:
        best[q][tgt] = (bs, cs_vars(cs, rs), rs + 1, re_)

HDR = ["array_id", "chrom", "strand", "array_pos", "N", "chrom_start", "chrom_end", "unit_len", "gap_to_next",
       "len_window_accepted_by", "ref_used", "ref_best_raw", "score_best", "score_second", "margin_frac",
       "assign", "cov_start", "cov_end", "n_var", "n_var_gene", "varlist", "varlist_arrayref",
       "verify"] + ["score_%s" % k for k in refs]
M = io.open(os.path.join(OUT, "g8map", AID + ".g8mr.tsv"), "w", encoding="utf8", newline=NL)
M.write(TAB.join(HDR) + NL)
N = len(meta)
cnt, acnt, rcnt = collections.Counter(), collections.Counter(), collections.Counter()
for k, (p, a, b, la, lb, ulen, accby) in enumerate(meta):
    h = best.get("u%d" % p, {})
    gap = (meta[k + 1][1] - b - 1) if k + 1 < len(meta) else -1
    cs_, ce_ = (R1 - b + 1, R1 - a + 1) if major == "-" else (L0 + a - 1, L0 + b - 1)
    d = dict(array_id=AID, chrom=CHROM, strand=major, array_pos=p, N=N, chrom_start=cs_, chrom_end=ce_,
             unit_len=ulen, gap_to_next=gap, len_window_accepted_by=accby)
    for x in refs:
        d["score_%s" % x] = h[x][0] if x in h else 0

    def call(key):
        if key not in h:
            return None
        bs, toks, c0, c1 = h[key]
        U = len(refs[key])
        t = [(pp - CTX, kk, vv) for (pp, kk, vv) in toks if CTX < pp <= CTX + U]
        return (bs, fmt(t), max(1, c0 - CTX), min(U, c1 - CTX), (c0 - CTX <= 1 and c1 - CTX >= U),
                sum(1 for (pp, kk, vv) in t if pp > U - GL))

    ar = call("arrayref")
    d["varlist_arrayref"] = ",".join(ar[1]) if ar and ar[1] else ("-" if ar else "NOALN")
    cand = sorted(((h[x][0], x) for x in GK if x in h), reverse=True)
    if not cand:
        d.update(dict(ref_used="NOALN", ref_best_raw="NOALN", assign="NOALN", verify="NOALN",
                      score_best=0, score_second=0, margin_frac="-", cov_start="NA", cov_end="NA",
                      n_var="NA", n_var_gene="NA", varlist="NOALN"))
        cnt["NOALN"] += 1
    else:
        s1, k1 = cand[0]
        s2 = cand[1][0] if len(cand) > 1 else 0
        frac = (s1 - s2) / float(DM) if DM > 0 else 1.0
        lv = "OK" if frac >= 0.50 else ("WEAK" if frac >= 0.20 else "AMBIG")
        used = "arrayref" if lv == "AMBIG" else k1
        cc = call(used) or ar
        d.update(dict(ref_used=used, ref_best_raw=k1, score_best=s1, score_second=s2,
                      margin_frac="%.3f" % frac, assign=lv,
                      cov_start=cc[2], cov_end=cc[3], n_var=len(cc[1]), n_var_gene=cc[5],
                      varlist=",".join(cc[1]) if cc[1] else "-",
                      verify="OK" if cc[4] else "NO_FULL"))
        cnt[d["verify"]] += 1
        acnt[lv] += 1
        rcnt[used] += 1
    M.write(TAB.join(str(d.get(x, "-")) for x in HDR) + NL)
M.close()
for x in (rp, qp, ctxp):
    try:
        os.remove(x)
    except OSError:
        pass
S = io.open(os.path.join(OUT, "g8map", AID + ".g8mr.info.tsv"), "w", encoding="utf8", newline=NL)
S.write(TAB.join(["array_id", "n_units", "n_break", "n_rescued_by_group_window", "verify_OK", "verify_NO_FULL",
                  "verify_NOALN", "assign_OK", "assign_WEAK", "assign_AMBIG", "ref_pair_Dmax",
                  "windows", "window_overlap", "union_window", "multiple_period_warning"]) + NL)
S.write(TAB.join([AID, str(N), str(brk), str(n_rescued), str(cnt["OK"]), str(cnt["NO_FULL"]), str(cnt["NOALN"]),
                  str(acnt["OK"]), str(acnt["WEAK"]), str(acnt["AMBIG"]), str(DM),
                  " ".join("%s:%d-%d(주기%d)" % (k, WIN[k][0], WIN[k][1], PERIOD[k]) for k in refs),
                  ovtxt, "%d-%d" % (ULO, UHI), " · ".join(warn) if warn else "없음"]) + NL)
S.close()
print("%s : 유닛 %d (끊김 %d · 그룹 창이 살린 유닛 %d) · OK %d NO_FULL %d NOALN %d · 배정 %s · 참조 %s%s"
      % (AID, N, brk, n_rescued, cnt["OK"], cnt["NO_FULL"], cnt["NOALN"], dict(acnt), dict(rcnt),
         (" · 배수함정 " + warn[0]) if warn else ""))
