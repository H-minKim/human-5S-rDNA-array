# -*- coding: utf-8 -*-
# G8v3 (2026-09-11): only one gene boundary changed from G8v2
#   blastn gene HSPs stop one or two bases early when the gene end base is a variant (50 of 57,878 HSPs measured in the 616 canonical sequences)
#   G8v2 used only the canonical coordinates of the HSP, so that boundary shifted -> extend to the full gene length with the query coordinates (qstart/qend) to set the boundary
#     + strand : start = sstart - (qstart - 1), end = send + (119 - qend)
#     - strand : start = send - (119 - qend)   , end = sstart + (qstart - 1)
#   work folder ~/g8v3. Everything else (orientation, cutting rule, variant calling, checks) is as in G8v2
# G8v2 (2026-09-11): recut from the latest canonical sequences (those used by merged v7) following METHOD_CONSENSUS_MAP_G8_REBUILD.md
#   Changes are only the four places of the METHOD §5-5 table plus two for recording. Cutting, variant and check methods are as in G8.
#   (1) SEQF  : new hap_seq_node from env var HSN (stops if missing, so it never silently runs on the old file)
#   (2) .fa   : new folder CMPD (default ~/cmp4d_v7). The old ~/cmp4d is left as is
#   (3) RCF   : removed, strand majority only (METHOD §4). The revcomp column of SUM = whether it was flipped by strand majority
#   (4) W     : ~/g8v2. The old ~/g8 is not overwritten
#   (5) slot_space column added (= assembly, METHOD §3-3). Slot numbers are as in the old hap_seq_node
#   (6) no canonical, no file, and - strand still dominant after flipping are counted, not silently skipped
# consensus map G8: recut assembly arrays from scratch (2026-09-10)
#
#  defect of G6 and the stopgap of G7
#    G6 : cut units at the gene start -> [gene 119bp][AluY..NTS2]
#         the canonical unit is [AluY..NTS2][gene], so a 119bp rotation offset.
#         minimap2 cannot wrap around, so canonical 2122-2240 was never covered.
#         result: map variants stopped at 2101 and there were 0 5S gene variants.
#    G7 : patched by rotating the cut units by 119bp afterwards. It worked, but
#         as a patch on top of a wrong cut it cannot be canonical.
#
#  G8 = cut in canonical phase from the start
#    canonical unit = [AluY 1-318][NTS1 319-1178][GArich 1179-1551][NTS2 1552-2121][gene 2122-2240]
#    => the canonical unit starts right after the gene ends.
#    => therefore unit boundaries are set at the end of the gene HSP.
#         unit i = (gene_end[i] + 1) .. gene_end[i+1]
#       No rotation is needed. A cut unit is already in canonical phase.
#
#  flanks are already handled (shared)
#    cmp4d/*.fa are flank-to-flank sequences cut by asm_units3.py with 50kb flank probes.
#    The boundaries of the .fa are the junctions. Flanks are not handled again here.
#
#  proximal no. 1 (PARTIAL_PROX) is excluded
#    from the junction to the end of the first gene is not a complete unit (measured median 1,516 bp).
#    The leftover piece at the end of the array likewise. Both are marked PARTIAL and
#    no variants are called on them.
#
#  checks are done with sequence
#    ① restore unit sequences from the variant list and compare character by character with the original
#    ② ends not covered by the aligner are filled by direct sequence comparison
#       (minimap2 does not extend when the end mismatches -> it trims that variant.
#        in measurement, 237 of 400 cases were real variants)
#    ③ unit boundaries are checked with gene HSPs. No bp/2240 conversion
import io
import os
import re
import sys
import subprocess
import collections

NL = "\n"
W = os.path.expanduser(os.environ.get("WDIR", "~/g8v4"))
os.makedirs(W, exist_ok=True)
R = os.path.expanduser("~/s5s_tuc/ref")
CANON = os.path.join(R, "unit_canon.fa")
GENE = os.path.join(R, "5S_gene_query.fa")
TH = os.environ.get("TH", "4")
TAG = os.environ.get("TAG", "x")
PRESET = os.environ.get("PRESET", "asm20")
LO = int(os.environ.get("LO", "2100"))          # lower bound of the allowed unit spacing
HI = int(os.environ.get("HI", "2400"))          # upper bound
FRAG = int(os.environ.get("FRAG", "5"))         # if the extension is larger than this, treat it as a fragment HSP
CMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
CS = re.compile(r"([:*+\-=])([A-Za-z0-9]+)")
NUM = re.compile(r"\d+")


def rc(s):
    return s.translate(CMP)[::-1]


def rd(p):
    return "".join(x.strip() for x in io.open(p, encoding="utf8").read().split(NL)[1:])


CANSEQ = rd(CANON).upper()
PAD = int(os.environ.get("PAD", "100"))           # margin before and after when cutting units
CTX = int(os.environ.get("CTX", "200"))           # context on both sides of the canonical unit (tandem, so the end and start of the same unit)
REFSEQ = CANSEQ[-CTX:] + CANSEQ + CANSEQ[:CTX]
CANREF = os.path.join(W, "canon_ctx.fa")
io.open(CANREF, "w", newline=NL).write(">canon_ctx" + NL + REFSEQ + NL)
if len(CANSEQ) != 2240:
    sys.stderr.write("정본유닛 길이가 2240 이 아니다: %d%s" % (len(CANSEQ), NL))
    raise SystemExit(1)
GENESEQ = rd(GENE).upper()
if len(GENESEQ) != 119:
    sys.stderr.write("유전자 질의 길이가 119 가 아니다: %d%s" % (len(GENESEQ), NL))
    raise SystemExit(1)
# check with sequence that the gene is at 2122..2240 within the canonical unit (do not assume it)
if CANSEQ[2121:2240] != GENESEQ:
    sys.stderr.write("정본 2122-2240 이 유전자와 다르다 — 중단%s" % NL)
    raise SystemExit(1)

SEQF = {}
for x in io.open(os.path.expanduser(os.environ["HSN"]), encoding="utf8"):
    f = x.rstrip(NL).split("\t")
    if len(f) >= 3 and f[0] != "sample":
        SEQF[(f[0], f[1])] = f[2]
RCF = {}
for i, x in enumerate([]):   # RCF not used (METHOD §4)
    if i == 0:
        continue
    f = x.rstrip(NL).split("\t")
    if len(f) >= 6:
        RCF[(f[0], f[1])] = (f[5] == "Y")
FLIP = set()


def gene_hits(seq, tag):
    """(start, end) of gene HSPs and the dominant strand. Coordinates are 1-based inclusive."""
    fa = os.path.join(W, "g_%s.fa" % tag)
    db = os.path.join(W, "gd_%s" % tag)
    io.open(fa, "w", newline=NL).write(">x" + NL + seq + NL)
    subprocess.run(["makeblastdb", "-in", fa, "-dbtype", "nucl", "-out", db],
                   capture_output=True)
    r = subprocess.run(["blastn", "-task", "blastn", "-query", GENE, "-db", db,
                        "-outfmt", "6 sstart send length pident qstart qend", "-dust", "no",
                        "-word_size", "9", "-evalue", "1e-5",
                        "-max_target_seqs", "100000000", "-num_threads", TH],
                       capture_output=True, text=True)
    plus = minus = 0
    hits = []
    for L in r.stdout.splitlines():
        g = L.split("\t")
        if len(g) < 6:
            continue
        a, b = int(g[0]), int(g[1])
        qs_, qe_ = int(g[4]), int(g[5])
        glen = len(GENESEQ)
        if a <= b:
            plus += 1
            s_, e_ = a - (qs_ - 1), b + (glen - qe_)
        else:
            minus += 1
            s_, e_ = b - (glen - qe_), a + (qs_ - 1)
        hits.append((max(1, s_), min(len(seq), e_), int(g[2]),
                     (qs_ - 1) + (glen - qe_)))   # fourth = extension
    for e in (".ndb", ".nhr", ".nin", ".not", ".nsq", ".ntf", ".nto"):
        try:
            os.remove(db + e)
        except OSError:
            pass
    os.remove(fa)
    # F method (2026-09-13 audit): replaced 1kb bins with interval overlap + fragment call
    #   bin defect: the ghost end 145074 and the real end 145147 differ by 73bp, but
    #                 fell into bins 144/145 so both survived, giving a false break + a unit with its boundary shifted by 73bp
    #   the extension rule itself is kept (some real ends are one or two bases short because of an end-base variant)
    keep = []
    for h in sorted(hits):
        if keep and h[0] <= keep[-1][1]:          # ① interval overlap -> keep only the longer
            if h[2] > keep[-1][2]:
                keep[-1] = h
        else:
            keep.append(h)
    keep.sort(key=lambda x: x[1])
    i = 0
    while i + 1 < len(keep):                      # ② neighbour with gap < LO : drop only the fragment side
        if keep[i + 1][1] - keep[i][1] < LO:
            fa_, fb_ = keep[i][3] > FRAG, keep[i + 1][3] > FRAG
            if fa_ and not fb_:
                keep.pop(i)
                i = max(0, i - 1)
                continue
            if fb_ and not fa_:
                keep.pop(i + 1)
                i = max(0, i - 1)
                continue
        i += 1
    return [(a_, b_, c_) for (a_, b_, c_, d_) in keep], (minus > plus)


def call_variants(seqs, tag):
    q = os.path.join(W, "q_%s.fa" % tag)
    with io.open(q, "w", newline=NL) as g:
        for i, s in enumerate(seqs):
            g.write(">U%d%s%s%s" % (i, NL, s, NL))
    r = subprocess.run(["minimap2", "-x", PRESET, "-c", "--cs=long", "-t", TH,
                        "--secondary=no", CANREF, q], capture_output=True, text=True)
    OUT = {}
    for L in r.stdout.splitlines():
        f2 = L.split("\t")
        if len(f2) < 12:
            continue
        i = int(f2[0][1:])
        if i in OUT:
            continue
        cs = ""
        for x in f2[12:]:
            if x.startswith("cs:Z:"):
                cs = x[5:]
        qs, qe, ts, te = int(f2[2]), int(f2[3]), int(f2[7]), int(f2[8])
        toks = []
        reb = []
        tp = ts
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
                if CTX <= tp < CTX + 2240:
                    toks.append("%d%s" % (tp - CTX + 1, v[1].upper()))
                reb.append(v[1].upper())
                tp += 1
            elif op == "-":
                if CTX <= tp < CTX + 2240:
                    toks.append("%dD%d" % (tp - CTX + 1, len(v)))
                tp += len(v)
            elif op == "+":
                if CTX <= tp < CTX + 2240:
                    toks.append("%dI%d" % (tp - CTX + 1, len(v)))
                reb.append(v.upper())
        # ends are filled by direct sequence, not trusting the aligner
        # no filling (G8v4): aligned with margin to the context-padded reference, so only the alignment result is used
        OUT[i] = dict(toks=sorted(set(toks)), qs=qs, qe=qe, ts=ts, te=te,
                      reb="".join(reb), fill=0, cov_e=(2240 if te >= CTX + 2240 else 0), cov_s=(0 if ts <= CTX else 1))
    os.remove(q)
    return OUT


def verify(seq, res):
    if res is None:
        return "NOALN", "정렬 없음"
    if res["reb"] != seq[res["qs"]:res["qe"]]:
        return "SEQ_MISMATCH", "복원 %d vs 원본 %d" % (len(res["reb"]), res["qe"] - res["qs"])
    if res["cov_e"] < 2240 or res["cov_s"] > 0:
        return "NO_FULL", "정본 %d..%d (+채움 %d)" % (res["ts"], res["te"], res["fill"])
    return ("OK_FILLED" if res["fill"] else "OK"), (
        "말단 %d bp 서열로 채움" % res["fill"] if res["fill"] else "")


HDR = ["sample", "slot", "N", "array_pos", "seg_type", "unit_start", "unit_end",
       "unit_len", "gap_to_next", "chain_ok", "n_var", "n_var_gene",
       "varlist", "verify", "verify_note", "slot_space"]
MAP = io.open(os.path.join(W, "CONSENSUS_MAP_G8_%s.tsv" % TAG), "w",
              encoding="utf8", newline=NL)
MAP.write("\t".join(HDR) + NL)
SUM = io.open(os.path.join(W, "CONSENSUS_SUM_G8_%s.tsv" % TAG), "w",
              encoding="utf8", newline=NL)
SUM.write("sample\tslot\tbp\tN_units\trevcomp\toff_prox\toff_dist\tn_break\n")
STAT = collections.Counter()

for line in io.open(sys.argv[1], encoding="utf8"):
    f = line.strip().split("\t")
    if len(f) < 2:
        continue
    smp, sl = f[0], f[1]
    fn = SEQF.get((smp, sl))
    if not fn:
        STAT["정본 없음"] += 1
        continue
    p = os.path.join(os.path.expanduser(os.environ.get("CMPD", "~/cmp4d_v7")), fn)
    if not os.path.exists(p):
        STAT["파일 없음"] += 1
        sys.stderr.write("  파일 없음 %s_%s %s%s" % (smp, sl, p, NL))
        continue
    a = rd(p).upper()
    if RCF.get((smp, sl)):
        a = rc(a)
    hits, isrc = gene_hits(a, "%s_%s_%s" % (TAG, smp, sl))
    if isrc:
        # if the orientation is flipped, flip the sequence and search again (do not flip coordinates by hand)
        a = rc(a)
        FLIP.add((smp, sl))
        hits, isrc = gene_hits(a, "%s_%s_%s_r" % (TAG, smp, sl))
        if isrc:
            STAT["뒤집은 뒤에도 - 가닥 우세"] += 1
    if len(hits) < 3:
        STAT["유전자 HSP 부족"] += 1
        continue
    # unit boundary = gene end -> next gene end  (canonical phase)
    ends = [e for (s, e, ln) in hits]
    bounds = []
    nbreak = 0
    for i in range(len(ends) - 1):
        g = ends[i + 1] - ends[i]
        if LO <= g <= HI:
            bounds.append((ends[i] + 1, ends[i + 1]))     # 1-based inclusive
        else:
            nbreak += 1
    if not bounds:
        STAT["사슬 없음"] += 1
        continue
    N = len(bounds)
    seqs = [a[max(0, s - 1 - PAD):min(len(a), e + PAD)] for (s, e) in bounds]   # margin before and after
    V = call_variants(seqs, "%s_%s_%s" % (TAG, smp, sl))
    off_prox = ends[0]                       # junction ~ end of the first gene
    off_dist = len(a) - ends[-1]             # end of the last gene ~ junction
    SUM.write("%s\t%s\t%d\t%d\t%s\t%d\t%d\t%d\n"
              % (smp, sl, len(a), N, "Y" if (smp, sl) in FLIP else "N",
                 off_prox, off_dist, nbreak))
    # proximal fragment (junction ~ end of the first gene), no variant calling
    MAP.write("\t".join([smp, sl, str(N), "0", "PARTIAL_PROX", "1", str(ends[0]),
                         str(ends[0]), "0", "-", "0", "0", "-", "SKIP", "온전유닛 아님", "assembly"]) + NL)
    STAT["PARTIAL_PROX"] += 1
    for i, (s, e) in enumerate(bounds):
        res = V.get(i)
        vd, note = verify(seqs[i], res)
        STAT[vd] += 1
        gap = (bounds[i + 1][0] - s) if i + 1 < N else 0
        if vd in ("OK", "OK_FILLED"):
            toks = res["toks"]
            ng = sum(1 for t in toks if int(NUM.match(t).group()) > 2121)
            STAT["유전자변이"] += ng
            vl = ",".join(toks) if toks else "-"
            MAP.write("\t".join([smp, sl, str(N), str(i + 1), "UNIT", str(s), str(e),
                                 str(e - s + 1), str(gap), "Y", str(len(toks)),
                                 str(ng), vl, vd, note, "assembly"]) + NL)
        else:
            MAP.write("\t".join([smp, sl, str(N), str(i + 1), "UNIT", str(s), str(e),
                                 str(e - s + 1), str(gap), "Y", "0", "0", "-",
                                 vd, note, "assembly"]) + NL)
    # distal fragment (end of the last gene ~ junction)
    if off_dist > 0:
        MAP.write("\t".join([smp, sl, str(N), str(N + 1), "PARTIAL_DIST",
                             str(ends[-1] + 1), str(len(a)), str(off_dist), "0", "-",
                             "0", "0", "-", "SKIP", "온전유닛 아님", "assembly"]) + NL)
        STAT["PARTIAL_DIST"] += 1
    sys.stderr.write("  %-10s 슬롯%s 유닛 %d · 끊김 %d%s" % (smp, sl, N, nbreak, NL))
    sys.stderr.flush()
MAP.close()
SUM.close()
print("consensus map G8")
for k, v in sorted(STAT.items()):
    print("  %-18s %d" % (k, v))
