# -*- coding: utf-8 -*-
# Cut ARRAY_ONLY reads into G8 units without alignment (local WSL, 2026-09-12)
#   ① raw reads : the saved fastq (aoreads/*.fastq.gz, raw 4-line check passed 12,862/12,862)
#   ② unit boundaries : find the 5S gene (119bp) in the read sequence with blastn, same settings G8v3 used to cut the canonical
#                 (-task blastn -word_size 9 -evalue 1e-5 -dust no, HSP extended to the full gene length)
#                 G8 unit = base after the previous gene end ~ next gene end ([spacer p][gene p+1])
#   ③ variants : same method as the G8v3 map and ao_g8cut (minimap2 asm20 --cs=long, end fill, reconstruction check)
#   ④ numbering : actual position, add the start number to the unit index within the read (long stretches with a missed gene are skipped by counting units from length).
#            (2026-09-12: the saved placement wrongly renumbered units after those dropped within a read, so use the actual position)
#            start number and direction = the best-matching position when fitting the read's unit variant chain to the G8v3 canonical map (that slot) (same principle as the original place(),
#            except missing units are also counted in the index). Slot and range come from the saved placement (±5 units). Agreement with the saved variant list is written only as a reference column
#   arguments : work folder (local_blast), cmp folder (AOJV*.tsv, maps), number of parallel jobs
import io
import os
import re
import sys
import glob
import gzip
import subprocess
import collections
import multiprocessing

NL = chr(10)
TAB = chr(9)
W = sys.argv[1]
CMPD = sys.argv[2]
NP = int(sys.argv[3]) if len(sys.argv) > 3 else 4
GENE = os.path.join(W, "ref", "5S_gene_query.fa")
CANON = os.path.join(W, "ref", "unit_canon.fa")
CMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
CS = re.compile(r"([-:*+=])([A-Za-z0-9]+)")
NUM = re.compile(r"[0-9]+")
LONGU = 3300      # if gene end to gene end is longer than this, a gene in between was missed (unit 2,240bp, inside the canonical unit length distribution)
SHORTU = 1500     # flag if shorter than this


def rd(p):
    return "".join(x.strip() for x in io.open(p, encoding="utf8").read().split(NL)[1:])


def rc(s):
    return s.translate(CMP)[::-1]


CANSEQ = rd(CANON).upper()
GLEN = len(rd(GENE))


def call_variants(seqs, tag, wd):
    q = os.path.join(wd, "q_%s.fa" % tag)
    with io.open(q, "w", newline=NL) as g:
        for i, s in enumerate(seqs):
            g.write(">U%d%s%s%s" % (i, NL, s, NL))
    r = subprocess.run(["minimap2", "-x", "asm20", "-c", "--cs=long", "-t", "1", "--secondary=no", CANON, q],
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
        qs, qe, ts, te = int(f2[2]), int(f2[3]), int(f2[7]), int(f2[8])
        toks, reb, tp = [], [], ts
        for m in CS.finditer(cs):
            op, v = m.group(1), m.group(2)
            if op == "=":
                reb.append(v.upper())
                tp += len(v)
            elif op == ":":
                n = int(v)
                reb.append(CANSEQ[tp:tp + n])
                tp += n
            elif op == "*":
                toks.append("%d%s" % (tp + 1, v[1].upper()))
                reb.append(v[1].upper())
                tp += 1
            elif op == "-":
                toks.append("%dD%d" % (tp + 1, len(v)))
                tp += len(v)
            elif op == "+":
                toks.append("%dI%d" % (tp + 1, len(v)))
                reb.append(v.upper())
        qseq = seqs[i]
        kt = 0
        while te + kt < 2240 and qe + kt < len(qseq):
            if qseq[qe + kt] != CANSEQ[te + kt]:
                toks.append("%d%s" % (te + kt + 1, qseq[qe + kt]))
            kt += 1
        if te + kt < 2240:
            toks.append("%dD%d" % (te + kt + 1, 2240 - (te + kt)))
        kh = 0
        while ts - 1 - kh >= 0 and qs - 1 - kh >= 0:
            if qseq[qs - 1 - kh] != CANSEQ[ts - 1 - kh]:
                toks.append("%d%s" % (ts - kh, qseq[qs - 1 - kh]))
            kh += 1
        OUT[i] = dict(toks=toks, qs=qs, qe=qe, reb="".join(reb), fill=kh + kt, cov_e=te + kt, cov_s=ts - kh)
    os.remove(q)
    return OUT


def verify(seq, res):
    if res is None:
        return "NOALN"
    if res["reb"] != seq[res["qs"]:res["qe"]]:
        return "SEQ_MISMATCH"
    if res["cov_s"] != 0 or res["cov_e"] != 2240:
        return "NO_FULL"
    return "OK_FILLED" if res["fill"] else "OK"


def spacer(toks):
    return set(t for t in toks if 0 < int(NUM.match(t).group()) <= 2121)


def jac(a, b):
    un = a | b
    return len(a & b) / float(len(un)) if un else 1.0


def init():
    global SAVED, CFGF, G6TO8
    SAVED = collections.defaultdict(dict)       # (cfg, smp, r16) -> (slot, {g6: tokens})
    for cfg, fn in (("AOJV3", "AOJV3.tsv"), ("AOJV5", "AOJV5.tsv"), ("AOJV8", "AOJV8.tsv"), ("AOJVG", "AOJV_GIAB.tsv")):
        for L in io.open(os.path.join(CMPD, fn), encoding="utf8"):
            r = L.rstrip(NL).split(TAB)
            if len(r) < 5:
                continue
            e = SAVED[(cfg, r[0])].setdefault(r[1], [r[2], {}])
            e[1][int(r[3])] = set() if r[4] == "-" else set(r[4].split(","))
    # G6 number -> G8 number : as in ao_g8cut, the G6 unit containing the middle of the G8 unit spacer (start+1000)
    def units(p, types):
        U = collections.defaultdict(list)
        f = io.open(p, encoding="utf8")
        ix = dict((n, i) for i, n in enumerate(f.readline().rstrip(NL).split(TAB)))
        for L in f:
            r = L.rstrip(NL).split(TAB)
            if r[ix["seg_type"]] in types:
                U[(r[0], r[1])].append((int(r[ix["unit_start"]]), int(r[ix["unit_end"]]), int(r[ix["array_pos"]])))
        return U
    U6 = units(os.path.join(CMPD, "CONSENSUS_MAP_G6C2.tsv"), ("UNIT", "UNIT_LAST"))
    U8 = units(os.path.join(CMPD, "g8v3", "CONSENSUS_MAP_G8_full.tsv"), ("UNIT",))
    global MAP8, MALL, G8TO6
    MAP8 = collections.defaultdict(dict)
    f = io.open(os.path.join(CMPD, "g8v3", "CONSENSUS_MAP_G8_full.tsv"), encoding="utf8")
    ix = dict((n, i) for i, n in enumerate(f.readline().rstrip(NL).split(TAB)))
    for L in f:
        r = L.rstrip(NL).split(TAB)
        if r[ix["seg_type"]] == "UNIT" and r[ix["verify"]] in ("OK", "OK_FILLED"):
            v = r[ix["varlist"]]
            MAP8[(r[0], r[1])][int(r[ix["array_pos"]])] = set() if v == "-" else set(v.split(","))
    MALL = dict((k, set().union(*v.values())) for k, v in MAP8.items())
    G6TO8 = {}
    for k, u8 in U8.items():
        u6 = sorted(U6.get(k, []))
        m = {}
        for s8, e8, p8 in u8:
            mid = s8 + 1000
            for s6, e6, p6 in u6:
                if s6 <= mid <= e6:
                    m[p6] = p8
                    break
        G6TO8[k] = m
    G8TO6 = dict((k, dict((b, a) for a, b in m.items())) for k, m in G6TO8.items())


def work(fp):
    smp, cfg, lab = os.path.basename(fp).split(".")[:3]
    wd = os.path.join(W, "tmp", smp + "." + cfg)
    os.makedirs(wd, exist_ok=True)
    sv = SAVED[(cfg, smp)]
    R = {}
    with gzip.open(fp, "rt", encoding="utf8") as fh:
        rec = []
        for x in fh:
            rec.append(x.rstrip(NL))
            if len(rec) == 4:
                rid = rec[0][1:].split()[0]
                R[rid] = rec[1]
                rec = []
    fa = os.path.join(wd, "reads.fa")
    with io.open(fa, "w", newline=NL) as g:
        for rid, s in R.items():
            g.write(">%s%s%s%s" % (rid, NL, s, NL))
    subprocess.run(["makeblastdb", "-in", fa, "-dbtype", "nucl", "-out", os.path.join(wd, "db")], capture_output=True)
    b = subprocess.run(["blastn", "-task", "blastn", "-query", GENE, "-db", os.path.join(wd, "db"),
                        "-outfmt", "6 sseqid sstart send length pident qstart qend", "-dust", "no",
                        "-word_size", "9", "-evalue", "1e-5", "-max_target_seqs", "100000", "-max_hsps", "100000"],
                       capture_output=True, text=True)
    HS = collections.defaultdict(list)
    for L in b.stdout.splitlines():
        g = L.split(TAB)
        a, e, qs_, qe_ = int(g[1]), int(g[2]), int(g[5]), int(g[6])
        if a <= e:
            HS[g[0]].append(("+", a - (qs_ - 1), e + (GLEN - qe_), int(g[3])))
        else:
            HS[g[0]].append(("-", e - (GLEN - qe_), a + (qs_ - 1), int(g[3])))
    out_u, out_r, seqs, meta = [], [], [], []
    for rid, s in sorted(R.items()):
        r16 = rid[:16]
        if r16 not in sv:
            continue
        slot, P = sv[r16]
        L0 = len(s)
        hs = HS.get(rid, [])
        np_ = sum(1 for h in hs if h[0] == "+")
        nm_ = len(hs) - np_
        st = "+" if np_ >= nm_ else "-"
        # canonical (G8) direction = the direction in which genes read on the + strand. If the - strand is the majority, treat the read as reverse complement and convert coordinates
        gi = []
        for h in hs:
            if h[0] != st:
                continue
            a, e = max(1, h[1]), min(L0, h[2])
            if st == "-":
                a, e = L0 - e + 1, L0 - a + 1
            gi.append((a, e, h[3]))
        gi.sort()
        dd = []                                   # for overlapping HSPs keep the longest one
        for x in gi:
            if dd and x[0] <= dd[-1][1]:
                if x[2] > dd[-1][2]:
                    dd[-1] = x
            else:
                dd.append(x)
        o = s if st == "+" else rc(s)
        ku = []
        kk = 0                                    # actual unit index within the read (where a missed gene merged units, skip by counting units from length)
        anom = collections.Counter()
        for i in range(len(dd) - 1):
            e0, e1 = dd[i][1], dd[i + 1][1]      # 1-based gene end
            ln = e1 - e0
            if ln > LONGU:                        # one or more genes missed -> do not cut the units in between, only advance the index
                kk += max(2, int(round(ln / 2240.0)))
                anom["긴구간(유전자못찾음)"] += 1
                continue
            if ln < SHORTU:                       # short stretch caught on a non-gene site -> only flag it, index advances by one
                anom["짧은구간"] += 1
            useq = o[e0:e1]                       # base after the previous gene end ~ next gene end
            f0, f1 = (e0, e1) if st == "+" else (L0 - e1, L0 - e0)
            ku.append((kk, len(seqs)))
            seqs.append(useq)
            meta.append([smp, rid, r16, slot, cfg, lab, "-" if st == "+" else "+rc", kk, f0, f1, len(useq)])
            kk += 1
        out_r.append([smp, r16, rid, slot, cfg, lab, str(len(dd)), str(len(ku)), st, ku, ";".join("%s=%d" % x for x in sorted(anom.items())) or "-"])
    V = call_variants(seqs, smp, wd) if seqs else {}
    CALL = []
    for i, u in enumerate(seqs):
        res = V.get(i)
        vd = verify(u, res)
        CALL.append((vd, sorted(set(res["toks"])) if (res and vd in ("OK", "OK_FILLED")) else None))
    RO, UO = [], []
    for rr in out_r:
        smp_, r16, rid, slot, cfg_, lab_, ng, nu, st, ku, anom = rr
        P = sv[r16][1]
        m86 = G6TO8.get((smp_, slot), {})
        m68 = G8TO6.get((smp_, slot), {})
        M8 = MAP8.get((smp_, slot), {})
        ma = MALL.get((smp_, slot), set())
        items = [(k, set(CALL[j][1])) for k, j in ku if CALL[j][1] is not None]
        p8s = [m86.get(g, g) for g in P]
        best = []
        if M8 and items:
            for dr in (1, -1):
                for base in range(min(p8s) - max(k for k, _ in items) - 5, max(p8s) + max(k for k, _ in items) + 6):
                    v = [jac(t & ma, M8[base + dr * k]) for k, t in items if (base + dr * k) in M8]
                    if len(v) >= max(2, len(items) // 2):
                        best.append((sum(v) / len(v), len(v), dr == 1, base, dr))
        best.sort(reverse=True)
        if not best:
            RO.append([smp_, r16, rid, slot, cfg_, lab_, ng, str(len(ku)), st, "-", "-", "-", "-", "NONUMBER", anom, "-"])
            continue
        s1, n1, _, base, dr = best[0]
        s2 = max([x[0] for x in best[1:] if (x[3], x[4]) != (base, dr)] or [-1])
        sv_sc = []
        for k, j in ku:
            g8 = base + dr * k
            g6 = m68.get(g8, g8)
            mt = meta[j]
            vd, toks = CALL[j]
            sc = "-"
            if toks is not None and g6 in P:
                sc = "%.3f" % jac(spacer(toks), spacer(P[g6]))
                sv_sc.append(float(sc))
            UO.append([mt[0], mt[1], mt[2], mt[3], mt[4], mt[5], mt[6], str(g8), str(g6), str(mt[8]), str(mt[9]),
                       vd, (",".join(toks) or "-") if toks is not None else "-", sc, "Y" if g6 in P else "N", str(mt[10])])
        RO.append([smp_, r16, rid, slot, cfg_, lab_, ng, str(len(ku)), st, str(dr), str(base), "%.3f" % s1, "%.3f" % s2, "OK", anom,
                   "%.3f" % (sum(sv_sc) / len(sv_sc)) if sv_sc else "-"])
    for fn in glob.glob(os.path.join(wd, "*")):
        os.remove(fn)
    os.rmdir(wd)
    return smp, cfg, UO, RO


if __name__ == "__main__":
    init()
    files = sorted(glob.glob(os.path.join(W, "aoreads", "*.fastq.gz")))
    OU = io.open(os.path.join(W, "AO_BLAST_UNITS.tsv"), "w", encoding="utf8", newline=NL)
    OR = io.open(os.path.join(W, "AO_BLAST_READS.tsv"), "w", encoding="utf8", newline=NL)
    with multiprocessing.Pool(NP, initializer=init) as pool:
        for smp, cfg, uu, rr in pool.imap_unordered(work, files):
            for x in uu:
                OU.write(TAB.join(x) + NL)
            for x in rr:
                OR.write(TAB.join(x) + NL)
            OU.flush()
            OR.flush()
            sys.stderr.write("  %s %s 유닛 %d · read %d%s" % (smp, cfg, len(uu), len(rr), NL))
    OU.close()
    OR.close()
    print("done")
