# -*- coding: utf-8 -*-
# Maps the transcribed-region table onto the unit frame used here
#   Input: 05_requests/Atlas_reply_2026-09-30/{CHM13_consensus_units.fa, transcribed_parts_for_drawing.tsv}
#   Method: the transcribed sequence is extracted and locally aligned to both strands of a doubled copy of our unit.
#          (Bio.Align.PairwiseAligner)
#   The whole consensus sequence is also aligned, to measure the phase offset between the two frames separately.
#   Coordinate frame: our canonical = the canon of P10_SUBUNIT_FRAMES.json
#            (RNU2, SST1 and DXZ4 are the identity; 5S is 1357-p reversed; D4Z4 is rotated by +336 / +826)
import io, json, os
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner

def _need(_v):
    """No defaults are set: the script stops if a variable is missing.
    The roots of the data this script reads are listed in the
    "Running this code" section of README.md."""
    import os as _o, sys as _s
    _x = _o.environ.get(_v)
    if not _x:
        _s.exit(
            "environment variable %s is not set; it is required by %s"
            % (_v, _o.path.basename(__file__))
            + "\n(a module stops the script that imports it)."
            + "\nThe deposited copy has no default; see README.md"
            + " section \"Running this code\"."
        )
    return _x


REP = os.path.join(_need("RDNA_OUT_ROOT"), "05_requests", "Atlas_reply_2026-09-30")
OURD = os.path.join(_need("RDNA_OUT_ROOT"), "02_output_tables", "transcribed_region_Atlas", "our_units")
FRJ = os.path.join(_need("RDNA_OUT_ROOT"), "02_output_tables", "P10", "P10_SUBUNIT_FRAMES.json")
OUT = os.path.join(_need("RDNA_OUT_ROOT"), "02_output_tables", "transcribed_region_Atlas")
CONS2OUR = {"RNA5S_CHM13cons": ("5S", "5S"),
            "RNU2_CHM13cons": ("RNU2", "RNU2"),
            "D4Z4_4q_CHM13cons": ("D4Z4_4q", "D4Z4_4q"),
            "D4Z4_10q_CHM13cons": ("D4Z4_10q", "D4Z4_10q"),
            "DXZ4_CHM13cons": ("DXZ4", "DXZ4"),
            "SST1c19_a_chr19_38.8Mb_CHM13cons": ("SST1_a", "SST1"),
            "SST1c19_c_chr19_39.9Mb_CHM13cons": ("SST1_c", "SST1")}
ARR2CONS = {"RNA5S": "RNA5S_CHM13cons", "RNU2": "RNU2_CHM13cons",
            "D4Z4_4q": "D4Z4_4q_CHM13cons", "D4Z4_10q": "D4Z4_10q_CHM13cons",
            "DXZ4": "DXZ4_CHM13cons",
            "SST1c19_a_chr19_38.8Mb": "SST1c19_a_chr19_38.8Mb_CHM13cons",
            "SST1c19_c_chr19_39.9Mb": "SST1c19_c_chr19_39.9Mb_CHM13cons"}


def rdfa(p):
    d, n = {}, None
    for L in io.open(p, encoding="utf8"):
        L = L.rstrip("\n")
        if L.startswith(">"):
            n = L[1:].split()[0]; d[n] = []
        elif n:
            d[n].append(L.strip())
    return {k: "".join(v).upper() for k, v in d.items()}


CONS = rdfa(os.path.join(REP, "CHM13_consensus_units.fa"))
OURS = {k: rdfa(os.path.join(OURD, k + ".fa"))[k]
        for k in ("5S", "RNU2", "D4Z4_4q", "D4Z4_10q", "DXZ4", "SST1_a", "SST1_c")}
FR = json.load(io.open(FRJ, encoding="utf8"))
AL = PairwiseAligner(mode="local", match_score=2, mismatch_score=-3,
                     open_gap_score=-5, extend_gap_score=-2)


def locate(q, unit):
    """local alignment of q to unit (doubled copy, both strands) -> (start, end, strand, score, coverage) 1-based"""
    L = len(unit)
    dbl = unit + unit
    best = None
    for strand, tgt in (("+", dbl), ("-", str(Seq(dbl).reverse_complement()))):
        a = AL.align(tgt, q)[0]
        if best is None or a.score > best[0]:
            best = (a.score, strand, a)
    score, strand, a = best
    qmap = {}
    for (ts, te), (qs, qe) in zip(*a.aligned):
        for k in range(qe - qs):
            qmap[qs + k] = ts + k
    ks = sorted(qmap)
    if not ks:
        return None
    t0, t1 = qmap[ks[0]], qmap[ks[-1]]
    if strand == "-":
        t0, t1 = 2 * L - 1 - t1, 2 * L - 1 - t0
    return t0 % L + 1, t1 % L + 1, strand, score / (2.0 * len(q)), len(ks) / float(len(q))


print("1) consensus sequence vs our unit: phase offset")
print("%-34s %-9s %6s %6s %8s %8s %7s" % ("consensus", "our unit", "cons", "ours", "cons1 위치", "strand", "score"))
for c, (ok, fam) in CONS2OUR.items():
    r = locate(CONS[c][:300], OURS[ok])
    if r:
        print("%-34s %-9s %6d %6d %8d %8s %7.3f"
              % (c, ok, len(CONS[c]), len(OURS[ok]), r[0], r[2], r[3]))
    else:
        print("%-34s %-9s not found" % (c, ok))

print()
print("2) transcribed region to our canonical coordinates")
rows = [("array_our", "feature", "grade", "cons_start", "cons_end", "cons_strand",
         "our_win_start", "our_win_end", "our_strand", "score", "coverage",
         "canon_start", "canon_end", "wraps", "note")]
print("%-9s %-24s %-5s %7s %7s %6s %7s %7s %7s %6s  %s"
      % ("array", "feature", "grade", "win_s", "win_e", "strand", "score", "cov", "can_s", "can_e", "wrap"))
with io.open(os.path.join(REP, "transcribed_parts_for_drawing.tsv"), encoding="utf8") as f:
    hdr = None
    for L in f:
        if L.startswith("#"):
            continue
        g = L.rstrip("\n").split("\t")
        if hdr is None:
            hdr = g
            continue
        arr, feat, s, e, st, grade, note = g[0], g[1], int(g[2]), int(g[3]), g[4], g[5], g[6]
        cname = ARR2CONS.get(arr)
        if cname is None:
            continue
        ok, fam = CONS2OUR[cname]
        q = CONS[cname][s - 1:e]
        r = locate(q, OURS[ok])
        if r is None:
            print("%-9s %-24s %-5s not found" % (fam, feat, grade))
            continue
        ws, we, strand, sc, cov = r
        canon = FR[fam]["canon"]
        f2 = lambda p: canon[min(max(p - 1, 0), len(canon) - 1)]
        cs, ce = f2(ws), f2(we)
        wrap = "yes" if ws > we else ""
        print("%-9s %-24s %-5s %7d %7d %6s %7.3f %6.3f %7d %7d  %s"
              % (fam, feat, grade, ws, we, strand, sc, cov, cs, ce, wrap))
        rows.append((fam, feat, grade, str(s), str(e), st, str(ws), str(we), strand,
                     "%.3f" % sc, "%.3f" % cov, str(cs), str(ce), wrap, note))
io.open(os.path.join(OUT, "TRANSCRIBED_FOR_PANEL_A.tsv"), "w", encoding="utf8",
        newline="\n").write("\n".join("\t".join(r) for r in rows) + "\n")
print()
print("tables written:", os.path.join(OUT, "TRANSCRIBED_FOR_PANEL_A.tsv"))
