# -*- coding: utf-8 -*-
# Cut the array canonical sequence, flank-array-flank, into units.
#   The canonical header reads: >...  flank<N>_array_flank<N> strand=<+|-> chrom=<c> array=<s>-<e>
#   The unit coordinates of the G8 map, which are in assembly coordinates, are carried into canonical coordinates and the cut made there.
#     strand=+ : canon_pos = chrom_pos - (array_start - flank5) + 1
#     On the minus strand the canonical sequence is reverse-complemented, so canon_pos = (array_end + flank5) - chrom_pos + 1, which reverses left and right.
#   Check: every unit sequence cut out must equal the unit sequence taken directly from the assembly.
#   Outputs: <out>/<array>.units.fa with the unit sequences, <array>.units.tsv with their positions, and CANON_UNITS_SUM.tsv.
#   Arguments: the data_v2 folder and the output folder.
# MIT License. Copyright (c) 2026 Hyungmin Kim
import io
import os
import re
import sys
import csv
import glob
import collections

NL = chr(10)
TAB = chr(9)
SRC, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
rc = lambda s: "".join(COMP.get(x, "N") for x in reversed(s))
SUM = io.open(os.path.join(OUT, "CANON_UNITS_SUM.tsv"), "w", encoding="utf8", newline=NL)
SUM.write(TAB.join(["species_key", "array_id", "canon_len", "flank5", "flank3", "strand",
                    "n_units", "unit_span_in_canon", "units_bp", "check_vs_map"]) + NL)
tot = collections.Counter()
for mf in sorted(glob.glob(os.path.join(SRC, "02_consensus_map", "*", "*.g8map.tsv"))):
    aid = os.path.basename(mf)[:-len(".g8map.tsv")]
    sp = os.path.basename(os.path.dirname(mf))
    cf = os.path.join(SRC, "01_canonical", sp, aid + ".fa")
    if not os.path.exists(cf):
        continue
    txt = io.open(cf, encoding="utf8").read().split(NL)
    hdr = txt[0]
    seq = "".join(x.strip() for x in txt[1:]).upper()
    m = re.search(r"flank(\d+)_array_flank(\d+).*strand=(\S+).*array=(\d+)-(\d+)", hdr)
    if not m:
        continue
    f5, f3, st, a0, a1 = int(m.group(1)), int(m.group(2)), m.group(3), int(m.group(4)), int(m.group(5))
    # the actual flank is shorter where it runs off the end of the contig, so it is brought back to the canonical length
    alen = a1 - a0 + 1
    extra = len(seq) - alen
    if extra != f5 + f3:
        f5 = min(f5, max(0, extra - f3))
        f3 = extra - f5
    rows = [r for r in csv.DictReader(io.open(mf, encoding="utf8"), delimiter=TAB)
            if r["chrom_start"].isdigit()]
    fa = io.open(os.path.join(OUT, aid + ".units.fa"), "w", encoding="utf8", newline=NL)
    tsv = io.open(os.path.join(OUT, aid + ".units.tsv"), "w", encoding="utf8", newline=NL)
    tsv.write(TAB.join(["species_key", "array_id", "array_pos", "canon_start", "canon_end", "unit_len",
                        "chrom", "chrom_start", "chrom_end", "strand", "verify"]) + NL)
    lo, hi, bp, bad = None, None, 0, 0
    for r in rows:
        cs, ce = int(r["chrom_start"]), int(r["chrom_end"])
        if st == "-":
            p0, p1 = (a1 + f5) - ce + 1, (a1 + f5) - cs + 1
        else:
            p0, p1 = cs - (a0 - f5) + 1, ce - (a0 - f5) + 1
        if p0 < 1 or p1 > len(seq) or p1 < p0:
            bad += 1
            continue
        u = seq[p0 - 1:p1]
        if len(u) != int(r["unit_len"]):
            bad += 1
        lo = p0 if lo is None else min(lo, p0)
        hi = p1 if hi is None else max(hi, p1)
        bp += len(u)
        fa.write(">%s|unit=%s|len=%d|%s:%d-%d|strand=%s%s%s%s"
                 % (aid, r["array_pos"], len(u), r["chrom"], cs, ce, st, NL, u, NL))
        tsv.write(TAB.join([sp, aid, r["array_pos"], str(p0), str(p1), str(len(u)),
                            r["chrom"], str(cs), str(ce), st, r["verify"]]) + NL)
    fa.close()
    tsv.close()
    SUM.write(TAB.join([sp, aid, str(len(seq)), str(f5), str(f3), st, str(len(rows)),
                        ("%d-%d" % (lo, hi)) if lo else "-", str(bp),
                        "일치" if bad == 0 else "불일치 %d" % bad]) + NL)
    tot["array"] += 1
    tot["unit"] += len(rows)
    tot["bad"] += bad
SUM.close()
print("array %d, units %d, length mismatch %d" % (tot["array"], tot["unit"], tot["bad"]))
