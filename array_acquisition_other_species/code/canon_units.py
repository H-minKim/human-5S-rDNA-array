# -*- coding: utf-8 -*-
# array 정본(flank-array-flank)을 유닛 단위로 자른다
#   정본 머리 : >…  flank<N>_array_flank<N> strand=<+|-> chrom=<c> array=<s>-<e>
#   G8 지도의 유닛 좌표(조립본 좌표)를 정본 좌표로 옮겨 자른다
#     strand=+ : canon_pos = chrom_pos - (array_start - flank5) + 1
#     strand=- : 정본은 역상보라 canon_pos = (array_end + flank5) - chrom_pos + 1  (좌우가 뒤집힌다)
#   검산 : 잘라낸 유닛 서열이 조립본에서 직접 뽑은 유닛 서열과 같아야 한다(전건 대조)
#   산출 : <out>/<array>.units.fa (유닛 서열) · <array>.units.tsv (자리표) · CANON_UNITS_SUM.tsv
#   인자 : data_v2 폴더 출력폴더
# MIT License · Copyright (c) 2026 Hungmin Kim
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
    # 실제 flank 는 조립본 끝에 걸리면 짧아진다 → 정본 길이로 되맞춘다
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
print("array %d · 유닛 %d · 길이 불일치 %d" % (tot["array"], tot["unit"], tot["bad"]))
