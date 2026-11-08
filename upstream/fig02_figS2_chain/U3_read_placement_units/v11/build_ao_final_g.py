# -*- coding: utf-8 -*-
# ARRAY_ONLY G8v3 final unit table (stage before deletion marking), local, only merges tables (2026-09-12)
#   base = 3rd-round cut (window alignment transfers canonical G8v3 boundaries onto the read, same as for junction), number = alignment position
#     - the 20 relabeled reads are reverted (the number is the real position. Forcing stored values is wrong, 66-sample measurement 48 vs 19)
#   supplement = blastn cut pieces that do not overlap 3rd-round pieces: numbered by the number difference to the nearest paired piece in the read,
#          and included if within the 3rd-round number range of that read (bnd = BLASTN). Read-end pieces outside the range are set aside (decision pending)
#   flag = reads where 3rd-round and blastn numbers disagree -> numbering column after the 17 columns 'UNCERTAIN' (others 'OK')
#   Output: per node AO_G8_UNITS_<node>.final0.tsv (same 17 columns as the 3rd-round table plus numbering), AO_READEND_EXTRA.tsv, summary
#   Arguments: cmp folder, local_blast folder
import io
import os
import re
import sys
import glob
import collections

NL = chr(10)
TAB = chr(9)
C, LB = sys.argv[1], sys.argv[2]
NT = {}
for p in glob.glob(os.path.join(C, "g8cut_mid", "run3", "AO_G8_READS_*.raw.tsv")):
    for L in io.open(p, encoding="utf8"):
        r = L.rstrip(NL).split(TAB)
        if len(r) > 10:
            NT[(r[0], r[1])] = r[10]


def reld(k):
    m = re.search(r"relabeled_d=([+-]\d+)", NT.get(k, ""))
    return int(m.group(1)) if m else 0


B = collections.defaultdict(list)
for L in io.open(os.path.join(LB, "full", "AO_BLAST_UNITS.tsv"), encoding="utf8"):
    r = L.rstrip(NL).split(TAB)
    if r[7].lstrip("-").isdigit() and r[8].lstrip("-").isdigit():
        B[(r[0], r[2])].append(r)
S = collections.Counter()
EX = io.open(os.path.join(C, "g8cut_mid", "final", "AO_READEND_EXTRA.tsv"), "w", encoding="utf8", newline=NL)
for p in sorted(glob.glob(os.path.join(C, "g8cut_mid", "run3", "AO_G8_UNITS_*.raw.tsv"))):
    host = os.path.basename(p).split("_")[3].replace(".raw.tsv", "")
    host = os.path.basename(p)[len("AO_G8_UNITS_"):-len(".raw.tsv")]
    T = collections.defaultdict(list)
    order = []
    for L in io.open(p, encoding="utf8"):
        r = L.rstrip(NL).split(TAB)
        k = (r[0], r[2])
        if k not in T:
            order.append(k)
        d = reld(k)
        if d and r[7].lstrip("-").isdigit() and r[8].lstrip("-").isdigit():
            r[7] = str(int(r[7]) - d)
            r[8] = str(int(r[8]) - d)
            S["번호옮기기 되돌린 유닛"] += 1
        T[k].append(r)
    O = io.open(os.path.join(C, "g8cut_mid", "final", "AO_G8_UNITS_%s.final0.tsv" % host), "w", encoding="utf8", newline=NL)
    for k in order:
        t = T[k]
        b = sorted(B.get(k, []), key=lambda x: int(x[9]))
        # pair: coordinates within 3bp
        off = {}
        unc = False
        for i, x in enumerate(b):
            f0, f1 = int(x[9]), int(x[10])
            m = [y for y in t if y[8].lstrip("-").isdigit() and abs(int(y[9]) - f0) <= 3 and abs(int(y[10]) - f1) <= 3]   # repaired spans have no G6 counterpart, so g6=None -> excluded from pair candidates (old data had 0 such rows, so results unchanged)
            if m:
                off[i] = int(m[0][8]) - int(x[8])
                if int(m[0][8]) != int(x[8]):
                    unc = True
        nums = [int(y[8]) for y in t if y[8].lstrip("-").isdigit()]
        lo, hi = (min(nums), max(nums)) if nums else (None, None)
        have = set(nums)
        tag = "UNCERTAIN" if unc else "OK"
        S["read"] += 1
        S["read 번호 불확실"] += unc
        for y in t:
            O.write(TAB.join(y[:17] + [tag]) + NL)
            S["3차 유닛"] += 1
        for i, x in enumerate(b):
            if i in off:
                continue
            f0, f1 = int(x[9]), int(x[10])
            if any(min(f1, int(y[10])) - max(f0, int(y[9])) > 0.5 * (f1 - f0) for y in t):
                S["blastn 조각 : 3차 조각과 겹침(3차 경계 씀)"] += 1
                continue
            if not off:
                S["blastn 조각 : 짝 없는 read(번호 근거 없음 · 뺌)"] += 1
                continue
            j = min(off, key=lambda z: abs(z - i))
            g6 = int(x[8]) + off[j]
            g8 = int(x[7]) + off[j]
            row = [x[0], x[1], x[2], x[3], x[4], x[5], x[6], str(g8), str(g6), x[9], x[10], "BLASTN", "BLASTN", x[11], x[12], "-", "-"]
            if lo is None or g6 < lo or g6 > hi:
                EX.write(TAB.join(row + [tag]) + NL)
                S["blastn 조각 : 3차 번호 범위 밖(read 끝 · 따로 뗌)"] += 1
                continue
            if g6 in have:
                S["blastn 조각 : 번호가 3차와 겹침(뺌)"] += 1
                continue
            have.add(g6)
            O.write(TAB.join(row + [tag]) + NL)
            S["blastn 조각 : 범위 안 보충"] += 1
    O.close()
EX.close()
for x, v in sorted(S.items()):
    print("  %-44s %d" % (x, v))
