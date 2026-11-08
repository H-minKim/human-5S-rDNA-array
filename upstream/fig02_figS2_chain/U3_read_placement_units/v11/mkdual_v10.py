# -*- coding: utf-8 -*-
# Regenerate the DUAL list, based on the new ARRAY_ONLY placement (2026-09-13)
#   DUAL = reads accepted as ARRAY_ONLY that are also in the junction assignment table (HAP_READ_FINAL).
#   They are removed from the ARRAY_ONLY side so that the same read is not counted twice in the combined set.
#   Same rule as the old dual*.py: union of HAP_READ_FINAL.tsv and .pairingspace.tsv, rows whose 5th field is not "-".
#   What changed: the AO read list is taken directly from the new AOJV table instead of the old aoreads*.txt.
#   The original DUAL is not touched. Output goes to a new folder.
import io
import os

H = "/home/ubuntu/"
NL = chr(10)
TAB = chr(9)
OUT = H + "dual_v10"
os.makedirs(OUT, exist_ok=True)


def load_hap(p):
    S = set()
    if not os.path.exists(p):
        return S
    for i, L in enumerate(io.open(p, encoding="utf8")):
        if i == 0:
            continue
        r = L.rstrip(NL).split(TAB)
        if len(r) >= 7 and r[4] != "-":
            S.add((r[0], r[1][:16]))
    return S


J = load_hap(H + "HAP_READ_FINAL.pairingspace.tsv") | load_hap(H + "HAP_READ_FINAL.tsv")
print("junction 배정표 read %d" % len(J))
for src, fn in (("DUAL20", "AOJV3.tsv"), ("DUAL5", "AOJV5.tsv"),
                ("DUAL8", "AOJV8.tsv"), ("DUALG", "AOJV_GIAB.tsv")):
    AO = set()
    p = H + "repro_in_v10/" + fn
    for L in io.open(p, encoding="utf8"):
        f = L.rstrip(NL).split(TAB)
        if len(f) >= 2:
            AO.add((f[0], f[1][:16]))
    both = AO & J
    o = io.open(os.path.join(OUT, src + ".tsv"), "w", encoding="utf8", newline=NL)
    for k in sorted(both):
        o.write(k[0] + TAB + k[1] + NL)
    o.close()
    old = H + "aog8/dual/" + src + ".tsv"
    no = sum(1 for _ in io.open(old, encoding="utf8")) if os.path.exists(old) else -1
    print("%-7s : AO read %5d · DUAL %4d (옛 %4d · 차이 %+d)" % (src, len(AO), len(both), no, len(both) - no))
