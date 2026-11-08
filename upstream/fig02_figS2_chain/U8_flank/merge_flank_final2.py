# -*- coding: utf-8 -*-
# Final flank merge v2, reassembly first + quality gate
#
#   Priority : hifiasm reassembly (v12) > v11 (HPRC assembly, read consensus)
#   But if the reassembled flank is bad, v11 is used. Failure modes confirmed by measurement :
#     HG02140 h1 3' - reassembly is 87.82% to the reference with 1,108bp overlap, read consensus is 99.42%/5,002bp
#     HG03902 h1 3' - 88.00%/1,100bp vs 99.44%/5,002bp     (the assembly continues into a different sequence at that spot)
#     NA19118 h2 5' - reassembly only 1,380bp
#   Gate (separately for each side)
#     - length >= MINLEN(5,000)
#     - overlap with the reference (CHM13 flank 5kb) >= MINBP(4,000), a short overlap means it continued elsewhere
#     - if identity is lower than v11 by MARGIN(2.0%p) or more, v11 is used
#     Slots that do not align to the reference at all (0 on both sides) cannot use the reference as a yardstick, so only length is checked.
#   Invariant : the array canonical sequence is not changed by a single character. The middle is compared character by character and discarded if different.
#   Each side is chosen separately : 5' from reassembly and 3' from read consensus can be mixed.
#
#   Output : ~/flank_final2/<canonical file name>, ~/FLANK_FINAL2.tsv
import io
import os
import subprocess
import collections

NL = chr(10)
TAB = chr(9)
H = os.path.expanduser("~")
V11 = os.path.join(H, "flank_v11")
V12 = os.path.join(H, "flank_v12")
CMPD = os.path.join(H, "cmp4d_v11")
REF = os.path.join(H, "s5s_tuc", "ref")
OUT = os.path.join(H, "flank_final2")
TMP = os.path.join(H, "mf2tmp")
os.makedirs(OUT, exist_ok=True)
os.makedirs(TMP, exist_ok=True)
MINLEN = 5000
MINBP = 4000
MARGIN = 2.0
CMPX = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def rc(s):
    return "".join(CMPX.get(c, "N") for c in reversed(s))


def rd(p):
    return "".join(x.strip() for x in io.open(p, encoding="utf8") if not x.startswith(">")).upper()


def hd(p):
    h = io.open(p, encoding="utf8").readline().rstrip(NL)
    return dict(x.split("=", 1) for x in h[1:].split("|") if "=" in x), h


def ident(a, b, tag):
    fa = os.path.join(TMP, "a_%s.fa" % tag)
    fb = os.path.join(TMP, "b_%s.fa" % tag)
    io.open(fa, "w", newline=NL).write(">a" + NL + a + NL)
    io.open(fb, "w", newline=NL).write(">b" + NL + b + NL)
    r = subprocess.run(["minimap2", "-x", "asm20", "-c", "--secondary=no", "-t", "2", fa, fb],
                       capture_output=True, text=True)
    best = (0, 0)
    for L in r.stdout.splitlines():
        f = L.split(TAB)
        if len(f) < 12:
            continue
        nm, bl = int(f[9]), int(f[10])
        if bl > best[1]:
            best = (nm, bl)
    for x in (fa, fb):
        try:
            os.remove(x)
        except OSError:
            pass
    return (100.0 * best[0] / best[1] if best[1] else 0.0), best[1]


FD = rd(os.path.join(REF, "FLANK_distal_50000.fa"))
FP = rd(os.path.join(REF, "FLANK_proximal_50000.fa"))
SEQF = {}
for L in io.open(os.path.join(H, "hap_seq_node_v8.tsv"), encoding="utf8"):
    r = L.rstrip(NL).split(TAB)
    if len(r) >= 3:
        SEQF[(r[0], r[1])] = r[2]
FLIP = {}
f = io.open(os.path.join(H, "g8canonv10", "CONSENSUS_SUM_G8_full.tsv"), encoding="utf8")
f.readline()
for L in f:
    r = L.rstrip(NL).split(TAB)
    FLIP[(r[0], r[1])] = (r[4] == "Y")

LOG = io.open(os.path.join(H, "FLANK_FINAL2.tsv"), "w", encoding="utf8", newline=NL)
LOG.write(TAB.join(["sample", "slot", "file", "src5", "len5", "id5", "bp5",
                    "src3", "len3", "id3", "bp3", "canon_bp", "mid_exact"]) + NL)
st = collections.Counter()
for (smp, slot), fn in sorted(SEQF.items()):
    canon = rd(os.path.join(CMPD, fn))
    fl = FLIP.get((smp, slot))
    cand = {}
    for tag, d in (("REASM", V12), ("V11", V11)):
        p = os.path.join(d, fn)
        if not os.path.exists(p):
            continue
        m, h = hd(p)
        l5, l3 = int(m.get("flank5", 0)), int(m.get("flank3", 0))
        full = rd(p)
        mid = full[l5:len(full) - l3] if l3 else full[l5:]
        if mid != canon:
            st["%s 가운데 불일치(버림)" % tag] += 1
            continue
        src = tag if tag == "REASM" else ("HPRC" if "HPRC" in h else "READCONS")
        cand[tag] = (src, full[:l5], full[len(full) - l3:] if l3 else "")
    if not cand:
        st["없음"] += 1
        LOG.write(TAB.join([smp, slot, fn, "-", "0", "-", "0", "-", "0", "-", "0",
                            str(len(canon)), "-"]) + NL)
        continue

    def score(side, seq):
        """Measure all four reference options (FD/FP x forward/reverse complement) and use the best.
           Inferring the reference correspondence from FLIP is wrong (2026-09-15 measurement: 514/515 misjudged)."""
        if len(seq) < MINLEN:
            return None
        q = seq[-5000:] if side == "5" else seq[:5000]
        best = None
        for nm, r0 in (("FD", FD[:5000]), ("FP", FP[-5000:])):
            for ori, rr in (("+", r0), ("-", rc(r0))):
                v, bp = ident(rr, q, "%s%s%s%s%s" % (smp, slot, side, nm, ori))
                if best is None or bp > best[1]:
                    best = (v, bp)
        return best

    picked = {}
    for side, idx in (("5", 1), ("3", 2)):
        opts = []
        for tag in ("REASM", "V11"):
            if tag not in cand:
                continue
            seq = cand[tag][idx]
            sc = score(side, seq)
            opts.append((tag, cand[tag][0], seq, sc))
        ok = [o for o in opts if o[3] and o[3][1] >= MINBP]
        if ok:
            # Reassembly first. But if it is lower than v11 by MARGIN or more, use v11.
            best = ok[0]
            if len(ok) > 1 and ok[0][0] == "REASM" and ok[1][3][0] - ok[0][3][0] >= MARGIN:
                best = ok[1]
                st["관문으로 교체 %s" % side] += 1
            picked[side] = best
        else:
            # Spot where the reference cannot serve as a yardstick, use the longer one
            longest = max(opts, key=lambda o: len(o[2])) if opts else None
            if longest and len(longest[2]) >= MINLEN:
                picked[side] = longest
                st["참조무관·길이선택 %s" % side] += 1
    f5 = picked.get("5")
    f3 = picked.get("3")
    s5 = f5[2] if f5 else ""
    s3 = f3[2] if f3 else ""
    new = s5 + canon + s3
    assert new[len(s5):len(s5) + len(canon)] == canon, "정본 손상"
    if not s5 and not s3:
        st["없음"] += 1
        continue
    io.open(os.path.join(OUT, fn), "w", newline=NL).write(
        ">%s|flank5=%d|flank3=%d|src5=%s|src3=%s%s%s%s"
        % (fn.replace(".fa", ""), len(s5), len(s3),
           f5[1] if f5 else "-", f3[1] if f3 else "-", NL, new, NL))
    st["채움"] += 1
    st["5:" + (f5[1] if f5 else "-")] += 1
    st["3:" + (f3[1] if f3 else "-")] += 1
    LOG.write(TAB.join([smp, slot, fn,
                        f5[1] if f5 else "-", str(len(s5)),
                        "%.2f" % f5[3][0] if (f5 and f5[3]) else "-", str(f5[3][1]) if (f5 and f5[3]) else "0",
                        f3[1] if f3 else "-", str(len(s3)),
                        "%.2f" % f3[3][0] if (f3 and f3[3]) else "-", str(f3[3][1]) if (f3 and f3[3]) else "0",
                        str(len(canon)), "Y"]) + NL)
    LOG.flush()
LOG.close()
print("슬롯 %d · 채움 %d · 없음 %d" % (len(SEQF), st["채움"], st["없음"]))
for k in sorted(st):
    if k.startswith(("5:", "3:", "관문", "참조무관")) or "불일치" in k:
        print("  %-24s %d" % (k, st[k]))
