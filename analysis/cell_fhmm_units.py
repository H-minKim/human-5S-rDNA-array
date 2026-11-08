# -*- coding: utf-8 -*-
"""HepG2 FiberHMM calls.bam -> 유닛 안 TSS 기준 rel 좌표 (2026-10-12)

자리표(어느 유닛인가)가 필요 없다 — read 안에서 유전자 닻을 찾아 유닛을 가른다.
좌표는 전부 분자(molecular) 좌표다 — FiberHMM @PG 가 그렇게 적어 놨다
  (coord=molecular · ns/nl/as/al/MA in molecular original-fiber coordinates).
rel = 0 은 HPRC 와 같은 규칙 : 유닛 안 마지막 GTCTACGGC 의 시작 = TSS.
창은 TSS +-1200 (HPRC READCOV 폭 2401 과 같게) · 관통은 TSS-1120..TSS+1119.

내는 것
  FHMM_UNITS_HEPG2.tsv    sample slot unit_key read_id cls rel_start rel_end size
  FHMM_READCOV_HEPG2.tsv  sample slot unit_key read_id rel_lo rel_hi
  UNIT_HEPG2.tsv          sample slot unit_key read_id n_cpg n_meth beta n_at n_m6a
  SITE_M6A_HEPG2.tsv      sample unit_key read_id rel 자리 목록 (쉼표) · 비호출기 방법용
"""
import io
import os
import sys
import collections
import pysam

NL = chr(10)
TB = chr(9)
_SELF = os.path.basename(os.path.abspath(__file__))


def _need(k):
    # No defaults are set: if a variable is missing the script stops with a readable message
    v = os.environ.get(k, "")
    if not v:
        raise SystemExit(
            "environment variable %s is not set; it is required by %s." % (k, _SELF)
            + chr(10)
            + "The deposited copy has no default; see README.md section 2")
    return v


W = _need("RDNA_CELLS_DIR")
LABEL = sys.argv[1]
BAMP = sys.argv[2]
HIA = sys.argv[3] if len(sys.argv) > 3 else LABEL
ANCH = "GTCTACGGC"
PAD = 1200
SPAN = 1120
MINGAP = 500          # 이보다 가까운 닻 둘은 같은 유닛 -> 뒤쪽을 쓴다(rfind 와 같음)
MLCUT = 128
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def rc(s):
    return "".join(COMP.get(c, "N") for c in reversed(s))


def parse_mm(seq, mm, ml):
    """MM/ML -> {('A','a'): [(pos, prob)], ...}  pos 는 seq(원본 방향) 좌표."""
    out = {}
    if not mm:
        return out
    blocks = [b for b in mm.rstrip(";").split(";") if b]
    k = 0
    for b in blocks:
        head, _, body = b.partition(",")
        base = head[0]
        modc = head[2:].rstrip(".?")
        deltas = [int(v) for v in body.split(",") if v != ""] if body else []
        idx = [i for i, c in enumerate(seq) if c == base]
        res = []
        j = 0
        for d in deltas:
            j += d
            if j >= len(idx):
                break
            p = idx[j]
            prob = ml[k] if k < len(ml) else 0
            res.append((p, prob))
            k += 1
            j += 1
        out[(base, modc)] = res
    return out


def parse_ma(ma):
    """MA:Z:<len>;nuc.QQQ:s-l,...;msp.:s-l,...;tf.QQQ:s-l,... -> {cls: [(s, e, size)]}"""
    D = {}
    parts = ma.split(";")
    rlen = int(parts[0]) if parts and parts[0].isdigit() else None
    for p in parts[1:]:
        if ":" not in p:
            continue
        nm, _, body = p.partition(":")
        cls = nm.split(".")[0]
        iv = []
        for t in body.split(","):
            if "-" not in t:
                continue
            a, _, b = t.partition("-")
            try:
                s, L = int(a), int(b)
            except ValueError:
                continue
            if L <= 0:
                continue
            iv.append((s, s + L - 1, L))
        D[cls] = iv
    return rlen, D


def anchors(seq):
    """닻 시작 자리 목록 (가까운 둘은 뒤쪽만)."""
    hits = []
    i = seq.find(ANCH)
    while i >= 0:
        hits.append(i)
        i = seq.find(ANCH, i + 1)
    keep = []
    for n, h in enumerate(hits):
        if n + 1 < len(hits) and hits[n + 1] - h < MINGAP:
            continue
        keep.append(h)
    return keep


def flip_iv(iv, L):
    return [(L - 1 - e, L - 1 - s, sz) for (s, e, sz) in iv]


def main():
    fu = io.open(os.path.join(W, "FHMM_UNITS_" + LABEL + ".tsv"), "w", encoding="utf8", newline=NL)
    fu.write(TB.join(["sample", "slot", "unit_key", "read_id", "cls",
                      "rel_start", "rel_end", "size"]) + NL)
    fc = io.open(os.path.join(W, "FHMM_READCOV_" + LABEL + ".tsv"), "w", encoding="utf8", newline=NL)
    fc.write(TB.join(["sample", "slot", "unit_key", "read_id", "rel_lo", "rel_hi"]) + NL)
    fb = io.open(os.path.join(W, "UNIT_" + LABEL + ".tsv"), "w", encoding="utf8", newline=NL)
    fb.write(TB.join(["sample", "slot", "unit_key", "read_id", "n_cpg", "n_meth",
                      "beta", "n_at", "n_m6a", "unit_len"]) + NL)
    fs = io.open(os.path.join(W, "SITE_M6A_" + LABEL + ".tsv"), "w", encoding="utf8", newline=NL)
    fs.write(TB.join(["sample", "unit_key", "read_id", "at_rel", "m6a_rel"]) + NL)

    STAT = collections.Counter()
    SEEN = set()
    MSPM = [0, 0]      # msp 안 (m6A, A/T)
    NUCM = [0, 0]
    ORI = collections.Counter()
    for smp, hia5 in ((LABEL, HIA),):
        bp = BAMP
        if not os.path.exists(bp):
            print("없음", bp)
            continue
        bam = pysam.AlignmentFile(bp, "rb", check_sq=False)
        for r in bam.fetch(until_eof=True):
            STAT[smp + "_read"] += 1
            if r.is_secondary or r.is_supplementary:
                continue
            if (smp, r.query_name) in SEEN:
                STAT[smp + "_dup_read"] += 1
                continue
            SEEN.add((smp, r.query_name))
            seq = r.query_sequence
            if not seq:
                continue
            if r.is_reverse:
                seq = rc(seq)          # 분자(원본) 방향으로 되돌린다
            L = len(seq)
            try:
                mm = r.get_tag("MM")
            except KeyError:
                continue
            try:
                ml = list(r.get_tag("ML"))
            except KeyError:
                ml = []
            try:
                ma = r.get_tag("MA")
            except KeyError:
                STAT[smp + "_noMA"] += 1
                continue
            _, SEG = parse_ma(ma)
            MODS = parse_mm(seq, mm, ml)

            # 방향 — 닻이 정방향이냐
            nf = seq.count(ANCH)
            nr = rc(seq).count(ANCH)
            if nr > nf:
                ORI["rev"] += 1
                seq = rc(seq)
                SEG = dict((c, flip_iv(v, L)) for c, v in SEG.items())
                # C+m 은 CpG 의 C 자리라 뒤집으면 짝 C 로 -1 (L-2-p). A/T 채널은 L-1-p.
                MODS = dict((k, [((L - 2 - p) if k == ("C", "m") else (L - 1 - p), q)
                                 for (p, q) in v])
                            for k, v in MODS.items())
            else:
                ORI["fwd"] += 1
            if max(nf, nr) < 2:
                STAT[smp + "_anchor<2"] += 1
                continue

            # m6A 자리(하드콜) · A/T 분모
            m6a = set()
            for key in (("A", "a"), ("T", "a")):
                for (p, q) in MODS.get(key, []):
                    if q >= MLCUT:
                        m6a.add(p)
            atset = [i for i, c in enumerate(seq) if c == "A" or c == "T"]
            # 5mC 하드콜 (CpG 자리만)
            cmeth = {}
            for (p, q) in MODS.get(("C", "m"), []):
                cmeth[p] = q
            cpg = set(i for i in range(L - 1) if seq[i] == "C" and seq[i + 1] == "G")

            AN = anchors(seq)
            for k, tss in enumerate(AN):
                lo, hi = tss - PAD, tss + PAD
                if tss - SPAN < 0 or tss + SPAN - 1 > L - 1:
                    continue
                STAT["unit_span"] += 1
                ukey = "%s|%s|%d" % (smp, r.query_name, k)
                fc.write(TB.join([smp, hia5, ukey, r.query_name,
                                  "%d" % (max(0, lo) - tss),
                                  "%d" % (min(L - 1, hi) - tss)]) + NL)
                # 발자국
                for cls in ("tf", "nuc", "msp"):
                    for (s, e, sz) in SEG.get(cls, []):
                        c = (s + e) // 2
                        if not (tss - SPAN <= c <= tss + SPAN - 1):
                            continue
                        fu.write(TB.join([smp, hia5, ukey, r.query_name, cls,
                                          "%d" % (s - tss), "%d" % (e - tss),
                                          "%d" % sz]) + NL)
                        STAT[cls] += 1
                # β — 유닛 = 앞 닻+119 .. 이 닻+118 (유전자가 유닛 끝)
                us = (AN[k - 1] + 119) if k > 0 else None
                ue = tss + 118
                if us is None or ue > L - 1:
                    nc = nm_ = 0
                    beta = ""
                    ulen = 0
                else:
                    ulen = ue - us + 1
                    nc = nm_ = 0
                    for i in range(us, ue + 1):
                        if i in cpg and i in cmeth:
                            nc += 1
                            if cmeth[i] >= MLCUT:
                                nm_ += 1
                    beta = ("%.4f" % (nm_ / nc)) if nc else ""
                nat = sum(1 for i in atset if tss - SPAN <= i <= tss + SPAN - 1)
                nm6 = sum(1 for i in m6a if tss - SPAN <= i <= tss + SPAN - 1)
                fb.write(TB.join([smp, hia5, ukey, r.query_name, "%d" % nc,
                                  "%d" % nm_, beta, "%d" % nat, "%d" % nm6,
                                  "%d" % ulen]) + NL)
                # 자리별 (비호출기 방법용) — 좁은 창만 적는다 (-250..+200)
                a0, a1 = tss - 250, tss + 200
                at_r = [i - tss for i in atset if a0 <= i <= a1]
                m6_r = sorted(i - tss for i in m6a if a0 <= i <= a1)
                fs.write(TB.join([smp, ukey, r.query_name,
                                  ",".join("%d" % v for v in at_r),
                                  ",".join("%d" % v for v in m6_r)]) + NL)
            # 검산 — msp 안 / nuc 안 m6A 비율
            for cls, acc in (("msp", MSPM), ("nuc", NUCM)):
                for (s, e, sz) in SEG.get(cls, []):
                    for i in range(max(0, s), min(L - 1, e) + 1):
                        if seq[i] == "A" or seq[i] == "T":
                            acc[1] += 1
                            if i in m6a:
                                acc[0] += 1
        bam.close()
    for f in (fu, fc, fb, fs):
        f.close()
    print("방향 :", dict(ORI))
    for k in sorted(STAT):
        print("  %-18s %d" % (k, STAT[k]))
    print()
    print("검산 — m6A/(A+T)  MSP 안 %.4f (%d/%d) · 뉴클 안 %.4f (%d/%d) · 배수 %.1f"
          % (MSPM[0] / max(MSPM[1], 1), MSPM[0], MSPM[1],
             NUCM[0] / max(NUCM[1], 1), NUCM[0], NUCM[1],
             (MSPM[0] / max(MSPM[1], 1)) / max(NUCM[0] / max(NUCM[1], 1), 1e-9)))


main()
