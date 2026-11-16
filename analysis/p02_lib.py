# -*- coding: utf-8 -*-
"""The MM/ML and CIGAR functions of p02_pilot_bam.py (verified 100% identical to the table) taken out as a module; same content."""
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def parse_mm(mm, ml, oseq):
    """original read coordinates -> {(base, strand, code): ({op: ML}, idx)}"""
    out = {}
    k = 0
    for blk in [b for b in mm.split(";") if b]:
        parts = blk.split(",")
        head = parts[0]
        base = head[0]
        strand = head[1]
        code = head[2:].rstrip(".?")
        skips = [int(s) for s in parts[1:]]
        idx = [i for i, c in enumerate(oseq) if c == base] if base != "N" else list(range(len(oseq)))
        d = {}
        j = -1
        for s in skips:
            j += s + 1
            if j < len(idx):
                d[idx[j]] = ml[k]
            k += 1
        out[(base, strand, code)] = (d, idx)
    return out


def seq2ref(r):
    """SEQ coordinate (0-based) -> canonical coordinate (1-based), M/=/X only"""
    s2r = {}
    rp, sp = r["pos"], 0
    for op, n in r["cigar"]:
        if op in "M=X":
            for q in range(n):
                s2r[sp + q] = rp + q + 1
            rp += n; sp += n
        elif op in "IS":
            sp += n
        elif op in "DN":
            rp += n
    return s2r
