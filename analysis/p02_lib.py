# -*- coding: utf-8 -*-
"""p02_pilot_bam.py 에서 검증된(표와 100% 동일) MM/ML·CIGAR 함수를 모듈로 뺀 것 — 내용 동일."""
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def parse_mm(mm, ml, oseq):
    """원본 read 좌표 -> {(base, strand, code): ({op: ML}, idx)}"""
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
    """SEQ 좌표(0기반) -> 정본 좌표(1기반), M/=/X 만"""
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
