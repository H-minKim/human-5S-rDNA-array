# -*- coding: utf-8 -*-
r"""3-way split of variant origin (canonical definition, fork 4815ef09 code/V2_apply.py)

   cohort = slots with >= 30 units (NS), f = number of slots carrying the token / NS, fT = 0.10
   multi-individual shared type  f >= fT                  <- 'shared between haplotypes' = inherited
   within-array multiple type    f <  fT & k >= 2 within that slot
   single-unit type              f <  fT & k == 1 within that slot

Caution (trio measurement): the 'single-unit type' is not de novo.
   28 of 29 in R1 (97%) were in the same unit of the transmitting parent, and read-supported de novo unit variants were 0
   (upper bound <= 2/447 units per generation). So this layer must be read not as "newly arisen variants" but as
   "where rare variants remain". Mutation rates can only come from trio and generation data.
"""
import collections

FT = 0.10
MINU = 30


def build(U, HK):
    """-> (hapf, token unit counts per slot), hapf[t] = f"""
    S = {k: [set(r[4]) for r in U[k]] for k in HK}
    KS = [k for k in HK if len(S[k]) >= MINU]
    NS = len(KS)
    hapn = collections.Counter()
    for k in KS:
        for t in set().union(*S[k]) if S[k] else set():
            hapn[t] += 1
    hapf = {t: c / NS for t, c in hapn.items()}
    kcnt = {}
    for k in HK:
        c = collections.Counter()
        for s in S[k]:
            c.update(s)
        kcnt[k] = c
    return hapf, kcnt, NS, len(KS)


def klass(t, k, hapf, kcnt):
    f = hapf.get(t, 0.0)
    if f >= FT:
        return "다개체 공유형"
    return "배열내 다중형" if kcnt[k].get(t, 0) >= 2 else "단일유닛형"
