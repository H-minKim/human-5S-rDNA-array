# -*- coding: utf-8 -*-
r"""Definition of the terminal-type unit, set in this one place only
   12 diagnostic variants = the 12 of Supplementary Figure S2b I (carriage rate difference >= 0.70 between units carrying and not carrying D1096:1, table Y5)
   Terminal type = a unit carrying TMIN (10) or more of these 12
   Earlier definition = D1096:1 alone; with the environment variable TERM_DEF=D1096 that definition runs
   Tokens are v4 tokens under the majority-allele hook (hook_tier). None of the 12 is the majority allele at a majority-basis (flipped) position, so the names are the same with or without the hook."""
import os
DIAG12 = ("D1096:1", "635CT", "1164GA", "740GA", "DB43:2", "DB16:17", "488GC", "DB46:10", "DB58:1", "R1752CG", "R405CG", "R347CG")
TMIN = 10
MODE = os.environ.get("TERM_DEF", "k10")
assert MODE in ("k10", "D1096"), "TERM_DEF 는 k10 또는 D1096"


def n_diag(toks):
    """Number of the 12 diagnostic variants present in a unit token collection (set or list)"""
    return sum(1 for t in DIAG12 if t in toks)


def is_term(toks):
    if MODE == "D1096":
        return "D1096:1" in toks
    return n_diag(toks) >= TMIN
