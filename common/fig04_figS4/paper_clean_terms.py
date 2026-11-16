# -*- coding: utf-8 -*-
# Dictionary of figure terms and labels for the manuscript
# The dictionary works on whole labels, so that no inflection is touched.
# The longer keys are applied first.
# The code, the column names and the stored tables are not changed; only the text printed on the figure is.
#
# RAW_TERMS is matched before the symbols are substituted, against the original literal, and turns an explanatory label into a noun phrase.
# TERMS is matched after the symbols are substituted, and replaces terminology.
import re

# -- RAW, the original literals: an axis label carrying an explanation or a verdict becomes a noun phrase ───────────────
RAW_TERMS = [
    # a defining formula belongs in the legend (fig44)
    ("purity = min(frac of molecules <=0.35, frac >=0.7)", "purity"),
    # a sentence-form axis label becomes a noun phrase (figA6 C5)
    ("how far the FIRST unit still is from the new state", "distance from the new state"),
    # the decision criterion belongs in the legend (figF6)
    ("unit layer: Ashman D  (>2 = two peaks)", "unit layer: Ashman D"),
    # what a point or a bar means belongs in the legend (figA6 C5)
    ("low-site count (median; bar = IQR)", "low-site count"),
]
# some of these labels can render differently in different environments, so they are also matched by regular expression
RAW_RE = [
    (re.compile(r"^(유닛\s*층\s*[:：]\s*Ashman D)\s*\(.*\)$"), r"\1"),
    (re.compile(r"^(꺼진\s*자리\s*수)\s*\(.*\)$"), r"\1"),
]

# -- the terminology, after the symbols are substituted ────────────────────────────────────────────────
TERMS = [
    ("쓸 수 있었던 슬롯 수 (로그)", "쓸 수 있었던 하플로타입 수 (로그)"),
    ("slots available (log)", "haplotypes available (log)"),
    ("슬롯 수", "하플로타입 수"),
    ("slots", "haplotypes"),
    ("slots", "haplotype"),
    ("slot", "haplotype"),
    ("HPRC v11", "HPRC"),
    ("1kGP v11", "1kGP"),
    ("v11 0.055", "0.055"),
    ("정본", "참조서열"),
    ("canonical frame", "reference frame"),
    ("NTS1", "spacer-1"),
    ("NTS2", "spacer-2"),
]
RAW_TERMS = sorted(RAW_TERMS, key=lambda kv: -len(kv[0]))
TERMS = sorted(TERMS, key=lambda kv: -len(kv[0]))


def apply_raw(t):
    """before symbol substitution: dictionary of original literals"""
    if not t:
        return t
    for a, b in RAW_TERMS:
        if a in t:
            t = t.replace(a, b)
    for rx, rep in RAW_RE:
        t = rx.sub(rep, t)
    return t


def apply_terms(t):
    """after symbol substitution: dictionary of terms"""
    if not t:
        return t
    for a, b in TERMS:
        if a in t:
            t = t.replace(a, b)
    return t
