# -*- coding: utf-8 -*-
# Strip the statistics out of the legends, tick labels and annotations.
# What is removed:
# a) statistical fragments: n=, N=, mean 0.65, median 0.047, p=, SD, the plus-minus sign
# b) a bare number in parentheses: (298 ind.), (n=1,094), (N=44 units)
# c) a panel number: a leading (2) or a circled digit
# What is kept: the series names, the group names, "top/bottom n%", and method arguments such as cut 0.3 or ML/255.
import re

# only a statistic carrying a number; the words mean and median are themselves quantity names and are kept
STATNUM = re.compile(r"(\bn\s*=\s*[\d,]+"
                     r"|\bN\s*=\s*[\d,]+"
                     r"|\b(?:mean|median|중앙값|평균|중앙)\s*[=:]?\s*\d[\d.,]*"
                     r"|\bp\s*[=<]\s*\d[\d.e-]*"
                     r"|\bSD\s*[=:]?\s*\d[\d.]*"
                     r"|±\s*\d[\d.]*)", re.I)
# a bare number in parentheses, carrying only a unit
BARECNT = re.compile(r"\s*[\(\[]\s*[Nn]?\s*=?\s*[\d,]+\s*"
                     r"(?:ind\.?|individuals|units|slots|haplotypes|molecules|read-units|reads)?"
                     r"\s*[\)\]]")
# a leading panel number
PANELNO = re.compile(r"^\s*(?:\(\d+\)|[①-⑩])\s*")
# parentheses holding a statistic or a summary: (median 6 units), (median, bars are the IQR), (median over 4 haplotypes)
# What is kept: (bottom 21.9%), (hard call, ML >= 128), (120 bins per unit), (bp), (log)
PARENSTAT = re.compile(r"\s*[\(\[][^()\[\]]*?"
                       r"(?:\bmean\b|\bmedian\b|중앙값|중앙|평균|\bn\s*=|\bN\s*=|\bIQR\b|\bSD\b|±)"
                       r"[^()\[\]]*?[\)\]]", re.I)


def make(fix_sym, rec):
    """build legend_label from the hook's fix_sym and _rec"""

    def legend_label(raw):
        if raw is None:
            return raw
        src = str(raw)
        t = PANELNO.sub("", src)
        t = BARECNT.sub("", t)
        # a newline is a separator as well: there are tick labels of the form "stratum/nn=21" (fig42 J)
        nl = "\n" if "\n" in t else "  "
        parts = re.split(r"\n|\s{2,}", t)
        keep, drop = [], []
        for p in parts:
            if not p.strip():
                continue
            (drop if STATNUM.search(p) else keep).append(p)
        if not keep:                      # anything that would be left empty once everything is stripped is left as it is
            keep, drop = [parts[0]], parts[1:]
        out = fix_sym(nl.join(keep))
        if drop or t != src:
            rec("legend_stat", "-", src)
        return out

    return legend_label


def strip_stat(s):
    """remove only the statistics part from annotations (text, annotate); value labels stay
       "12 / n=42" -> "12"   ,   "median 0.0473" -> ""   ,   "0.37x" -> unchanged"""
    t = str(s)
    if not STATNUM.search(t):
        return t, False
    out = BARECNT.sub("", t)
    out = STATNUM.sub("", out)
    out = re.sub(r"[ \t]{2,}", " ", out)
    out = re.sub(r"\n{2,}", "\n", out)
    out = out.strip(" \t\n,;:()")
    return out, True
