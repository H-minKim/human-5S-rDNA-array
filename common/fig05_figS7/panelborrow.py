# -*- coding: utf-8 -*-
"""Borrow a panel block as it is from the original figure code.

Why : the figure is assembled from existing panels.
     A copy (1) introduces transcription errors and
     (2) does not follow when the original is corrected.
     So the block is lifted from the original file and executed in place.

Block boundary : from `ax = SUB("<letter>")` to the line starting with `tag(ax, "<letter>"` (inclusive).
              If `tag(` spans several lines, lines are taken until the parenthesis closes.
"""
import io
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
_CACHE = {}


def _src(fn):
    if fn not in _CACHE:
        _CACHE[fn] = io.open(os.path.join(HERE, fn), encoding="utf8").read().split(chr(10))
    return _CACHE[fn]


def grab(fn, letter):
    """return the source text of the block of panel <letter> in the original file fn."""
    L = _src(fn)
    start = None
    pat_s = 'ax = SUB("%s")' % letter
    pat_t = 'tag(ax, "%s"' % letter
    for i, ln in enumerate(L):
        if ln.strip().startswith(pat_s):
            start = i
            break
    if start is None:
        raise KeyError("%s: panel %s start (%s) not found" % (fn, letter, pat_s))
    end = None
    for j in range(start, len(L)):
        if L[j].lstrip().startswith(pat_t):
            # a tag(...) call can span several lines, up to the line where the bracket closes
            depth = 0
            for k in range(j, len(L)):
                depth += L[k].count("(") - L[k].count(")")
                if depth <= 0:
                    end = k
                    break
            break
    if end is None:
        raise KeyError("%s: panel %s end (%s) not found" % (fn, letter, pat_t))
    import textwrap
    # in the original it may sit inside an if block (e3 of F1), so the indentation is removed
    return textwrap.dedent(chr(10).join(L[start:end + 1]))


def borrow(fn, letter, ax, ns, subst=None):
    """Draw the block on ax. ns = the namespace to execute in (usually globals()).

    The `ax = SUB("x")` line is changed to use the ax passed in.
    subst = [(old text, new text), ...] can alter text inside the block further (use with care).
    """
    code = grab(fn, letter)
    code = re.sub(r'^ax = SUB\("%s"\)' % re.escape(letter),
                  "ax = _BORROW_AX", code, count=1, flags=re.M)
    for a, b in (subst or []):
        if a not in code:
            raise KeyError("%s(%s): no text to replace: %r" % (fn, letter, a))
        code = code.replace(a, b)
    g = dict(ns)
    g["_BORROW_AX"] = ax
    exec(compile(code, "<borrow %s:%s>" % (fn, letter), "exec"), g)
    return code


def lines(fn, letter):
    return grab(fn, letter).count(chr(10)) + 1
