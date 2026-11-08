# -*- coding: utf-8 -*-
"""패널 블록을 원본 그림 코드에서 그대로 빌려 온다 (2026-10-09).

왜 : 새 그림(fig5)은 기존 패널을 모아 만드는 것이다.
     베껴 쓰면 ⑴ 전사 오류가 생기고 ⑵ 원본이 고쳐져도 따라오지 않는다
     (기억 `staleness-must-follow-the-chain`).
     ⇒ 원본 파일에서 블록을 떠서 그 자리에서 실행한다. 사슬이 끊기지 않는다.

블록의 경계 : `ax = SUB("<글자>")` 로 시작해 `tag(ax, "<글자>"` 로 시작하는 줄까지(포함).
              `tag(` 가 여러 줄이면 괄호가 닫히는 줄까지 더 먹는다.
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
    """원본 파일 fn 에서 패널 <letter> 블록의 소스 글자를 돌려준다."""
    L = _src(fn)
    start = None
    pat_s = 'ax = SUB("%s")' % letter
    pat_t = 'tag(ax, "%s"' % letter
    for i, ln in enumerate(L):
        if ln.strip().startswith(pat_s):
            start = i
            break
    if start is None:
        raise KeyError("%s : 패널 %s 의 시작(%s)을 못 찾았다" % (fn, letter, pat_s))
    end = None
    for j in range(start, len(L)):
        if L[j].lstrip().startswith(pat_t):
            # tag(...) 가 여러 줄일 수 있다 — 괄호가 닫히는 줄까지
            depth = 0
            for k in range(j, len(L)):
                depth += L[k].count("(") - L[k].count(")")
                if depth <= 0:
                    end = k
                    break
            break
    if end is None:
        raise KeyError("%s : 패널 %s 의 끝(%s)을 못 찾았다" % (fn, letter, pat_t))
    import textwrap
    # 원본에서 `if ... :` 안에 들어가 있을 수 있다(F1 의 e3) ⇒ 들여쓰기를 푼다
    return textwrap.dedent(chr(10).join(L[start:end + 1]))


def borrow(fn, letter, ax, ns, subst=None):
    """블록을 ax 에 그린다. ns = 실행에 쓸 이름 공간(보통 globals()).

    `ax = SUB("x")` 줄은 넘겨받은 ax 를 쓰도록 바꾼다.
    subst = [(옛글자, 새글자), ...] 로 블록 안 글자를 더 고칠 수 있다(조심해서 쓸 것).
    """
    code = grab(fn, letter)
    code = re.sub(r'^ax = SUB\("%s"\)' % re.escape(letter),
                  "ax = _BORROW_AX", code, count=1, flags=re.M)
    for a, b in (subst or []):
        if a not in code:
            raise KeyError("%s(%s) : 바꿀 글자가 없다 — %r" % (fn, letter, a))
        code = code.replace(a, b)
    g = dict(ns)
    g["_BORROW_AX"] = ax
    exec(compile(code, "<borrow %s:%s>" % (fn, letter), "exec"), g)
    return code


def lines(fn, letter):
    return grab(fn, letter).count(chr(10)) + 1
