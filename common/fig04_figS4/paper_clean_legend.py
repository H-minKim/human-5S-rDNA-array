# -*- coding: utf-8 -*-
# 범례,눈금,주석에서 통계만 뗀다 (논문용 그림 공통 규칙 §1, 2026-10-07)
# 빼는 것
# ⒜ 통계 토막 : n=, N=, mean 0.65, median 0.047, 중앙값 0.4, p=, SD, ±
# ⒝ 괄호 속 맨 수치 : (298 ind.), (n=1,094), (N=44 units)
# ⒞ 패널 번호 : 앞머리의 (2), ②
# 남기는 것 : 계열 이름, 무리 이름, "상위/하위 n%", 방법 인자(cut 0.3, ML/255)
import re

# 숫자가 붙은 통계만 (낱말 'mean','median' 자체는 양 이름이라 남긴다)
STATNUM = re.compile(r"(\bn\s*=\s*[\d,]+"
                     r"|\bN\s*=\s*[\d,]+"
                     r"|\b(?:mean|median|중앙값|평균|중앙)\s*[=:]?\s*\d[\d.,]*"
                     r"|\bp\s*[=<]\s*\d[\d.e-]*"
                     r"|\bSD\s*[=:]?\s*\d[\d.]*"
                     r"|±\s*\d[\d.]*)", re.I)
# 괄호 속이 맨 수치(단위만 붙은 것)
BARECNT = re.compile(r"\s*[\(\[]\s*[Nn]?\s*=?\s*[\d,]+\s*"
                     r"(?:ind\.?|individuals|units|slots|haplotypes|molecules|read-units|reads)?"
                     r"\s*[\)\]]")
# 앞머리 패널 번호
PANELNO = re.compile(r"^\s*(?:\(\d+\)|[①-⑩])\s*")
# 괄호 안이 통계,요약인 것 : (median 6 units), (중앙값, 막대는 IQR), (median over 4 haplotypes)
# 남기는 것 : (bottom 21.9%), (hard call, ML >= 128), (120 bins per unit), (bp), (log)
PARENSTAT = re.compile(r"\s*[\(\[][^()\[\]]*?"
                       r"(?:\bmean\b|\bmedian\b|중앙값|중앙|평균|\bn\s*=|\bN\s*=|\bIQR\b|\bSD\b|±)"
                       r"[^()\[\]]*?[\)\]]", re.I)


def make(fix_sym, rec):
    """훅의 fix_sym·_rec 을 받아 legend_label 을 만든다"""

    def legend_label(raw):
        if raw is None:
            return raw
        src = str(raw)
        t = PANELNO.sub("", src)
        t = BARECNT.sub("", t)
        # 줄바꿈도 구분기호다 - "stratum\nn=21" 꼴 눈금 라벨이 있다 (fig42 J)
        nl = "\n" if "\n" in t else "  "
        parts = re.split(r"\n|\s{2,}", t)
        keep, drop = [], []
        for p in parts:
            if not p.strip():
                continue
            (drop if STATNUM.search(p) else keep).append(p)
        if not keep:                      # 다 떼면 빈칸이 되는 것은 그대로 둔다
            keep, drop = [parts[0]], parts[1:]
        out = fix_sym(nl.join(keep))
        if drop or t != src:
            rec("legend_stat", "-", src)
        return out

    return legend_label


def strip_stat(s):
    """주석(text·annotate)에서 통계 토막만 뺀다 — 값 라벨은 남는다
       "12 / n=42" -> "12"   ·   "median 0.0473" -> ""   ·   "0.37x" -> 그대로"""
    t = str(s)
    if not STATNUM.search(t):
        return t, False
    out = BARECNT.sub("", t)
    out = STATNUM.sub("", out)
    out = re.sub(r"[ \t]{2,}", " ", out)
    out = re.sub(r"\n{2,}", "\n", out)
    out = out.strip(" \t\n,;:()")
    return out, True
