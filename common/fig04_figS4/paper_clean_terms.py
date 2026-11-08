# -*- coding: utf-8 -*-
# Dictionary of figure terms and labels for the manuscript
# 사전은 문장(라벨) 단위 - 조사,어미를 건드리지 않는다.
# 긴 키부터 적용한다.
# 코드,열 이름,저장된 표는 바꾸지 않는다. 바꾸는 것은 그림에 찍히는 글자뿐.
#
# RAW_TERMS : 기호를 바꾸기 전(원본 리터럴 그대로)에 맞춘다 - 설명이 붙은 라벨을 명사구로
# TERMS : 기호를 바꾼 뒤에 맞춘다 - 용어 교체
import re

# ── RAW (원본 리터럴) : 설명,판정이 붙은 축 라벨을 명사구로 ───────────────
RAW_TERMS = [
    # 정의식은 legend 로 (fig44)
    ("purity = min(frac of molecules <=0.35, frac >=0.7)", "purity"),
    # 문장형 축 라벨 -> 명사구 (figA6 C5)
    ("how far the FIRST unit still is from the new state", "distance from the new state"),
    # 판정 기준은 legend 로 (figF6) - 한글판은 아래 RAW_RE 가 받는다
    ("unit layer: Ashman D  (>2 = two peaks)", "unit layer: Ashman D"),
    # 점,막대의 뜻은 legend 로 (figA6 C5)
    ("low-site count (median; bar = IQR)", "low-site count"),
]
# 한글 '유닛 층' 은 글자 조합이 환경마다 달라 보일 수 있어 정규식으로도 받는다
RAW_RE = [
    (re.compile(r"^(유닛\s*층\s*[:：]\s*Ashman D)\s*\(.*\)$"), r"\1"),
    (re.compile(r"^(꺼진\s*자리\s*수)\s*\(.*\)$"), r"\1"),
]

# ── 용어 (기호 치환 뒤) ────────────────────────────────────────────────
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
    """기호 치환 전에 — 원본 리터럴 사전"""
    if not t:
        return t
    for a, b in RAW_TERMS:
        if a in t:
            t = t.replace(a, b)
    for rx, rep in RAW_RE:
        t = rx.sub(rep, t)
    return t


def apply_terms(t):
    """기호 치환 뒤에 — 용어 사전"""
    if not t:
        return t
    for a, b in TERMS:
        if a in t:
            t = t.replace(a, b)
    return t
