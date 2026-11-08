# -*- coding: utf-8 -*-
r"""Variant classification without filtering: every variant in the data falls into exactly one group

 12 substitution types, each kept separate (C>T and G>A further split by CpG context) = 14 groups
 deletion/insertion x inside/outside repeat block x length (1bp, 2-9bp, >=10bp)  =  8 groups
 a form that fits nowhere                                                = 'other'
 => 23 groups in total. In the full survey (VTINV) 'other' was 0.
"""
import re

SUBS = ["C>T", "G>A", "C>A", "G>T", "C>G", "G>C", "T>C", "A>G", "T>G", "A>C", "T>A", "A>T"]
_SNV = re.compile(r"^(R?)(\d+)([ACGT])([ACGT])$")
_D = re.compile(r"^D(\d+):(\d+)$")
_I = re.compile(r"^I(\d+):(\d+):")
_DB = re.compile(r"^DB(\d+):(\d+)$")
_IB = re.compile(r"^IB(\d+):(\d+)$")


def _len3(n):
    return "1bp" if n == 1 else ("2-9bp" if n < 10 else ">=10bp")


# drawing order: deamination -> oxidation -> other substitutions -> indel -> other
GROUPS = (["C>T|CpG", "C>T|비CpG", "G>A|CpG", "G>A|비CpG"]
          + ["C>A", "G>T"]
          + ["C>G", "G>C", "T>C", "A>G", "T>G", "A>C", "T>A", "A>T"]
          + ["결실 1bp", "결실 2-9bp", "결실 >=10bp",
             "삽입 1bp", "삽입 2-9bp", "삽입 >=10bp",
             "블록 결실", "블록 삽입", "기타"])
GI = {g: i for i, g in enumerate(GROUPS)}

# colors: related groups share a color family
COLOR = {
    "C>T|CpG": "#7b241c", "C>T|비CpG": "#c0392b", "G>A|CpG": "#cd6155", "G>A|비CpG": "#e6b0aa",
    "C>A": "#b9770e", "G>T": "#f5b041",
    "C>G": "#1a5276", "G>C": "#2874a6", "T>C": "#5499c7", "A>G": "#a9cce3",
    "T>G": "#117a65", "A>C": "#45b39d", "T>A": "#7d6608", "A>T": "#d4ac0d",
    "결실 1bp": "#512e5f", "결실 2-9bp": "#76448a", "결실 >=10bp": "#a569bd",
    "삽입 1bp": "#145a32", "삽입 2-9bp": "#1e8449", "삽입 >=10bp": "#58d68d",
    "블록 결실": "#5d6d7e", "블록 삽입": "#aeb6bf", "기타": "#000000"}


def iscpg_factory(CAN):
    def iscpg(p):
        c = CAN[p - 1]
        return ((c == "C" and p < 2240 and CAN[p] == "G")
                or (c == "G" and p >= 2 and CAN[p - 2] == "C"))
    return iscpg


def vgroup(x, iscpg):
    """one variant token -> one group (always exactly one)"""
    m = _DB.match(x)
    if m:
        return "블록 결실"
    m = _IB.match(x)
    if m:
        return "블록 삽입"
    m = _SNV.match(x)
    if m:
        p, r, a = int(m.group(2)), m.group(3), m.group(4)
        s = r + ">" + a
        if s in ("C>T", "G>A"):
            return s + ("|CpG" if iscpg(p) else "|비CpG")
        return s if s in GI else "기타"
    m = _D.match(x)
    if m:
        return "결실 " + _len3(int(m.group(2)))
    m = _I.match(x)
    if m:
        return "삽입 " + _len3(int(m.group(2)))
    return "기타"
