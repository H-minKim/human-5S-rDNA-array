# -*- coding: utf-8 -*-
r"""n and filter notation for paper figures. Builds a one-line text and places it at the bottom of the figure
Format (paper convention):
   n = <unit> (<individuals>), inclusion: <filter>, cohort: <source composition>
All numbers are measured values (`NCOUNT_by_result.tsv`, 2026-10-01).
"""

# ---- measured denominators (output of code/NCOUNT_hprc.py) ----
COHORT = {
    564: dict(ind=305, hprc=558, giab=6, pub=374, hifi=185, cons2=5, filt="primary"),
    562: dict(ind=304, hprc=556, giab=6, pub=374, hifi=183, cons2=5, filt="유닛 >= 30"),
    549: dict(ind=303, hprc=543, giab=6, pub=366, hifi=178, cons2=5, filt="유닛 >= 40"),
    502: dict(ind=298, hprc=496, giab=6, pub=333, hifi=164, cons2=5, filt="유닛 >= 60"),
    492: dict(ind=295, hprc=486, giab=6, pub=328, hifi=159, cons2=5, filt="유닛 >= 61 (30+MID+30)"),
}
PRESET_KO = "primary (생식세포 · 단계 B/E · 끊김 · HG002_2 제외)"
PRESET_EN = "primary set (germline; stages B/E, break slots and HG002_2 excluded)"

def ko(n, extra="", unit="하플로타입"):
    c = COHORT[n]
    f = c["filt"]
    f = PRESET_KO if f == "primary" else (f + " · " + PRESET_KO)
    s = ("n = %d %s (%d 개체)  ·  포함 기준 : %s  ·  "
         "코호트 : HPRC %d + GIAB %d (HPRC 공개 조립본 %d · 우리 hifiasm %d · read 합의2 %d)"
         % (n, unit, c["ind"], f, c["hprc"], c["giab"], c["pub"], c["hifi"], c["cons2"]))
    return s + ("  ·  " + extra if extra else "")

def en(n, extra="", unit="haplotypes"):
    c = COHORT[n]
    fm = {"primary": "", "유닛 >= 30": "arrays with >=30 units",
         "유닛 >= 40": "arrays with >=40 units", "유닛 >= 60": "arrays with >=60 units",
         "유닛 >= 61 (30+MID+30)": "arrays with >=61 units (30 + middle + 30)"}[c["filt"]]
    f = PRESET_EN if not fm else (fm + "; " + PRESET_EN)
    s = ("n = %d %s (%d individuals)  ·  inclusion: %s  ·  "
         "cohort: HPRC %d + GIAB %d (HPRC public assemblies %d; our hifiasm %d; read consensus %d)"
         % (n, unit, c["ind"], f, c["hprc"], c["giab"], c["pub"], c["hifi"], c["cons2"]))
    return s + ("  ·  " + extra if extra else "")

def put(fig, text, y=0.004, size=8.2, color="#333333"):
    """Places the text as one line at the very bottom of the figure"""
    fig.text(0.5, y, text, ha="center", va="bottom", fontsize=size, color=color)
