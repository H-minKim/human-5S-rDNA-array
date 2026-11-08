# -*- coding: utf-8 -*-
# 그림 3 : 종마다 스페이서 종류 전부 + 각 종류의 유닛 내부 구조 + 문헌 대조 결과
#   입력 : lit/SPACER_TYPES.tsv · lit/UNIT_VERDICT.tsv
#   자 : 0~2,400 bp 선형. 그보다 긴 유닛은 끝을 톱니로 끊고 실제 길이를 적는다
# MIT License · Copyright (c) 2026 Hungmin Kim
import io
import os
import collections

NL = chr(10)
TAB = chr(9)
B = os.path.dirname(os.path.abspath(__file__))
L = os.path.join(B, "lit")
SPN = {"human_CHM13v2": "사람 (CHM13v2)", "chimp_mPanTro3": "침팬지", "bonobo_mPanPan1": "보노보",
       "gorilla_mGorGor1": "고릴라", "orang_sumatran_mPonAbe1": "오랑우탄 (수마트라)",
       "orang_bornean_mPonPyg2": "오랑우탄 (보르네오)", "siamang_mSymSyn1": "시아망",
       "rhesus_T2T_MMU8v2": "붉은털원숭이", "mouse_C57BL6J_T2T": "마우스",
       "zfinch_v03_dip": "금화조", "zebrafish_GRCz12tu": "제브라피시", "arabidopsis_ColCC": "애기장대"}
ORDER = ["human_CHM13v2", "chimp_mPanTro3", "bonobo_mPanPan1", "gorilla_mGorGor1",
         "orang_sumatran_mPonAbe1", "orang_bornean_mPonPyg2", "siamang_mSymSyn1",
         "rhesus_T2T_MMU8v2", "mouse_C57BL6J_T2T", "zfinch_v03_dip", "zebrafish_GRCz12tu",
         "arabidopsis_ColCC"]


def rdtsv(p):
    f = io.open(p, encoding="utf8")
    hd = f.readline().rstrip(NL).split(TAB)
    return hd, [dict(zip(hd, l.rstrip(NL).split(TAB))) for l in f if l.strip()]


_, ST = rdtsv(os.path.join(L, "SPACER_TYPES.tsv"))
_, VD = rdtsv(os.path.join(L, "UNIT_VERDICT.tsv"))
V = dict((r["array_id"], r) for r in VD)
# 종 > 종류 묶기
TY = collections.defaultdict(lambda: collections.defaultdict(list))
for r in ST:
    TY[r["species"]][int(r["spacer_type"])].append(r)

MAXBP = 2400.0
BW = 520.0
X0 = 270.0
ROWH = 26
GAPH = 14
HEADH = 24
GENE = "#B85C1E"      # 유전자 119 bp
SPC = "#2C6E9B"       # 스페이서
SPC2 = "#7FA8C4"      # 스페이서 (문헌 검증 안 된 종류)
INK = "#1F2328"
MUT = "#6B7280"
GRID = "#E3E6EA"

rows = []   # (kind, ...)
for sp in ORDER:
    if sp not in TY:
        continue
    rows.append(("head", sp))
    for t in sorted(TY[sp]):
        rows.append(("type", sp, t, TY[sp][t]))
H = 118 + sum((HEADH + 6) if r[0] == "head" else ROWH for r in rows) + 40
W = 1064
o = []
o.append('<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d" font-family="Malgun Gothic, Segoe UI, sans-serif">' % (W, H, W, H))
o.append('<rect width="%d" height="%d" fill="#FFFFFF"/>' % (W, H))
o.append('<text x="24" y="34" font-size="17" font-weight="700" fill="%s">5S rDNA — 종마다 스페이서 종류 전부와 유닛 내부 구조</text>' % INK)
o.append('<text x="24" y="54" font-size="11.5" fill="%s">종류 = 유전자 119 bp 를 떼고 스페이서만 비교(blastn 일치 80%% 이상 · 긴 쪽의 50%% 이상 정렬 · 단일연결). 막대는 유닛 하나 = [스페이서][유전자].</text>' % MUT)
o.append('<text x="24" y="70" font-size="11.5" fill="%s">진한 파랑 = 공개 반복단위와 전장 일치한 종류 · 연한 파랑 = 대조할 같은 종 공개 단위가 없거나 부분만 상동. 2,400 bp 를 넘는 유닛은 끝을 톱니로 끊었다.</text>' % MUT)
# 자
ty = 100
for bp in (0, 500, 1000, 1500, 2000, 2400):
    x = X0 + BW * bp / MAXBP
    o.append('<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" stroke="%s" stroke-width="1"/>' % (x, ty, x, H - 30, GRID))
    o.append('<text x="%.1f" y="%d" font-size="10" fill="%s" text-anchor="middle">%s</text>' % (x, ty - 6, MUT, ("%d bp" % bp) if bp else "0"))
y = ty + 10
for r in rows:
    if r[0] == "head":
        sp = r[1]
        nA = len(ST and [x for x in ST if x["species"] == sp])
        nT = len(TY[sp])
        cps = sum(int(x["copies"]) for x in ST if x["species"] == sp)
        y += HEADH
        o.append('<text x="24" y="%d" font-size="13" font-weight="700" fill="%s">%s</text>' % (y, INK, SPN.get(sp, sp)))
        o.append('<text x="%d" y="%d" font-size="11" fill="%s">array %d · copy %s · 스페이서 종류 %d</text>'
                 % (206, y, MUT, nA, "{:,}".format(cps), nT))
        y += 6
        continue
    sp, t, mem = r[1], r[2], r[3]
    y += ROWH
    ul = int(mem[0]["unit_len"])
    slen = int(mem[0]["spacer_len"])
    lens = sorted(set(int(m["unit_len"]) for m in mem))
    cps = sum(int(m["copies"]) for m in mem)
    lit = [V[m["array_id"]]["lit_verdict"] for m in mem if m["array_id"] in V]
    ok = any(x == "공개 단위와 전장 일치" for x in lit)
    col = SPC if ok else SPC2
    shown = min(ul, MAXBP)
    bw = BW * shown / MAXBP
    gw = BW * min(119, shown) / MAXBP
    o.append('<rect x="%.1f" y="%.1f" width="%.1f" height="13" rx="3" fill="%s"/>' % (X0, y - 10, max(bw - gw, 1), col))
    o.append('<rect x="%.1f" y="%.1f" width="%.1f" height="13" rx="3" fill="%s"/>' % (X0 + bw - gw, y - 10, gw, GENE))
    if ul > MAXBP:
        xe = X0 + BW
        pts = " ".join("%.1f,%.1f" % (xe - 5 + (4 if i % 2 else 0), y - 10 + i * 3.25) for i in range(5))
        o.append('<polyline points="%s" fill="none" stroke="#FFFFFF" stroke-width="2"/>' % pts)
    lab = "종류 %d" % t
    o.append('<text x="%d" y="%.1f" font-size="11" font-weight="600" fill="%s">%s</text>' % (36, y, INK, lab))
    o.append('<text x="%d" y="%.1f" font-size="10.5" fill="%s">array %d · copy %s</text>' % (92, y, MUT, len(mem), "{:,}".format(cps)))
    txt = "%s bp = 스페이서 %s + 유전자 119" % (("{:,}".format(ul)) if len(lens) == 1 else "%d-%d" % (lens[0], lens[-1]), "{:,}".format(slen))
    o.append('<text x="%.1f" y="%.1f" font-size="10.5" fill="%s">%s</text>' % (X0 + BW + 10, y, INK if ok else MUT, txt))
    y += 0
o.append('<text x="24" y="%d" font-size="10.5" fill="%s">붉은털원숭이 종류 4·5 와 제브라피시 긴 종류들은 유닛 하나가 아니라 유전자 사이가 실제로 수 kb 인 것이다 — 06_aux/lit_reference/README.md 참조.</text>' % (H - 12, MUT))
o.append('</svg>')
io.open(os.path.join(B, "fig3.svg"), "w", encoding="utf8", newline=NL).write(NL.join(o))
print("fig3.svg · 행 %d" % len(rows))
