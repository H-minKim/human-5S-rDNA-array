# -*- coding: utf-8 -*-
# Fig 1 A — 유전체 방향 유닛(색만) → X 교차 → 우리가 쓰는 유닛(이름·눈금)
#            + 세 곳 모두 전사 방향(TSS 꺾인 화살표) + chr1 왼쪽 = 1
#
# (마지막 판 · 손그림)
#     *"우선은 뒤집힌 상태의 그림을 그리되(영역 표시 없이 색깔로만 영역 표시),
#       그 아래에 우리가 쓰는 유닛형태 두기. 모든 경우에 transcription 방향 표시하기"*
#     ⇒ 층 셋 : 배열 → (확대) 유전체 방향 유닛 → X → 우리 틀 유닛
#     ⇒ X 가 크게 보인다 — 위아래 두 상자의 폭이 같기 때문이다
#        (앞 판에서 X 가 안 보였던 까닭은 위 꼭지가 6.4 pt 뿐이었기 때문)
#
#   왜 둘의 방향이 반대인가 (실측 2026-10-26)
#     · CHM13 `chr1:227,740,000-228,030,000` blastn : 5S 유전자 127 중 126 이 minus 가닥
#       (이웃 유전자 간격 2,237-2,240 bp ⇒ 우리 2,240 bp 틀과 맞는다)
#     · flank 창 : proximal `chr1:227,704,289-227,744,288`(낮음) ·
#                 distal  `chr1:228,026,263-228,066,262`(높음)
#       ⇒ 전사 = chr1 좌표 감소 = distal -> proximal = 그림에서 오른쪽 -> 왼쪽
#     · 우리 틀 : AluY 1-318 … 5'-flanking 2032-2121 · 5S gene 2122-2240 (TSS 2122)
#       ⇒ 틀 증가 = 전사 방향 ⇒ 틀을 1->2240 왼->오른쪽으로 그리면 유전체와 반대가 된다
#     · 원고 : 379 행 *"Right to left is the direction of transcription"* ·
#             422-423 행 *"the proximal end of the array is drawn on the left"*
#     ⇒ 유전체 방향 유닛은 gene 왼쪽 · AluY 오른쪽 · 전사 왼쪽
#        우리 틀 유닛은   AluY 왼쪽 · gene 오른쪽 · 전사 오른쪽
#
#   입력이 현행 산출물이다 — 일러스트레이터 원본 `fig1.pdf`(md5 150ddf302139)가
#     `$RDNA_DATA_ROOT`·`$WORK_ROOT`·`$REPO_ROOT_MIRROR` 깊이 6 에서 0 건이라 본 글을 다시 돌릴 수 없다.
#     ⇒ 이 글은 덧칠이다. 원본이 돌아오면 본 글에 합쳐야 한다.
import os
import sys
import fitz

sys.stdout.reconfigure(encoding="utf-8")

SRC = os.environ["FIG1_CURRENT_PDF"]
OUT = os.environ["FIG1_MARKED_PDF"]

# ── 실측 자리 (현행 PDF crop 좌표 · page 423.5 x 559.84) ────────────────────
CHR_X0 = 37.00
BAR_Y0, BAR_Y1 = 43.74, 52.00        # 배열 막대
SRC_X0, SRC_X1 = 228.52, 234.92      # 확대되는 칸
BOX_X0, BOX_X1 = 138.33, 346.09      # 유닛 블록 좌우 (두 층 공통)
ZFILL = (0.949, 0.961, 0.976)
EDGE = (0.29, 0.478, 0.69)
INK = (0.102, 0.102, 0.102)
GREY = (0.361, 0.380, 0.412)
LEAD = (0.690, 0.416, 0.416)
XLINE = (0.302, 0.349, 0.447)         # X 교차선 — 더 진하게 (전 0.55/0.60/0.68)
FS_EL, FS_TICK, FS_RU = 4.70, 4.70, 5.10
UTOT = 2240.0
MIR = BOX_X0 + BOX_X1

# 새 세로 배치 (배열 아래 52.0 ~ 패널 B 의 "B" 91.69 사이를 나눈다)
FAN_Y0, FAN_Y1 = 56.60, 62.00        # 확대선 (배열 아래 화살표 자리를 비운다)
GEN_Y0, GEN_Y1 = 62.00, 68.30        # 유전체 방향 유닛 (색만 · 6.30 pt)
XB_Y0, XB_Y1 = 72.30, 77.70          # X 교차 띠 — 유전체유닛 화살표(68.3->70.9) 아래에서 시작
XB_IN = 8.0                           # X 선 양 끝을 모서리에서 들인다 (왼쪽 화살표와 겹치지 않게)
OUR_Y0, OUR_Y1 = 77.70, 84.00        # 우리 틀 유닛 (이름 들어감 · 6.30 pt = 위 층과 같게)
TICK_BASE = 91.50                    # 눈금 바닥선 (패널 B 잉크는 y 95.4 부터 — 실측)

ELEM = [("AluY", 1, 318, (0.725, 0.804, 0.894)),
        ("spacer-1", 319, 1178, (0.890, 0.925, 0.965)),
        ("GA-rich", 1179, 1551, (0.953, 0.875, 0.722)),
        ("spacer-2", 1552, 2031, (0.890, 0.925, 0.965)),
        ("5'-flanking", 2032, 2121, (0.835, 0.878, 0.812)),
        ("5S rDNA gene", 2122, 2240, (0.918, 0.722, 0.682))]
INSIDE = {"AluY", "spacer-1", "GA-rich", "spacer-2"}
TICKS = [1, 319, 1179, 1552, 2240]


def xf(bp):
    """틀 좌표 -> 우리 틀 유닛의 x (1 이 왼쪽)"""
    return BOX_X0 + (BOX_X1 - BOX_X0) * (bp - 1) / UTOT


def xm(bp):
    """틀 좌표 -> 유전체 방향 유닛의 x (거울)"""
    return MIR - xf(bp)


def tlen(t, sz, bold=False):
    return fitz.get_text_length(t, fontname="hebo" if bold else "helv", fontsize=sz)


def bent(page, x_tss, y_base, rise, run, w=0.85, head=2.30, hw=1.10, color=INK):
    """TSS 꺾인 화살표 — 세로로 섰다가 전사 방향으로 (run<0 이면 왼쪽).
    머리는 채운 삼각형 (선 두 개짜리 V 는 작은 크기에서 깨진다)."""
    yt = y_base - rise
    s = 1.0 if run > 0 else -1.0
    xe = x_tss + run                 # 화살촉 끝
    xs = xe - s * head               # 샤프트가 끝나고 머리가 시작하는 자리
    sh = page.new_shape()
    sh.draw_polyline([fitz.Point(x_tss, y_base), fitz.Point(x_tss, yt), fitz.Point(xs, yt)])
    sh.finish(color=color, width=w, closePath=False, lineJoin=1)
    sh.commit()
    tr = page.new_shape()
    tr.draw_polyline([fitz.Point(xe, yt), fitz.Point(xs, yt - hw),
                      fitz.Point(xs, yt + hw), fitz.Point(xe, yt)])
    tr.finish(color=color, fill=color, width=0.0, closePath=True)
    tr.commit()


def unit_row(page, y0, y1, xpos, labels, outside_left):
    """유닛 한 줄 — xpos(bp)->x · labels=True 면 이름을 넣는다"""
    page.draw_rect(fitz.Rect(BOX_X0 - 1.0, y0 - 1.0, BOX_X1 + 1.0, y1 + 1.0),
                   color=None, fill=(1, 1, 1))
    for nm, a, b, col in ELEM:
        xa, xb = sorted((xpos(a), xpos(b + 1)))
        page.draw_rect(fitz.Rect(xa, y0, xb, y1), color=None, fill=col)
    page.draw_rect(fitz.Rect(BOX_X0, y0, BOX_X1, y1), color=EDGE, width=0.8)
    for nm, a, b, col in ELEM[1:]:
        page.draw_line(fitz.Point(xpos(a), y0), fitz.Point(xpos(a), y1), color=EDGE, width=0.5)
    if not labels:
        return
    base = y0 + (y1 - y0) / 2.0 + FS_EL * 0.36
    for nm, a, b, col in ELEM:
        if nm not in INSIDE:
            continue
        cx = (xpos(a) + xpos(b + 1)) / 2.0
        page.insert_text((cx - tlen(nm, FS_EL) / 2.0, base), nm,
                         fontname="helv", fontsize=FS_EL, color=INK)
    # 좁은 둘은 바깥에 이음선으로
    for (bp_in, txt, dy) in ((2045, "5'-flanking", -1.9), (2190, "5S rDNA gene", +2.6)):
        bx = xpos(bp_in)
        ex = (BOX_X1 + 7.0) if not outside_left else (BOX_X0 - 7.0)
        ey = (y0 + y1) / 2.0 + dy * 1.35
        page.draw_line(fitz.Point(bx, y0 if dy < 0 else y1), fitz.Point(ex, ey),
                       color=LEAD, width=0.70)
        tx = (ex + 2.2) if not outside_left else (ex - 2.2 - tlen(txt, FS_EL))
        page.insert_text((tx, ey + FS_EL * 0.36), txt, fontname="helv",
                         fontsize=FS_EL, color=INK)


doc = fitz.open(SRC)
page = doc[0]
print("입력 %s (%.1f x %.1f pt)" % (os.path.basename(SRC), page.rect.width, page.rect.height))

# ── ① 옛 패널 A 아래 절반(확대선·유닛·눈금·바깥라벨·repeat unit)을 지운다 ──────
for r in (fitz.Rect(60.0, 52.30, 410.0, 90.0),
          fitz.Rect(64.0, 69.5, 120.0, 78.5)):
    page.draw_rect(r, color=None, fill=(1, 1, 1))

# ── ② 확대선 (배열의 한 칸 -> 유전체 방향 유닛) ────────────────────────────
fan = page.new_shape()
fan.draw_polyline([fitz.Point(SRC_X0, FAN_Y0), fitz.Point(SRC_X1, FAN_Y0),
                   fitz.Point(BOX_X1, FAN_Y1), fitz.Point(BOX_X0, FAN_Y1),
                   fitz.Point(SRC_X0, FAN_Y0)])
fan.finish(color=None, fill=ZFILL)
fan.commit()

# ── ③ 유전체 방향 유닛 — 색만 (이름·눈금 없음) ──────────────────────────
unit_row(page, GEN_Y0, GEN_Y1, xm, labels=False, outside_left=True)

# ── ④ X 교차 — 위아래 상자의 폭이 같아 교차가 크게 보인다 ──────────────
for (xa, xb) in ((BOX_X0 + XB_IN, BOX_X1 - XB_IN), (BOX_X1 - XB_IN, BOX_X0 + XB_IN)):
    page.draw_line(fitz.Point(xa, XB_Y0), fitz.Point(xb, XB_Y1),
                   color=XLINE, width=0.55, dashes="[1.8 1.5] 0")   # 점선 · 더 얇게

# ── ⑤ 우리 틀 유닛 — 이름·눈금 ──────────────────────────────────────
unit_row(page, OUR_Y0, OUR_Y1, xf, labels=True, outside_left=False)
for v in TICKS:
    cx = xf(v)
    w = tlen(str(v), FS_TICK)
    if v == 2240:
        cx -= 8.60
    page.insert_text((cx - w / 2.0, TICK_BASE), str(v), fontname="helv",
                     fontsize=FS_TICK, color=GREY)
page.insert_text((BOX_X0 - 7.5 - tlen("repeat unit, 2,240 bp", FS_RU, True),
                  TICK_BASE), "repeat unit, 2,240 bp",
                 fontname="hebo", fontsize=FS_RU, color=INK)

# ── ⑥ 전사 방향 — 세 곳 모두 ─────────────────────────────────────────
#   배열 : 한 칸 안에서 틀 1 은 오른쪽 끝 (전사가 proximal=왼쪽 쪽이므로)
TSS_BAR = SRC_X1 - (2121.0 / UTOT) * (SRC_X1 - SRC_X0)
bent(page, TSS_BAR, BAR_Y1, -2.30, -(SRC_X1 - SRC_X0))  # 막대 아래 · 가로 6.40 pt
page.insert_text((TSS_BAR + 2.2, 56.00), "transcription",
                 fontname="helv", fontsize=FS_EL, color=INK)
#   유전체 방향 유닛 : TSS = xm(2122) · 왼쪽
bent(page, xm(2122), GEN_Y1, -2.60, -7.0)   # 상자 아래
#   우리 틀 유닛 : TSS = xf(2122) · 오른쪽 (틀이 유전체와 반대이므로)
bent(page, xf(2122), OUR_Y1, -2.25, 7.0)    # 상자 아래

# ── ⑦ chr1 왼쪽 끝이 1 번 자리 ───────────────────────────────────────
page.draw_line(fitz.Point(CHR_X0, 22.4), fitz.Point(CHR_X0, 24.6), color=INK, width=0.7)
page.insert_text((CHR_X0 - tlen("1", 6.0) / 2.0, 30.2), "1",
                 fontname="helv", fontsize=6.0, color=INK)

tmp = OUT + ".tmp"
doc.save(tmp, garbage=3, deflate=True)
doc.close()
os.replace(tmp, OUT)
print("저장 :", OUT)
print("유전체유닛 y %.1f-%.1f · X 띠 %.1f-%.1f · 우리유닛 %.1f-%.1f · 눈금 %.1f"
      % (GEN_Y0, GEN_Y1, XB_Y0, XB_Y1, OUR_Y0, OUR_Y1, TICK_BASE))
