# -*- coding: utf-8 -*-
# Fig 1 — 사용자가 논문 규격으로 배치한 fig1.ai 위의 패치
#   원본(.ai)에서 매번 새로 만든다 — 덧칠이 쌓이지 않는다.
#
#   ① 8pt 미만 글자를 8pt 로 (눈금 6 · 줄임표 2)
#   ② 어색한 글자 : "Hifiasm -ont" -> "hifiasm --ont" · "array pos" -> "array position"
#   ③ v11 참조서열 표 (코호트 + ONT read 출처)
#   ④ 겹침 고침
#        - 눈금 2122/2240 이 겹쳐 2240 이 두 번째 줄로 내려가 있었다 -> 2122 를 빼고 한 줄로
#        - 표 제목이 조립본 막대(y197-206)를 파고들었다 -> 표를 아래로, 줄 간격을 좁게
#        - "array position 1 ... N" 이 array 막대(y271.6-288.6) 위에 얹혀 있었다 -> 막대 위로
#        - "5S gene" 앞 2.4pt 잇는 선이 하이픈처럼 붙어 보였다 -> 지움
#   ⑤ 위아래 정렬 : 같은 줄인데 어긋난 글자들을 한 바닥선으로
#        - 두 갈래의 "5S rDNA array" (왼 230.27 / 오른 232.66)
#        - 좌우 "Partial unit" (295.97 / 296.77)
#        - "5S gene" 과 "repeat unit, 2,240 bp"
#   heredoc 으로 만들지 말 것 — 역슬래시가 뭉개진다
import fitz, os, sys
sys.stdout.reconfigure(encoding="utf-8")

# -- (repo) paths come from the environment --------------------------------
#   기본값이 분석 때 쓴 경로와 같다 ⇒ 작업 폴더에서는 전과 똑같이 돈다.
#   그림 내용은 한 글자도 바뀜이 없다 — png 화소로 맞대어 확인한다.
WORK_ROOT = os.environ["WORK_ROOT"]
OUT_DIR = os.environ.get("OUT_DIR",
                         os.path.join(WORK_ROOT, "04_figures", "paper_figures_final_2026-10-07",
                                      "01_selected_for_paper"))
HERE = os.path.dirname(os.path.abspath(__file__))
# `fig1.pdf` 는 코드가 만드는 것이 아니라 입력이다 —
#   사용자가 Illustrator 로 배치한 원본(.ai)을 pdf 로 내보낸 것(md5 150ddf302139).
#   자료와 함께 기탁한다. `FIG1_SCHEMATIC_PDF` 로 자리를 가리킬 수 있다.
SRC = os.environ.get("FIG1_SCHEMATIC_PDF", os.path.join(HERE, "fig1.pdf"))
OUT = os.path.join(OUT_DIR, "Fig1_pipeline_EN.pdf")

doc = fitz.open(SRC)
page = doc[0]
FH, FB = fitz.Font("helv"), fitz.Font("hebo")
# ── 일러스트레이터 안내선 4 개 지우기 (사용자 2026-10-19) ──
#   판 밖까지 뻗는 선이다 — x 72.0/526.0 (y -21~862) · y 70.9/771.0 (x -15~610).
#   그린 내용은 x 85.7~519.0 · y 82.6~666 이라 네 줄은 빈 여백에만 있다.
#   ⇒ 가는 띄로 그림요소를 지우면 안내선만 간다. (확인 : 지운 뒤 판을 가로지르는 선 0)
_GUIDES = [fitz.Rect(70.5, -30.0, 74.5, 870.0), fitz.Rect(524.5, -30.0, 528.5, 870.0),
           fitz.Rect(-20.0, 69.4, 615.0, 72.4), fitz.Rect(-20.0, 769.5, 615.0, 772.5),
           # 고아 눈금 하나 더 — 폭 0 · 높이 2.9 pt · 내용(509.4)보다 바깥
           fitz.Rect(530.0, 426.0, 536.0, 434.0)]
for _g in _GUIDES:
    page.add_redact_annot(_g)
page.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE,
                      graphics=fitz.PDF_REDACT_LINE_ART_REMOVE_IF_TOUCHED,
                      text=fitz.PDF_REDACT_TEXT_NONE)
_left = [d for d in page.get_drawings()
         if d["rect"].height > page.rect.height * 0.95 or d["rect"].width > page.rect.width * 0.95]
print("안내선 지움 — 남은 판 가로지르는 선 : %d 개" % len(_left))
assert not _left, "안내선이 남았다"


INK, GREY = (0.102, 0.102, 0.102), (0.36, 0.38, 0.41)
DESC = 0.299   # 다시 넣은 글자의 descent 보정 (실측 : 8pt 에서 2.39pt)



def rgb(c):
    if isinstance(c, int):
        return ((c >> 16 & 255) / 255.0, (c >> 8 & 255) / 255.0, (c & 255) / 255.0)
    return INK


def spans():
    out = {}
    for b in page.get_text("dict")["blocks"]:
        for l in b.get("lines", []):
            for s in l["spans"]:
                t = s["text"].strip()
                if t:
                    out[(t, round(s["bbox"][0], 1), round(s["bbox"][1], 1))] = dict(
                        rect=fitz.Rect(s["bbox"]), size=round(s["size"], 2),
                        font=s["font"], color=s.get("color", 0))
    return out


BEFORE = spans()

# ── b 의 왼쪽 머리글에서 `HPRC` 를 뺀다 (사용자 2026-10-11) ────────────
#   까닭 : 오른쪽 짝이 `Individuals without a public assembly` 라 `HPRC` 가 한쪽에만 붙어
#          "HPRC 조립본이 없는 사람" 과 "공개 조립본이 없는 사람" 이 다른 말처럼 읽힌다.
#   자리 : 왼쪽맞춤으로 두면 단 하나의 기둥 가운데에서 어긋난다 -> 옛 가운데를 지킨다.
_IHA = "Individuals with a public HPRC assembly"
_IHB = "Individuals with a public assembly"
_IHX = 133.9 + FB.text_length(_IHA, 9.0) / 2.0 + FB.text_length(_IHB, 9.0) / 2.0   # 오른끝

# (글자, x0) -> (새 글자, 바닥선 보정 dy, 오른쪽맞춤 x 또는 None, 새 크기 또는 None)
CHANGE = {
    # dy -2.4 : 다시 넣은 글자만 descent 보정만큼 내려앉아 오른쪽 짝(y 176.1)과
    #           바닥선이 어긋났다(실측 178.5). 짝이므로 같은 줄로 되돌린다.
    (_IHA, 133.9): (_IHB, -2.4, _IHX, None),
    # 패널 글자 a·b 만 9.0 pt 였다 — c~e 는 11.5 pt (2026-10-13 실측).
    #   한 그림 안에서 패널 글자 크기가 두 가지면 안 된다 -> 둘 다 11.5 pt · 왼끈은 c 와 같은 85.7
    ("a", 84.0): ("a", 0.0, 85.7 + FB.text_length("a", 11.5), 11.5),
    ("b", 87.1): ("b", 0.0, 85.7 + FB.text_length("b", 11.5), 11.5),
    # 유닛 모식도는 아래에서 6 구획으로 다시 그린다 -> 옛 글자는 지우기만 한다
    ("1", 226.3):    (None, 0, None, None),
    ("319", 254.6):  (None, 0, None, None),
    ("1179", 335.4): (None, 0, None, None),
    ("1552", 370.7): (None, 0, None, None),
    ("2122", 424.8): (None, 0, None, None),
    ("2240", 436.1): (None, 0, None, None),
    ("AluY", 235.4): (None, 0, None, None),
    ("spacer-1", 286.7): (None, 0, None, None),
    ("GA-rich", 344.2): (None, 0, None, None),
    ("spacer-2", 388.0): (None, 0, None, None),
    ("\u2026", 376.0): ("...", 0.0, None, 8.0),
    ("\u2026", 297.7): ("...", 0.0, None, 8.0),
    ("Hifiasm -ont", 403.4): ("hifiasm --ont", 0.0, None, None),
    ("array pos 1 \u2026 N", 269.0): (None, 0, None, None),   # 사용자 : 이 글자 지운다
    # the numbering order is reversed
    #   실측 : `array_pos 1` 은 원위다 — `hg002` 를 다리로 다른 파이프라인의 유전체 좌표와
    #   대어 보니 유닛 길이열이 뒤집혀 맞는다(pat 68/68 · mat 101/103 ↔ 같은 차례 21·22).
    #   ⇒ 그림 규약(근위는 왼쪽)은 지키고 번호만 오른쪽(원위)에서 시작하게 적는다.
    #   글자 수가 같아 자리는 안 움직인다. 코드 `scratchpad\close_direction.py`
    # 같은 줄로 맞출 짝은 둘 다 다시 넣는다 (한쪽만 넣으면 descent 만큼 어긋난다)
    # 첫 글자가 제 막대 가운데에서 3.03 pt 오른쪽이었다(실측
    #   막대 180.37~256.82 가운데 218.60 · 글자 가운데 221.63) -> 가운데로
    ("5S rDNA array", 195.8): ("5S rDNA array", 0.0,
                               218.60 + FH.text_length("5S rDNA array", 8.0) / 2.0, None),
    ("5S rDNA array", 192.2): ("5S rDNA array", 0.0, None, None),     # 왼쪽 갈래 (기준)
    ("5S rDNA array", 366.2): ("5S rDNA array", -2.39, None, None),   # 오른쪽 갈래
    ("Partial unit", 115.7): ("Partial unit", 0.0, None, None),       # 왼쪽 (기준)
    ("Partial unit", 472.2): ("Partial unit", -0.80, None, None),     # 오른쪽
    ("repeat unit, 2,240 bp", 137.0): ("repeat unit, 2,240 bp", 0.0, None, None),  # 기준
    ("5S gene", 444.4): (None, 0, None, None),
}
TICK = {"1", "319", "1179", "1552", "2240"}

tg = {}
for (t, x, y), v in BEFORE.items():
    if (t, round(x, 1)) in CHANGE:
        tg[(t, round(x, 1))] = v
miss = [k for k in CHANGE if k not in tg]
assert not miss, "못 찾음 %s" % miss

for v in tg.values():
    page.add_redact_annot(v["rect"])
page.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE,
                      graphics=fitz.PDF_REDACT_LINE_ART_NONE,
                      text=fitz.PDF_REDACT_TEXT_REMOVE)

# 안전망 : 의도치 않게 깎인 이웃 글자를 되살린다
AFTER = spans()
gone = []
for (t, x, y), v in list(AFTER.items()):
    for (t0, x0, y0), v0 in BEFORE.items():
        if (t0, round(x0, 1)) in CHANGE:
            continue
        if abs(x - x0) < 0.6 and abs(y - y0) < 0.6 and t != t0:
            gone.append(((t0, x0, y0), v0))
            page.add_redact_annot(v["rect"])
for (t, x, y), v in BEFORE.items():
    if (t, round(x, 1)) in CHANGE:
        continue
    if (t, x, y) not in AFTER and not any(g[0] == (t, x, y) for g in gone):
        gone.append(((t, x, y), v))
if gone:
    page.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE,
                          graphics=fitz.PDF_REDACT_LINE_ART_NONE,
                          text=fitz.PDF_REDACT_TEXT_REMOVE)
    for (t, x, y), v in gone:
        page.insert_text((v["rect"].x0, v["rect"].y1 - DESC * v["size"]), t,
                         fontname="hebo" if "Bold" in v["font"] else "helv",
                         fontsize=v["size"], color=rgb(v["color"]))
print("되살린 글자 %d : %s" % (len(gone), [g[0][0] for g in gone]))

for k, (new, dy, rightx, newsz) in CHANGE.items():
    if new is None:
        continue
    v = tg[k]
    fn = "hebo" if "Bold" in v["font"] else "helv"
    f = FB if fn == "hebo" else FH
    sz = newsz if newsz else v["size"]
    r = v["rect"]
    if rightx is not None:
        x = rightx - f.text_length(new, sz)
    elif k[0] in TICK or k[0] == "2122":
        x = (r.x0 + r.x1) / 2.0 - f.text_length(new, sz) / 2.0
    else:
        x = r.x0
    page.insert_text((x, r.y1 - DESC * sz + dy), new, fontname=fn, fontsize=sz, color=rgb(v["color"]))

# "5S gene" 앞 2.4pt 잇는 선 — 하이픈처럼 붙어 보인다. 흰 선으로 덮는다
page.draw_line(fitz.Point(440.9, 154.38), fitz.Point(444.6, 154.38),
               color=(1, 1, 1), width=1.6)

# ── 유닛 모식도를 합의된 6 구획 으로 다시 그린다 ─────────────────
# : spacer-2 와 5S gene 사이에 5'-flanking 이 있는 구조가 약속된 것
# boundary source : `s01_element_bounds.tsv` (TSS 2,123)
#   주의 : 그 표의 유전자 끝은 2,241 이라 합이 2,241 bp 가 된다. 우리 유닛은 2,240 bp 이므로
#          마지막 구획 끝만 2,240 으로 두었다(1 bp 차이는 사용자 확인 대기).
UX0, UX1 = 227.23, 439.80          # 기존 유닛 상자의 좌우 끝 (그대로 쓴다)
UY0, UY1 = 149.90, 158.85
# 우리 틀 : 유닛 2,240 bp · TSS 2122 (합계 정확히 2,240)
#    표(2,241bp · TSS 2123)와는 GA 뒤가 1 bp 밀린 판 차이 이지 오류가 아니다
#   (근거 figlib.py:46 · mk_C_offstate.py:32)
ELEM = [("AluY", 1, 318, (0.725, 0.804, 0.894)),
        ("spacer-1", 319, 1178, (0.890, 0.925, 0.965)),
        ("GA-rich", 1179, 1551, (0.953, 0.875, 0.722)),
        ("spacer-2", 1552, 2031, (0.890, 0.925, 0.965)),
        ("5'-flanking", 2032, 2121, (0.835, 0.878, 0.812)),
        ("5S rDNA gene", 2122, 2240, (0.918, 0.722, 0.682))]
UTOT = 2240.0
EDGE = (0.29, 0.478, 0.69)


def uxp(bp):
    return UX0 + (UX1 - UX0) * (bp - 1) / UTOT


# 옛 상자·글자를 흰색으로 덮는다 (글자는 위에서 이미 지웠다)
page.draw_rect(fitz.Rect(UX0 - 1.2, UY0 - 1.2, UX1 + 1.2, UY1 + 1.2),
               color=None, fill=(1, 1, 1))
for nm, a, b, col in ELEM:
    r = fitz.Rect(uxp(a), UY0, uxp(b + 1), UY1)
    page.draw_rect(r, color=None, fill=col)
page.draw_rect(fitz.Rect(UX0, UY0, UX1, UY1), color=EDGE, width=0.8)
for nm, a, b, col in ELEM[1:]:
    page.draw_line(fitz.Point(uxp(a), UY0), fitz.Point(uxp(a), UY1), color=EDGE, width=0.5)

# 이름 : 넓은 구획은 상자 안, 좁은 둘은 오른쪽 바깥에 이음선으로
# 상자 높이가 8.95 pt 뿐이라 8 pt 글자의 머리가 위로 삐져나왔다(실측 1.08 pt).
#   머리(cap 5.736)와 꼬리(desc 1.68)가 양쪽에서 같은 만큼 들어가는 기준선을 쓴다.
BASE = (max(UY0 + 0.15 + 0.717 * 8.0, UY1 - 0.15 - 0.21 * 8.0)
        + min(UY1 - 0.15 - 0.21 * 8.0, UY0 + 0.15 + 0.717 * 8.0)) / 2.0
assert UY0 <= BASE - 0.717 * 8.0 and BASE + 0.21 * 8.0 <= UY1, (
    '구획 이름이 상자를 벗어난다 : 머리 %.2f (UY0 %.2f) · 꼬리 %.2f (UY1 %.2f)'
    % (BASE - 0.717 * 8.0, UY0, BASE + 0.21 * 8.0, UY1))
for nm, a, b, col in ELEM:
    w = uxp(b + 1) - uxp(a)
    tw = FH.text_length(nm, 8.0)
    if tw + 2.0 < w:
        page.insert_text(((uxp(a) + uxp(b + 1)) / 2.0 - tw / 2.0, BASE), nm,
                         fontname="helv", fontsize=8.0, color=INK)
OUTX = UX1 + 7.0
for k, (nm, a, b, col) in enumerate([e for e in ELEM if FH.text_length(e[0], 8.0) + 2.0 >= uxp(e[2] + 1) - uxp(e[1])]):
    yy = UY0 - 1.0 + k * 9.0
    page.draw_line(fitz.Point((uxp(a) + uxp(b + 1)) / 2.0, UY1 if k else UY0),
                   fitz.Point(OUTX - 1.5, yy - 2.2), color=(0.55, 0.57, 0.6), width=0.5)
    page.insert_text((OUTX, yy), nm, fontname="helv", fontsize=8.0, color=INK)

# 눈금 : 좁은 둘의 경계는 글자가 겹치므로 넓은 경계만
TICKY = UY1 + 7.2
for bp in (1, 319, 1179, 1552, 2240):
    lab = "%d" % bp
    tw = FH.text_length(lab, 8.0)
    x = uxp(bp) if bp != 2240 else UX1
    page.insert_text((min(x - tw / 2.0, UX1 - tw), TICKY), lab,
                     fontname="helv", fontsize=8.0, color=(0.36, 0.38, 0.41))
print("유닛 모식도 : 6 구획으로 다시 그림 (우리 틀 2,240 · TSS 2122 · 5'-flanking 2,032-2,121)")

# ── v11 참조서열 표는 뺐다 (사용자 2026-10-10 : *"b 의 array references 를
#     e 와 합쳐줘 · 그림이 아니라 표로"*) -> 새 패널 c `mk_Fig1_panelC_TABLE.py`
#   여기서 그리지 않는다. 원본 .ai 에도 그 표는 없다(내가 그려 넣었던 것이다).

# ── 주황 array 막대 정렬 ──
#   실측 : 회색 받침 x0 가 256.83 / 257.58 / 259.25 로 2.42 pt 어긋나 있었다.
#          (폭은 셋 다 같다 — 회색 87.24 · 주황 76.45) ⇒ 자리만 맞춘다.
#   방법 : 어긋난 둘을 흰색으로 덮고 같은 크기로 다시 그린다. 원본 .ai 는 건드리지 않는다.
GX, GW, GH = 256.83, 87.24, 8.60          # 회색 받침 : 기준 x0 · 폭 · 높이
OIN_L, OIN_R, OIN_T = 5.01, 5.79, 0.25    # 주황이 회색 안으로 들어간 양 (실측)
OH = 7.86
GREYB, ORANGE = (0.89, 0.89, 0.89), (0.91, 0.64, 0.39)
for y0, oldx in ((242.85, 257.58), (260.75, 259.25)):
    page.draw_rect(fitz.Rect(oldx - 1.0, y0 - 1.0, oldx + GW + 1.0, y0 + GH + 1.0),
                   color=None, fill=(1, 1, 1))
    page.draw_rect(fitz.Rect(GX, y0, GX + GW, y0 + GH),
                   color=None, fill=GREYB, radius=0.30)
    page.draw_rect(fitz.Rect(GX + OIN_L, y0 + OIN_T, GX + GW - OIN_R, y0 + OIN_T + OH),
                   color=None, fill=ORANGE, radius=0.32)
print("주황 막대 정렬 : x0 257.58 / 259.25 -> %.2f (폭 %.2f 그대로)" % (GX, GW))

# ── "5S rDNA arrays" 를 그 무리의 설명으로 옮긴다 (사용자 ① 지시) ──
#   지금 x 179.2-234.7 y 234.7 — 표 오른쪽 끝(174.3)과 5 pt 밖에 안 떨어져 붙어 보인다.
#   무리는 y 231.4-269.4 (세로 가운데 250.4) 이고 왼쪽 끝이 GX=256.83 이다.
#   이 무리를 정리한 것이 패널 c 이므로 같은 이름으로 쓴다.
#     패널 c 제목 = "Array references and cohort tiers" ⇒ 여기 라벨도 `Array references`.
_LB = "Array references"
_LW = FH.text_length(_LB, 8.0)
page.add_redact_annot(fitz.Rect(177.0, 230.0, 237.0, 241.0))
page.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE)
_LX = GX - 6.0 - _LW                       # 무리 왼쪽에서 6 pt 띄운다
page.insert_text((_LX, 250.4 + 8.0 * 0.5 - 0.299 * 8.0), _LB,
                 fontname="helv", fontsize=8.0, color=INK)
print("라벨 '5S rDNA arrays' : x %.1f-%.1f y 250.4 (표 끝 174.3 에서 %.1f pt 띄움)"
      % (_LX, _LX + _LW, _LX - 174.3))

# ── 확대 표시 : 주황 막대 -> 파란 unit 열 ──
#   실측 : 출처 막대 x 256.83-344.08 (아래끝 y 269.35) · 파란 열 x 120.50-514.71 (윗끝 y 285.0)
#          기존 쐐기 꼭짓점은 x 279.92-328.06 으로 막대보다 좁아 출처가 안 읽혔다.
ZX0, ZX1, ZY = 256.83, 344.08, 269.35          # 출처 막대 (맨 아래 것)
AX0, AX1, AY = 120.50, 514.71, 285.00          # 파란 unit 열
ZCOL = (0.62, 0.68, 0.76)
# 선을 그으면 "Proximal flank"(x 115.7-177.2) · "Distal flank"(x 471-519) 글자를
#   가로지른다 (y 272.5-280.5 에서 선이 x 170-174 를 지난다 — 실측).
#   ⇒ 선은 쓰지 않는다. 출처 막대에 테두리만 둘러 "이것을 확대했다" 를 못박고,
#     아래로 잇는 일은 이미 있는 옅은 쐐기(y 271.9-288.6)가 한다.
_ = (AX0, AX1, AY)          # 재어 둔 값 — 나중에 쐐기를 다시 그릴 때 쓴다
page.draw_rect(fitz.Rect(ZX0 - 0.9, 260.75 - 0.9, ZX1 + 0.9, ZY + 0.9),
               color=ZCOL, width=0.7, radius=0.22)
print("확대 표시 : 막대 x %.1f-%.1f -> 열 x %.1f-%.1f (선 2 + 테두리)" % (ZX0, ZX1, AX0, AX1))

# ── ② 확대 영역을 연한 하늘색 한 덩이로 (사용자 2026-10-10) ──────────
#   지금은 쐐기 두 장(x 124.5-279.9 · 328.1-509.4 · y 271.9-288.6)이라
#   삼각형 둘로 보인다. array 를 확대해서 본 것이므로 하나로 이어야 한다.
#   ⇒ 출처 막대(x 256.8-344.1 아래끝 269.4) 에서 파란 열(x 120.5-514.7 윗끝 285.0) 로
#     사다리꼴 하나를 그 쐐기와 같은 색으로 채운다.
ZFILL = (0.949, 0.961, 0.976)            # 기존 쐐기와 같은 색 (실측)
_zq = [fitz.Point(256.83, 269.9), fitz.Point(344.08, 269.9),
       fitz.Point(514.71, 288.6), fitz.Point(120.50, 288.6)]   # 288.6 = 옛 쐐기 아래끝 (284.6 이면 흰 틈이 남는다)
_sh = page.new_shape()
_sh.draw_polyline(_zq + [_zq[0]])
_sh.finish(color=None, fill=ZFILL, closePath=True)
_sh.commit(overlay=False)                 # 글자·막대 아래에 깔린다
print("확대 영역 : 사다리꼴 하나로 채움 (막대 256.8-344.1 -> 열 120.5-514.7)")

# ── ③ flank 로 뻗는 read 를 더 길게 (사용자 2026-10-10) ───────────────
#   실측 구조 : 회색(flank 쪽 9.6 pt) + 보라(array 쪽) + 흰 눈금
#   보라가 18-48 pt 인데 array 안 read(초록)는 57-147 pt 라 너무 짧아 보인다.
PUR, GRY = (0.44, 0.34, 0.57), (0.66, 0.66, 0.66)
#   길이는 array 안 read(초록 57-147 pt)와 견줄 만하게. 더 늘리면 라벨·초록을 침범한다(실측).
LEFTR = [(314.2, 143.3, 175.0), (316.7, 145.6, 193.0), (319.1, 152.4, 208.0),
         (321.5, 171.8, 222.0), (324.0, 173.7, 212.0)]
RIGHTR = [(325.4, 461.4, 420.0), (327.8, 478.1, 437.0), (330.2, 490.1, 452.0),
          (332.7, 492.2, 444.0), (335.1, 494.1, 460.0)]
for _y, _old, _new in LEFTR:
    page.draw_rect(fitz.Rect(125.0, _y - 1.5, max(_old, _new) + 1.5, _y + 1.5),
                   color=None, fill=(1, 1, 1))
    page.draw_line(fitz.Point(125.3, _y), fitz.Point(_new, _y), color=PUR, width=1.7)
for _y, _old, _new in RIGHTR:
    page.draw_rect(fitz.Rect(min(_old, _new) - 1.5, _y - 1.5, 509.7, _y + 1.5),
                   color=None, fill=(1, 1, 1))
    page.draw_line(fitz.Point(_new, _y), fitz.Point(509.4, _y), color=PUR, width=1.7)
# 유닛 눈금(흰 선)을 다시 찍는다 — 2,240 bp 가 파란 열에서 차지하는 폭으로
#   사용자 : *"늘인 부분에 흰 선이 없다"* — 실측으로 끊겨 있었다.
#   까닭 : 옛 식이 `_old`·`_new` 비교로 구간을 골라 왼쪽 무리 세 줄이 174.1 에서 멈췄다.
#   ⇒ read 마다 그 read 가 실제로 그려진 구간(왼쪽 125.3~_new · 오른쪽 _new~509.4)에
#     처음부터 끝까지 긋는다. 눈금 자리는 array 격자(120.50 + k·_UW)에 그대로 맞춘다.
_UW = (514.71 - 120.50) / 125.0           # CHM13 array 125 유닛 (내 실측)
for _y, _old, _new in LEFTR + RIGHTR:
    _a, _b = (125.3, _new) if _y < 325 else (_new, 509.4)   # 그려진 구간 그대로
    _a, _b = min(_a, _b), max(_a, _b)
    _k = int((_a - 120.50) / _UW) + 1
    _n = 0
    while True:
        _t = 120.50 + _k * _UW
        if _t >= _b - 0.4:
            break
        page.draw_line(fitz.Point(_t, _y - 1.5), fitz.Point(_t, _y + 1.5),
                       color=(1, 1, 1), width=0.55)
        _k += 1
        _n += 1
    print("   눈금 y %.1f : %2d 개 (x %.1f ~ %.1f)" % (_y, _n, _a, _b))
# 라벨 `Reads reaching into the flank` 가 read 띠(y 314-335)와 같은 높이라 묻힌다.
#   -> `Partial unit`(y 298-311) 과 같은 줄의 가운데 빈칸으로 올린다.
_LAB = "Reads reaching into the flank"
_lw = page.search_for(_LAB)
if _lw:
    _r = _lw[0]
    page.add_redact_annot(fitz.Rect(_r.x0 - 1, _r.y0 - 1, _r.x1 + 1, _r.y1 + 1))
    # `graphics` 를 안 주면 기본값이 닿은 밑그림을 지운다.
    #   이 라벨은 read 띠 한가운데(y 314~335)에 있었으므로, 바로 위에서 그린
    #   흰 유닛 눈금이 그 상자에 닿아 같이 지워졌다 — 왼쪽 2·3 번째 줄이
    #   x 174.1 뒤로 민둥이 된 까닭이다(사용자 2026-10-19 지적 · 로그는 21·26 개를
    #   그렸다고 적는데 결과물엔 17 개뿐이었다 ⇒ 그린 뒤에 지워진 것이다).
    #   ⇒ 글자만 지운다.
    page.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE,
                          graphics=fitz.PDF_REDACT_LINE_ART_NONE,
                          text=fitz.PDF_REDACT_TEXT_REMOVE)
    _fs, _bw = 9.0, FB.text_length(_LAB, 9.0)
    page.insert_text(((120.50 + 514.71) / 2.0 - _bw / 2.0, 309.0 - 0.299 * _fs),
                     _LAB, fontname="hebo", fontsize=_fs, color=(0.42, 0.33, 0.55))
    print("라벨 '%s' : y %.1f -> 309.0 · 가운데 정렬" % (_LAB, _r.y0))
else:
    print("label not found; position left as it is")
print("flank read : 보라 10 개를 길게 (왼쪽 143-174 -> 206-268 · 오른쪽 461-494 -> 344-404)")

# ── 모식도를 423.5 pt 로 (사용자 2026-10-19) ──
#   지금 85.7~519.0 = 433.3 pt · 패널(c~f)은 이미 423.5 라 모식도만 줄인다.
#   글자 규격화 앞에 해야 한다 — 뒤에 하면 4.7 이 4.59 로 내려간다.
_sx0, _sx1 = 85.7, 519.0
_S = 423.5 / (_sx1 - _sx0)
_sy0 = 82.6
_tx, _ty = 85.9 - _S * _sx0, _sy0 - _S * _sy0
_buf = doc.tobytes()
_src = fitz.open("pdf", _buf)
_nd = fitz.open()
_W, _H = page.rect.width, page.rect.height
page = _nd.new_page(width=_W, height=_H)
page.show_pdf_page(fitz.Rect(_tx, _ty, _tx + _S * _W, _ty + _S * _H), _src, 0)
_src.close()
doc.close()
doc = _nd
print("모식도 축소 : 433.3 -> 423.5 pt (배율 %.6f)" % _S)

# ── 모식도(a·b) 글자를 규격으로 ──
#   눈금 4.7 · 축 이름 5.1 · legend/주석 4.5 · 패널 글자 8.2 대문자
#   `y < 375` 만 — 그 아래는 matplotlib 패널이라 이미 규격이다.
#   바닥선은 원래 크기로 쟴 값을 쓴다 — 새 크기로 재면 글자가 위로 뜨다.
# 크기 목록으로 걸렀더니 하나를 놓쳤다 — page.get_text("dict") 은 8.495,
#   page.get_texttrace() 는 같은 글자를 8.02 로 적는다(변환행렬 해석이 다르다).
#   ⇒ 목록이 아니라 구간으로 가른다. (이 주석은 한 번 셸에 먹혔다 —
#      백틱이 든 글을 bash -c 로 넘기면 명령 치환으로 사라진다. Write/Edit 로 쓸 것)
def _specsz(v):
    return 8.2 if v >= 10.5 else (5.1 if v >= 8.8 else 4.7)
_UP = {"a": "A", "b": "B"}
_todo = []
for _b in page.get_text("dict")["blocks"]:
    for _l in _b.get("lines", []):
        for _sp in _l["spans"]:
            _t = _sp["text"].strip()
            if not _t:
                continue
            _r = fitz.Rect(_sp["bbox"])
            if _r.y1 >= 375.0:
                continue
            _k = round(_sp["size"], 2)
            if abs(_k - 4.5) < 0.05 or abs(_k - 4.7) < 0.05 or abs(_k - 5.1) < 0.05 or abs(_k - 8.2) < 0.05:
                continue
            _todo.append((_r, _t, _k, _sp.get("color", 0), _sp["font"]))
for _r, _t, _k, _c, _f in _todo:
    page.add_redact_annot(_r)
page.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE,
                      graphics=fitz.PDF_REDACT_LINE_ART_NONE,
                      text=fitz.PDF_REDACT_TEXT_REMOVE)
for _r, _t, _k, _c, _f in _todo:
    _nz = _specsz(_k)
    _bold = "Bold" in _f or "bold" in _f
    _F = FB if _bold else FH
    _new = _UP.get(_t, _t) if _nz == 8.2 else _t
    _base = _r.y1 - DESC * _k                 # 원래 크기로 재 바닥선
    if _nz == 8.2:
        _x = _r.x0                            # 패널 글자 — 왼끕 유지
    else:
        _x = (_r.x0 + _r.x1) / 2.0 - _F.text_length(_new, _nz) / 2.0
    page.insert_text((_x, _base), _new, fontname="hebo" if _bold else "helv",
                     fontsize=_nz, color=rgb(_c))
print("모식도 글자 규격화 : %d 개" % len(_todo))

# ── 패널 c + d : 하플로타입 수준만 · 한 줄 ──
#   만드는 코드 : mk_Fig1_panelCD.py  (먼저 돌려 두어야 한다)
_PT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_Fig1_panelC_TABLE.pdf")
if os.path.exists(_PT):
    _st = fitz.open(_PT)
    _pwt, _pht = _st[0].rect.width, _st[0].rect.height
    # d·e·f 를 위로 (사용자 2026-10-19) — 모식도 아래끝을 재서 바로 밑에 둔다.
    #   고정 378.0 을 쓰면 축소 뒤에 틈이 더 벌어진다.
    _sb = max([d["rect"].y1 for d in page.get_drawings() if d["rect"].y1 < 370.0]
              + [max(c[3][3] for c in sp["chars"] if chr(sp["chars"][0][0]))
                 for sp in page.get_texttrace()
                 if sp["chars"] and max(c[3][3] for c in sp["chars"]) < 370.0])
    _CT = _sb + 10.0
    print("\ubaa8\uc2dd\ub3c4 \uc544\ub798\ub05d %.1f \u2192 \ud328\ub110 c \uc704 %.1f" % (_sb, _CT))
    page.show_pdf_page(fitz.Rect(85.9, _CT, 85.9 + _pwt, _CT + _pht), _st, 0)
    _st.close()
    _CTOP = _CT + _pht + 8.0
    print("패널 c(표) 넣음 : y %.0f-%.0f (%.1f x %.1f pt)" % (_CT, _CT + _pht, _pwt, _pht))
else:
    _CTOP = 378.0
    print("패널 c 표 없음 — mk_Fig1_panelC_TABLE.py 를 먼저 돌리십시오")

_PC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_Fig1_panelCD.pdf")
if os.path.exists(_PC):
    _src = fitz.open(_PC)
    _pw, _ph = _src[0].rect.width, _src[0].rect.height
    _X0 = 85.9                      # 423.5 pt 내용폭을 A4 가운데로 ((595.3-423.5)/2)
    _TOP = _CTOP                    # 패널 c(표) 아래
    page.show_pdf_page(fitz.Rect(_X0, _TOP, _X0 + _pw, _TOP + _ph), _src, 0)
    _src.close()
    print("패널 c+d 넣음 : y %.0f-%.0f (%.1f x %.1f pt)" % (_TOP, _TOP + _ph, _pw, _ph))
else:
    print("패널 c+d 파일 없음 — mk_Fig1_panelCD.py 를 먼저 돌리십시오")

# 옛 패널 e(상자·화살표 그림)는 뺐다 — 사용자 2026-10-10 :
#   *"e 와 b 는 그림이 아니라 표로 합쳐줘"* ⇒ 새 패널 c `mk_Fig1_panelC_TABLE.py`
#   `mk_Fig1_panelE_COHORT.py` 와 `_Fig1_panelE.pdf` 는 지우지 않고 남겨 둔다(옛 판 보존).

# Fig1 은 확대하지 않는다 (사용자 2026-10-08 : *"fig1 은 이미 180 규격에 맞잖아"*)
#   이 그림은 사용자가 Illustrator 로 논문 규격에 맞춰 직접 만든 것이다.
#   다른 8 장처럼 판을 180 mm 로 바꾸는 작업을 여기에 적용하면 안 된다.
#   (내가 잰 내용 폭은 453 pt = 159.8 mm 였으나, 규격 판단은 사용자 것이다.)

# ── 내용 크기로 자르기 ──────────────────────
#   까닭 : 다른 그림은 전부 `423.5 x N` 인데 `Fig 1` 만 A4 라, 나란히 놓으면
#           `Fig 1` 만 작아 보인다 (위아래로 빈 띠가 크다).
#   `set_cropbox` 를 쓴다 — 다시 그리지도 배율을 바꾸지도 않고 보이는 상자만 바꾼다.
#     ⇒ 폭이 흔들릴 수 없다. `bbox_inches="tight"` 는 여백을 다시 계산해 폭이 흔들린다
# Besides, this figure is not a matplotlib output.
#   `x0` 는 재서 넣는다(85.9 로 박지 않는다) · 폭은 `x0 + 423.5` 로 떨어지게 둔다.
_CMAR = 2.0                      # 위아래 숨
_cx, _cy = [], []
for _d in page.get_drawings():
    _r = _d["rect"]
    if _r.width < 560 and _r.height < 800 and _r.x0 > -10:
        _cx += [_r.x0, _r.x1]; _cy += [_r.y0, _r.y1]
for _s in page.get_texttrace():
    _ch = [c for c in _s["chars"] if chr(c[0]).strip()]
    if _ch:
        _cx += [min(c[3][0] for c in _ch), max(c[3][2] for c in _ch)]
        _cy += [min(c[3][1] for c in _ch), max(c[3][3] for c in _ch)]
_CX0, _CY0, _CY1 = min(_cx), min(_cy), max(_cy)
assert abs((max(_cx) - _CX0) - 423.5) < 0.2, "내용 폭이 423.5 가 아니다 : %.2f" % (max(_cx) - _CX0)
page.set_cropbox(fitz.Rect(_CX0, _CY0 - _CMAR, _CX0 + 423.5, _CY1 + _CMAR))
print("내용 크기로 자름 : %.1f x %.1f pt  (내용 y %.1f~%.1f · 숨 %.1f)"
      % (page.rect.width, page.rect.height, _CY0, _CY1, _CMAR))

# 뷰어가 열고 있으면 덮어쓰기가 막힌다 -> 임시로 쓰고 바꿔치기
_tmp = OUT + ".tmp"
doc.save(_tmp, garbage=3, deflate=True)
doc.close()
os.replace(_tmp, OUT)
print("저장 :", OUT)
