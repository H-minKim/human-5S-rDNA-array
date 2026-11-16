# -*- coding: utf-8 -*-
# Fig 1 - patches applied on top of fig1.ai, which the author laid out to the journal's format
#   Built fresh from the original (.ai) every time, so patches do not accumulate.
#
#   (1) text below 8 pt raised to 8 pt (6 tick labels and 2 ellipses)
#   (2) awkward wording: "Hifiasm -ont" -> "hifiasm --ont"; "array pos" -> "array position"
#   (3) the v11 reference-sequence table (cohort and the source of the ONT reads)
#   (4) overlaps fixed
#        - the 2122 and 2240 ticks overlapped and 2240 had dropped to a second line -> 2122 removed and the rest kept on one line
#        - the table caption ran into the assembly bar (y197-206) -> the table moved down and the line spacing tightened
#@DROP
#@DROP
#@DROP
#@DROP
#@DROP
#@DROP
#@DROP
import fitz, os, sys
sys.stdout.reconfigure(encoding="utf-8")

# -- (repo) paths come from the environment --------------------------------
#   The roots are required; the output and table directories fall back to this file's folder.
#@DROP
WORK_ROOT = os.environ["WORK_ROOT"]
OUT_DIR = os.environ.get("OUT_DIR",
                         os.path.join(WORK_ROOT, "04_figures", "paper_figures_final_2026-10-07",
                                      "01_selected_for_paper"))
HERE = os.path.dirname(os.path.abspath(__file__))
# fig1.pdf is an input to this script, not something it produces:
#   it is the hand-laid-out Illustrator original exported to pdf (md5 150ddf302139).
#   It is deposited with the data, and FIG1_SCHEMATIC_PDF can point at it.
SRC = os.environ.get("FIG1_SCHEMATIC_PDF", os.path.join(HERE, "fig1.pdf"))
OUT = os.path.join(OUT_DIR, "Fig1_pipeline_EN.pdf")

doc = fitz.open(SRC)
page = doc[0]
FH, FB = fitz.Font("helv"), fitz.Font("hebo")
# -- remove the four Illustrator guide lines --
#   They run past the page edge: x 72.0 / 526.0 over y -21..862, and y 70.9 / 771.0 over x -15..610.
#   The drawn content spans x 85.7..519.0 and y 82.6..666, so all four lie in empty margin.
#   Deleting drawing elements along a thin band therefore removes the guides alone; afterwards no line crosses the page.
_GUIDES = [fitz.Rect(70.5, -30.0, 74.5, 870.0), fitz.Rect(524.5, -30.0, 528.5, 870.0),
           fitz.Rect(-20.0, 69.4, 615.0, 72.4), fitz.Rect(-20.0, 769.5, 615.0, 772.5),
           # One more orphan tick: width 0, height 2.9 pt, outside the content edge at 509.4.
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
DESC = 0.299   # Descent correction for re-inserted text; measured as 2.39 pt at 8 pt.



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

# -- drop HPRC from the left-hand heading of panel b --
#   Its counterpart on the right reads "Individuals without a public assembly", so HPRC on one side only
#          makes the two headings read as different claims.
#   Left-aligning it would move it off the centre of its column, so the original centre is kept.
_IHA = "Individuals with a public HPRC assembly"
_IHB = "Individuals with a public assembly"
_IHX = 133.9 + FB.text_length(_IHA, 9.0) / 2.0 + FB.text_length(_IHB, 9.0) / 2.0   # right edge

# (text, x0) -> (new text, baseline correction dy, right-aligned x or None, new size or None)
CHANGE = {
    # dy -2.4: only the re-inserted text sinks by the descent correction, putting its baseline
    #           out of line with its counterpart at y 176.1 (measured 178.5). The pair is returned to one baseline.
    (_IHA, 133.9): (_IHB, -2.4, _IHX, None),
    # The panel letters a and b were 9.0 pt while c to e were 11.5 pt.
    #   One figure must not carry two panel-letter sizes, so both become 11.5 pt with the left edge at 85.7, matching c.
    ("a", 84.0): ("a", 0.0, 85.7 + FB.text_length("a", 11.5), 11.5),
    ("b", 87.1): ("b", 0.0, 85.7 + FB.text_length("b", 11.5), 11.5),
    # The unit schematic is redrawn below in six segments, so the old text is only deleted here.
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
    ("array pos 1 \u2026 N", 269.0): (None, 0, None, None),   #@DROP
    # the numbering order is reversed
    #   Measured: array_pos 1 is the distal end. Bridging through hg002 to the genome coordinates of another
    #   pipeline, the unit-length series matches when reversed (paternal 68/68, maternal 101/103).
    #   The figure convention that proximal is on the left is kept, and only the numbering starts from the right, the distal end.
    #@DROP
    # A pair that must share a baseline is re-inserted together; inserting one alone offsets it by the descent.
    # The first label sat 3.03 pt right of the centre of its bar
    #   (bar 180.37..256.82, centre 218.60; text centre 221.63), so it is centred.
    ("5S rDNA array", 195.8): ("5S rDNA array", 0.0,
                               218.60 + FH.text_length("5S rDNA array", 8.0) / 2.0, None),
    ("5S rDNA array", 192.2): ("5S rDNA array", 0.0, None, None),     # left branch (the baseline reference)
    ("5S rDNA array", 366.2): ("5S rDNA array", -2.39, None, None),   # right branch
    ("Partial unit", 115.7): ("Partial unit", 0.0, None, None),       # left (the baseline reference)
    ("Partial unit", 472.2): ("Partial unit", -0.80, None, None),     # right
    ("repeat unit, 2,240 bp", 137.0): ("repeat unit, 2,240 bp", 0.0, None, None),  # reference
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

# A safety net that restores a neighbouring glyph removed by accident.
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

# The 2.4 pt leader before "5S gene" reads as a hyphen, so it is covered with a white line.
page.draw_line(fitz.Point(440.9, 154.38), fitz.Point(444.6, 154.38),
               color=(1, 1, 1), width=1.6)

# -- redraw the unit schematic in the agreed six segments --
# : the agreed structure places 5'-flanking between spacer-2 and the 5S gene
# boundary source : `s01_element_bounds.tsv` (TSS 2,123)
#   Note: that table ends the gene at 2,241, summing to 2,241 bp. The unit here is 2,240 bp, so
#          only the last segment ends at 2,240.
UX0, UX1 = 227.23, 439.80          # the left and right edges of the existing unit box, used unchanged
UY0, UY1 = 149.90, 158.85
# This frame: a 2,240 bp unit with the TSS at 2122; the segments sum to exactly 2,240.
#    Against the published table (2,241 bp, TSS 2123) everything after the GA is shifted by 1 bp; this is a difference of release, not an error.
#   (see figlib.py:46 and mk_C_offstate.py:32)
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


# The old box and labels are covered in white; the labels were already deleted above.
page.draw_rect(fitz.Rect(UX0 - 1.2, UY0 - 1.2, UX1 + 1.2, UY1 + 1.2),
               color=None, fill=(1, 1, 1))
for nm, a, b, col in ELEM:
    r = fitz.Rect(uxp(a), UY0, uxp(b + 1), UY1)
    page.draw_rect(r, color=None, fill=col)
page.draw_rect(fitz.Rect(UX0, UY0, UX1, UY1), color=EDGE, width=0.8)
for nm, a, b, col in ELEM[1:]:
    page.draw_line(fitz.Point(uxp(a), UY0), fitz.Point(uxp(a), UY1), color=EDGE, width=0.5)

# Labels: the wide segments are labelled inside the box, the two narrow ones outside on the right with a leader.
# The box is only 8.95 pt high, so the cap of 8 pt text overshot the top by 1.08 pt.
#   The baseline is chosen so that the cap height (5.736) and the descender (1.68) inset equally from both edges.
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

# Ticks: the labels would overlap at the two narrow boundaries, so only the wide boundaries are ticked.
TICKY = UY1 + 7.2
for bp in (1, 319, 1179, 1552, 2240):
    lab = "%d" % bp
    tw = FH.text_length(lab, 8.0)
    x = uxp(bp) if bp != 2240 else UX1
    page.insert_text((min(x - tw / 2.0, UX1 - tw), TICKY), lab,
                     fontname="helv", fontsize=8.0, color=(0.36, 0.38, 0.41))
print("유닛 모식도 : 6 구획으로 다시 그림 (우리 틀 2,240 · TSS 2122 · 5'-flanking 2,032-2,121)")

# -- The array-references table is not drawn here. It became panel c,
#     produced by mk_Fig1_panelC_TABLE.py.
#   The Illustrator original does not contain that table either.

# -- align the orange array bars --
#   Measured: the grey backing bars start at x0 256.83 / 257.58 / 259.25, 2.42 pt out of line.
#          All three are the same width (grey 87.24, orange 76.45), so only the position is corrected.
#   The two misaligned bars are covered in white and redrawn at the same size; the Illustrator original is untouched.
GX, GW, GH = 256.83, 87.24, 8.60          # grey backing bar: reference x0, width, height
OIN_L, OIN_R, OIN_T = 5.01, 5.79, 0.25    # how far the orange bar is inset within the grey one (measured)
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

# -- move "5S rDNA arrays" so that it labels the group it describes --
#   It sits at x 179.2-234.7, y 234.7, only 5 pt from the right edge of the table at 174.3, so it reads as attached to it.
#   The group spans y 231.4-269.4 (vertical centre 250.4) and starts at x 256.83.
#   Panel c tabulates this same group, so the label uses the same wording:
#     panel c is headed "Array references and cohort tiers", so this label is "Array references".
_LB = "Array references"
_LW = FH.text_length(_LB, 8.0)
page.add_redact_annot(fitz.Rect(177.0, 230.0, 237.0, 241.0))
page.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE)
_LX = GX - 6.0 - _LW                       # 6 pt clear of the left edge of the group
page.insert_text((_LX, 250.4 + 8.0 * 0.5 - 0.299 * 8.0), _LB,
                 fontname="helv", fontsize=8.0, color=INK)
print("라벨 '5S rDNA arrays' : x %.1f-%.1f y 250.4 (표 끝 174.3 에서 %.1f pt 띄움)"
      % (_LX, _LX + _LW, _LX - 174.3))

# -- the zoom indication, from the orange bar down to the blue unit row --
#   Measured: the source bar spans x 256.83-344.08 with its lower edge at y 269.35; the blue row spans x 120.50-514.71 with its top at y 285.0.
#          The existing wedge was narrower than the bar, at x 279.92-328.06, so the source was not legible.
ZX0, ZX1, ZY = 256.83, 344.08, 269.35          # the source bar, the lowest one
AX0, AX1, AY = 120.50, 514.71, 285.00          # the blue unit row
ZCOL = (0.62, 0.68, 0.76)
# A drawn line would cross the labels "Proximal flank" (x 115.7-177.2) and "Distal flank" (x 471-519):
#   measured, the line passes x 170-174 over y 272.5-280.5.
#   So no line is drawn: the source bar is outlined to show what was magnified,
#     and the existing pale wedge over y 271.9-288.6 carries the connection downwards.
_ = (AX0, AX1, AY)          # measured values, used below when the wedge is redrawn
page.draw_rect(fitz.Rect(ZX0 - 0.9, 260.75 - 0.9, ZX1 + 0.9, ZY + 0.9),
               color=ZCOL, width=0.7, radius=0.22)
print("확대 표시 : 막대 x %.1f-%.1f -> 열 x %.1f-%.1f (선 2 + 테두리)" % (ZX0, ZX1, AX0, AX1))

# -- 2) make the zoom region one pale blue shape --
#   As two wedges (x 124.5-279.9 and 328.1-509.4, y 271.9-288.6) it reads
#   as two triangles. It shows one magnified array, so it must be continuous.
#   A single trapezium is filled in the same colour as the wedge, from the source bar
#     (x 256.8-344.1, lower edge 269.4) to the blue row (x 120.5-514.7, top 285.0).
ZFILL = (0.949, 0.961, 0.976)            # the same colour as the existing wedge (measured)
_zq = [fitz.Point(256.83, 269.9), fitz.Point(344.08, 269.9),
       fitz.Point(514.71, 288.6), fitz.Point(120.50, 288.6)]   # 288.6 is the lower edge of the old wedge; 284.6 would leave a white gap.
_sh = page.new_shape()
_sh.draw_polyline(_zq + [_zq[0]])
_sh.finish(color=None, fill=ZFILL, closePath=True)
_sh.commit(overlay=False)                 # drawn beneath the text and bars
print("확대 영역 : 사다리꼴 하나로 채움 (막대 256.8-344.1 -> 열 120.5-514.7)")

# -- 3) lengthen the reads that reach into the flank --
#   Measured structure: grey on the flank side (9.6 pt), purple on the array side, and white ticks.
#   The purple runs 18-48 pt while the reads inside the array (green) run 57-147 pt, so it looks too short.
PUR, GRY = (0.44, 0.34, 0.57), (0.66, 0.66, 0.66)
#   The length is made comparable to the green reads (57-147 pt); any longer and it intrudes on the labels and the green, as measured.
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
# The unit ticks (white lines) are redrawn at the width 2,240 bp occupies in the blue row.
#   The ticks stopped short on the lengthened part.
#   The earlier rule selected the interval by comparing the old and new extents, which stopped three rows of the left group at 174.1.
#   Instead the ticks are drawn across the interval each read actually occupies (left 125.3 to _new, right _new to 509.4),
#     end to end, with the tick positions held on the array grid at 120.50 + k * _UW.
_UW = (514.71 - 120.50) / 125.0           # the CHM13 array holds 125 units (measured)
for _y, _old, _new in LEFTR + RIGHTR:
    _a, _b = (125.3, _new) if _y < 325 else (_new, 509.4)   # the drawn interval, unchanged
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
# The label "Reads reaching into the flank" sits at the same height as the read band (y 314-335) and is lost in it.
#   It is raised into the empty middle of the same line as "Partial unit" (y 298-311).
_LAB = "Reads reaching into the flank"
_lw = page.search_for(_LAB)
if _lw:
    _r = _lw[0]
    page.add_redact_annot(fitz.Rect(_r.x0 - 1, _r.y0 - 1, _r.x1 + 1, _r.y1 + 1))
    # Without an explicit graphics argument the default also erases any artwork the box touches.
    #   This label lay in the middle of the read band (y 314-335), so the white unit ticks drawn
    #   just above touched its box and were erased with it, which is why the second and third rows
    #   on the left were bare beyond x 174.1: the log recorded 21 and 26 ticks drawn
    #   while only 17 survived in the output.
    #   So only the text is deleted.
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

# -- scale the schematic to 423.5 pt --
#   It spans 85.7-519.0 = 433.3 pt while panels c to f are already 423.5, so only the schematic is reduced.
#   This must precede the type-size normalisation; afterwards 4.7 pt would come out as 4.59.
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

# -- normalise the type sizes in the schematic (panels a and b) --
#   tick labels 4.7, axis titles 5.1, legend and annotation 4.5, panel letters 8.2 capitals
#   Only y < 375: below that the panels come from matplotlib and are already to specification.
#   The baseline is taken from the original size; measuring it at the new size lifts the text.
# Selecting by a list of sizes misses glyphs: page.get_text("dict") reports 8.495
#   where page.get_texttrace() reports 8.02 for the same glyph, because they read the transform differently.
#   So the selection is by range, not by a list.
#@DROP
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
    _base = _r.y1 - DESC * _k                 # the baseline measured at the original size
    if _nz == 8.2:
        _x = _r.x0                            # panel letters: the left edge is preserved
    else:
        _x = (_r.x0 + _r.x1) / 2.0 - _F.text_length(_new, _nz) / 2.0
    page.insert_text((_x, _base), _new, fontname="hebo" if _bold else "helv",
                     fontsize=_nz, color=rgb(_c))
print("모식도 글자 규격화 : %d 개" % len(_todo))

# -- panels c and d: haplotype level only, on one row --
#   Produced by mk_Fig1_panelCD.py, which must be run first.
_PT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_Fig1_panelC_TABLE.pdf")
if os.path.exists(_PT):
    _st = fitz.open(_PT)
    _pwt, _pht = _st[0].rect.width, _st[0].rect.height
    # Panels d, e and f are raised: the lower edge of the schematic is measured and they are placed directly beneath it.
    #   A fixed 378.0 would leave a larger gap once the schematic is reduced.
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
    _X0 = 85.9                      # centre the 423.5 pt content width on A4: (595.3 - 423.5) / 2
    _TOP = _CTOP                    # below panel c, the table
    page.show_pdf_page(fitz.Rect(_X0, _TOP, _X0 + _pw, _TOP + _ph), _src, 0)
    _src.close()
    print("패널 c+d 넣음 : y %.0f-%.0f (%.1f x %.1f pt)" % (_TOP, _TOP + _ph, _pw, _ph))
else:
    print("패널 c+d 파일 없음 — mk_Fig1_panelCD.py 를 먼저 돌리십시오")

# The former panel e, a box-and-arrow drawing, was replaced by panel c,
#   a table produced by mk_Fig1_panelC_TABLE.py.
#   mk_Fig1_panelE_COHORT.py and _Fig1_panelE.pdf are kept rather than deleted.

# Figure 1 is not rescaled.
#   It was laid out by hand in Illustrator to the journal's width.
#   The 180 mm page conversion applied to the other figures must not be applied here.
#@DROP

# -- crop to the content --
#   Every other figure is 423.5 x N while Figure 1 alone is A4, so side by side
#           Figure 1 looks smaller, with wide empty bands above and below.
#   set_cropbox is used: it changes only the visible box, redrawing nothing and rescaling nothing,
#     so the width cannot shift. bbox_inches="tight" recomputes the margins and the width does shift.
# Besides, this figure is not a matplotlib output.
#   x0 is measured rather than hard-coded at 85.9, and the width falls out as x0 + 423.5.
_CMAR = 2.0                      # vertical breathing space
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

# A viewer holding the file open blocks the overwrite, so it is written to a temporary file and swapped in.
_tmp = OUT + ".tmp"
doc.save(_tmp, garbage=3, deflate=True)
doc.close()
os.replace(_tmp, OUT)
print("저장 :", OUT)
