# -*- coding: utf-8 -*-
# Figure 1A: the unit in genome orientation (colour only), an X crossover, then the unit in the frame used here (with labels and ticks),
#            with the direction of transcription marked in all three (a bent arrow at the TSS) and position 1 at the left of chr1.
#
# (the final, hand-drawn version)
#     The upper row shows the flipped orientation with the segments marked by colour alone,
#       the row beneath shows the unit in the frame used here, and the direction of transcription is marked in every case.
#     Three layers: the array, then the magnified unit in genome orientation, the X, and the unit in this frame.
#     The X reads clearly because the boxes above and below it are the same width.
#        (Earlier the upper vertex was only 6.4 pt wide and the X was invisible.)
#
#   Why the two orientations are opposite, measured:
#     - blastn over CHM13 chr1:227,740,000-228,030,000: 126 of the 127 5S genes lie on the minus strand
#       (neighbouring genes are 2,237-2,240 bp apart, matching the 2,240 bp frame used here).
#     - Flank windows: proximal chr1:227,704,289-227,744,288 (low)
#                 and distal chr1:228,026,263-228,066,262 (high).
#       So transcription runs towards decreasing chr1 coordinates, distal to proximal, right to left in the figure.
#     - This frame: AluY 1-318 ... 5'-flanking 2032-2121, 5S gene 2122-2240 (TSS at 2122).
#       Increasing frame coordinate is the direction of transcription, so drawing the frame 1 to 2240 left to right opposes the genome.
#     - The manuscript states that right to left is the direction of transcription, and
#             that the proximal end of the array is drawn on the left.
#     So in genome orientation the gene is on the left, AluY on the right, and transcription runs leftwards;
#        in the frame used here AluY is on the left, the gene on the right, and transcription runs rightwards.
#
#   The input is the current output, not the Illustrator original: fig1.pdf (md5 150ddf302139)
#     is absent from the data and work roots to a depth of 6, so the main script cannot be re-run.
#     This script therefore patches the output; once the original is available it should be folded into the main script.
import os
import sys
import fitz

sys.stdout.reconfigure(encoding="utf-8")

SRC = os.environ["FIG1_CURRENT_PDF"]
OUT = os.environ["FIG1_MARKED_PDF"]

# -- measured positions, in the crop coordinates of the current PDF, page 423.5 x 559.84 --
CHR_X0 = 37.00
BAR_Y0, BAR_Y1 = 43.74, 52.00        # the array bar
SRC_X0, SRC_X1 = 228.52, 234.92      # the cell being magnified
BOX_X0, BOX_X1 = 138.33, 346.09      # left and right edges of the unit block, shared by both layers
ZFILL = (0.949, 0.961, 0.976)
EDGE = (0.29, 0.478, 0.69)
INK = (0.102, 0.102, 0.102)
GREY = (0.361, 0.380, 0.412)
LEAD = (0.690, 0.416, 0.416)
XLINE = (0.302, 0.349, 0.447)         # the X crossover lines, darkened from 0.55 / 0.60 / 0.68
FS_EL, FS_TICK, FS_RU = 4.70, 4.70, 5.10
UTOT = 2240.0
MIR = BOX_X0 + BOX_X1

# the vertical layout, dividing the space between the array at 52.0 and the "B" of panel B at 91.69
FAN_Y0, FAN_Y1 = 56.60, 62.00        # the zoom lines, leaving room for the arrow below the array
GEN_Y0, GEN_Y1 = 62.00, 68.30        # the unit in genome orientation, colour only, 6.30 pt
XB_Y0, XB_Y1 = 72.30, 77.70          # the X band, starting below the arrow of the genome-orientation unit (68.3 -> 70.9)
XB_IN = 8.0                           # the ends of the X are inset from the corners so they do not overlap the left arrow
OUR_Y0, OUR_Y1 = 77.70, 84.00        # the unit in this frame, with labels, 6.30 pt to match the layer above
TICK_BASE = 91.50                    # the tick baseline; panel B's ink begins at y 95.4 (measured)

ELEM = [("AluY", 1, 318, (0.725, 0.804, 0.894)),
        ("spacer-1", 319, 1178, (0.890, 0.925, 0.965)),
        ("GA-rich", 1179, 1551, (0.953, 0.875, 0.722)),
        ("spacer-2", 1552, 2031, (0.890, 0.925, 0.965)),
        ("5'-flanking", 2032, 2121, (0.835, 0.878, 0.812)),
        ("5S rDNA gene", 2122, 2240, (0.918, 0.722, 0.682))]
INSIDE = {"AluY", "spacer-1", "GA-rich", "spacer-2"}
TICKS = [1, 319, 1179, 1552, 2240]


def xf(bp):
    """Frame coordinate -> x in the unit of our frame (1 is at the left)"""
    return BOX_X0 + (BOX_X1 - BOX_X0) * (bp - 1) / UTOT


def xm(bp):
    """Frame coordinate -> x in the genome-orientation unit (mirrored)"""
    return MIR - xf(bp)


def tlen(t, sz, bold=False):
    return fitz.get_text_length(t, fontname="hebo" if bold else "helv", fontsize=sz)


def bent(page, x_tss, y_base, rise, run, w=0.85, head=2.30, hw=1.10, color=INK):
    """The bent TSS arrow, vertical first and then along transcription (to the left when run<0).
    The head is a filled triangle, because a V of two strokes breaks up at this size."""
    yt = y_base - rise
    s = 1.0 if run > 0 else -1.0
    xe = x_tss + run                 # the tip of the arrowhead
    xs = xe - s * head               # where the shaft ends and the head begins
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
    """One row of units. xpos(bp)->x, and labels=True writes the element names."""
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
    # the two narrow segments are labelled outside with a leader
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
print("input %s (%.1f x %.1f pt)" % (os.path.basename(SRC), page.rect.width, page.rect.height))

# -- 1) delete the lower half of the old panel A: the zoom lines, unit, ticks, outside labels and "repeat unit" --
for r in (fitz.Rect(60.0, 52.30, 410.0, 90.0),
          fitz.Rect(64.0, 69.5, 120.0, 78.5)):
    page.draw_rect(r, color=None, fill=(1, 1, 1))

# -- 2) the zoom lines, from one cell of the array to the unit in genome orientation --
fan = page.new_shape()
fan.draw_polyline([fitz.Point(SRC_X0, FAN_Y0), fitz.Point(SRC_X1, FAN_Y0),
                   fitz.Point(BOX_X1, FAN_Y1), fitz.Point(BOX_X0, FAN_Y1),
                   fitz.Point(SRC_X0, FAN_Y0)])
fan.finish(color=None, fill=ZFILL)
fan.commit()

# -- 3) the unit in genome orientation, colour only, with no labels or ticks --
unit_row(page, GEN_Y0, GEN_Y1, xm, labels=False, outside_left=True)

# -- 4) the X crossover; the boxes above and below are the same width, so the crossing reads clearly --
for (xa, xb) in ((BOX_X0 + XB_IN, BOX_X1 - XB_IN), (BOX_X1 - XB_IN, BOX_X0 + XB_IN)):
    page.draw_line(fitz.Point(xa, XB_Y0), fitz.Point(xb, XB_Y1),
                   color=XLINE, width=0.55, dashes="[1.8 1.5] 0")   # dashed, and thinner

# -- 5) the unit in this frame, with labels and ticks --
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

# -- 6) the direction of transcription, in all three places --
#   In the array, frame position 1 is at the right edge of a cell, because transcription runs towards the proximal end on the left.
TSS_BAR = SRC_X1 - (2121.0 / UTOT) * (SRC_X1 - SRC_X0)
bent(page, TSS_BAR, BAR_Y1, -2.30, -(SRC_X1 - SRC_X0))  # below the bar, 6.40 pt wide
page.insert_text((TSS_BAR + 2.2, 56.00), "transcription",
                 fontname="helv", fontsize=FS_EL, color=INK)
#   In genome orientation the TSS (2122) is on the left.
bent(page, xm(2122), GEN_Y1, -2.60, -7.0)   # below the box
#   In this frame the TSS (2122) is on the right, because the frame opposes the genome.
bent(page, xf(2122), OUR_Y1, -2.25, 7.0)    # below the box

# -- 7) position 1 is at the left end of chr1 --
page.draw_line(fitz.Point(CHR_X0, 22.4), fitz.Point(CHR_X0, 24.6), color=INK, width=0.7)
page.insert_text((CHR_X0 - tlen("1", 6.0) / 2.0, 30.2), "1",
                 fontname="helv", fontsize=6.0, color=INK)

tmp = OUT + ".tmp"
doc.save(tmp, garbage=3, deflate=True)
doc.close()
os.replace(tmp, OUT)
print("saved :", OUT)
print("genome-orientation unit y %.1f-%.1f, X band %.1f-%.1f, our unit %.1f-%.1f, ticks %.1f"
      % (GEN_Y0, GEN_Y1, XB_Y0, XB_Y1, OUR_Y0, OUR_Y1, TICK_BASE))
