# -*- coding: utf-8 -*-
# Rebuilds the two panels of S5, for which there is no drawing code, from the original in one pass
#
# The svg is taken as canonical. The svg coordinates equal the pdf coordinates,
# both 1:1 and from the top left (the centre y of the 6 loci-fraction labels is 109.6, 222.3 and so on, with no digit differing).
# So the whole text layer of the pdf, between BT and ET, is stripped out and rewritten at the positions the svg gives.
#
# What it does:
# 1) cut out axes_33, the long plot at the bottom right, with the tags balanced.
# 2) remove the text inside that plot: "high mode ...", "low mode ...", and that plot's axis names.
# 3) keep only the first and last tick of each axis, chosen by position rather than from a list of values, the ranges differing between axes.
# 4) set the remaining text to T px, which is 3.5 pt on paper.
# 5) panel A: move the row labels ("spacer-1 + ...") from the left to the right
# and shift the whole panel by the empty left margin, without which the frame fitting drags it off the page.
import html
import io
import os
import re
import shutil

import pymupdf


def _need(k):
    """The deposited copy keeps no default: it stops rather than read the wrong place."""
    v = os.environ.get(k)
    if not v:
        raise SystemExit("environment variable %s is not set (the "
                         "deposited copy has no default; see README.md,"
                         " section 'Environment variables')" % k)
    return v


# The default base-14 font, Helvetica, has no beta glyph and would silently render "beta (hardcall)" without it.
# So a real Arial is embedded and used. The widths must also be measured in Arial for the positions to match, the svg being in Arial.
AR = os.path.join("C:" + os.sep, "Windows", "Fonts", "arial.ttf")
ARB = os.path.join("C:" + os.sep, "Windows", "Fonts", "arialbd.ttf")
FONT = {"ari": pymupdf.Font(fontfile=AR), "arib": pymupdf.Font(fontfile=ARB)}

SP = os.path.dirname(os.path.abspath(__file__))
# The four original panels are in the repository, under figS5/inputs_panels.
# search order : the S42_SRC environment value, then inside the repository, then the original.
_R = os.path.dirname(os.path.dirname(SP))
SRC = (os.environ.get("S42_SRC")
       or os.path.join(_R, "figS5", "inputs_panels"))
if not os.path.isdir(SRC):
    SRC = os.path.join(_need("RDNA_ARRAY_METH_DIR"), "03_figures",
                       "00_final_set_2026-10-16", "panels")
# The output folder is taken from PANEL_OUT in the environment.
DST = os.path.join(SP, os.environ.get("PANEL_OUT", "panelsplit_spec"), "_panels")
NL = chr(10)
os.makedirs(DST, exist_ok=True)
D1 = "FIG_D1_CODING32__panel_a"
D5 = "FIG_D5_a4_persite_tone_hardcall32"
T = float(os.environ.get("S42_T", "10.46"))      # the type size inside the panel, in px
TEXT = re.compile(r"<text ([^>]*)>(.*?)</text>", re.S)
GRE = re.compile(r"<g\b|</g>")


def _rgb(h):
    """`#rrggbb` -> (r, g, b) 0~1"""
    return tuple(int(h[1 + 2 * i:3 + 2 * i], 16) / 255.0 for i in range(3))


def cut_block(s, tag):
    """cut out only the group `<g id="tag">`, counting matching tags"""
    i = s.find('<g id="%s">' % tag)
    if i < 0:
        raise SystemExit("not found: " + tag)
    j, depth = i, 0
    while True:
        m = GRE.search(s, j)
        if not m:
            raise SystemExit("pairs do not match: " + tag)
        depth += 1 if m.group(0) == "<g" else -1
        j = m.end()
        if depth == 0:
            break
    return s[:i] + s[j:], j - i


def attr(a, name, default=None):
    m = re.search(r'%s="([^"]*)"' % name, a)
    return m.group(1) if m else default


def style(a, key, default=None):
    m = re.search(r"%s:\s*([^;\"]+)" % key, a)
    return m.group(1).strip() if m else default


def texts(s):
    out = []
    for m in TEXT.finditer(s):
        a, t = m.group(1), html.unescape(m.group(2)).strip()
        if not t:
            continue
        out.append(dict(span=m.span(), raw=m.group(0), x=float(attr(a, "x", "0")),
                        y=float(attr(a, "y", "0")),
                        size=float(style(a, "font-size").replace("px", "")),
                        anchor=style(a, "text-anchor", "start"),
                        bold="700" in (style(a, "font-weight", "") or ""),
                        rot="rotate(-90" in (attr(a, "transform", "") or ""),
                        text=t))
    return out


# text replacements
# a) the x-axis name of B: "beta (hardcall)" becomes "Unit methylation beta", the hard call being a matter for the legend
# b) a middle dot becomes a comma
WORD = ((chr(0x3b2) + " (hardcall)", "Unit methylation " + chr(0x3b2)),
        (" " + chr(0xb7) + " ", ", "))


def _word(t):
    for a_, b_ in WORD:
        t = t.replace(a_, b_)
    return t


def _resize(m):
    r = re.sub(r"font-size: [\d.]+px", "font-size: %gpx" % T, m.group(0), count=1)
    return r[:r.index(">") + 1] + _word(r[r.index(">") + 1:r.rindex("<")])         + r[r.rindex("<"):]


# ──────────────────────────────── panel B (D5)
s = io.open(os.path.join(SRC, D5 + ".svg"), encoding="utf8").read()
n_before = s.count("<text")
s, ncut = cut_block(s, "axes_33")
print("D5: axes_33 %d characters cut out" % ncut)

keep = set()
ax = [m.start() for m in re.finditer(r'<g id="axes_\d+">', s)] + [len(s)]
for a0, a1 in zip(ax, ax[1:]):
    for kind in ("xtick", "ytick"):
        pos = [m.start() for m in re.finditer(r'<g id="%s_\d+">' % kind, s[a0:a1])]
        if not pos:
            continue
        for p in (pos[0], pos[-1]):
            tm = TEXT.search(s, a0 + p)
            if tm and tm.start() < a1:
                keep.add(tm.span())

drop = []
for t in texts(s):
    if t["text"].startswith(("high mode", "low mode")):
        drop.append(t)
    elif abs(t["size"] - 7.5) < 0.1:          # the axis names of the removed plot
        drop.append(t)
    elif abs(t["size"] - 6.5) < 0.2 and t["span"] not in keep:
        drop.append(t)
out = s
for t in sorted(drop, key=lambda z: -z["span"][0]):
    out = out[:t["span"][0]] + out[t["span"][1]:]
out = TEXT.sub(_resize, out)
io.open(os.path.join(DST, D5 + ".svg"), "w", encoding="utf8", newline=NL).write(out)
print("D5: texts %d -> %d (thinned %d), size %.2f px" %
      (n_before, out.count("<text"), len(drop), T))

# ──────────────────────────────── panel A (D1)
s1 = io.open(os.path.join(SRC, D1 + ".svg"), encoding="utf8").read()
W1 = float(re.search(r'<svg[^>]*width="([\d.]+)pt"', s1).group(1))


def poly(m):
    xs, ys = [float(m.group(1))], [float(m.group(2))]
    for a, b in re.findall(r"L\s*([-\d.]+)\s+([-\d.]+)", m.group(3)):
        xs.append(float(a))
        ys.append(float(b))
    return xs, ys


# the x window is cut at the end of the unit
# The bars end at 2240 bp, x 853.9, while the axis frame runs to 2600 bp, x 959.2.
# A window wider than the unit makes a region that does not exist look as though it does, so the frame, the ticks and the axis name are all moved.
i1 = s1.find('<g id="axes_1">')
fx, fy = poly(re.search(r'<path d="M\s*([-\d.]+)\s+([-\d.]+)\s*(L[^"]*?)"',
                        s1[i1:i1 + 4000], re.S))
XR0, FY0, FY1 = max(fx), min(fy), max(fy)
XL0 = min(fx)
bars = []
for m in re.finditer(r'<path d="M\s*([-\d.]+)\s+([-\d.]+)\s*(L[^"]*?)z', s1, re.S):
    bx, _by = poly(m)
    if max(bx) < XR0 - 1:
        bars.append(max(bx))
XCUT = max(bars)
print("D1: x window %.1f -> %.1f (data %.0f -> %.0f bp)"
      % (XR0, XCUT, (XR0 - XL0) / 0.29340, (XCUT - XL0) / 0.29340))
# a) pull in the right line of the frame and the background
blk = s1[i1:]
blk = re.sub(r"(?<![\d.])%s(?![\d.])" % re.escape("%g" % XR0), "%.4f" % XCUT, blk)
s1 = s1[:i1] + blk
# b) remove the tick groups outside the cut, a tick line and its number forming one group
while True:
    gone = None
    for m in re.finditer(r'<g id="(xtick_\d+)">', s1):
        tm = TEXT.search(s1, m.start())
        if tm and float(attr(tm.group(1), "x", "0")) > XCUT + 0.5:
            gone = m.group(1)
            break
    if not gone:
        break
    s1, _c = cut_block(s1, gone)
    print("D1: ticks outside the window %s dropped" % gone)
# c) remove the x-axis name, "array unit coordinate (v11 canonical; ...)", and
# the panel title, "a) which part of the array unit ...".
# (v11 is a release name.)
# Instead, which position belongs to which segment is shown by a colour band below the axis, as in Figure 4F; see e) below.
# Deleted from the back forwards: deleting from the front shifts the later offsets and cuts the wrong place.
for _t in sorted([t for t in texts(s1)
                  if "coordinate" in t["text"] or t["text"].startswith("a) which part")],
                 key=lambda z: -z["span"][0]):
    s1 = s1[:_t["span"][0]] + s1[_t["span"][1]:]
    print("D1: dropped text %r" % _t["text"][:46])

# d) the bar colours are set to the paper's segment palette, matching the panels above and below
# (the source panel has its own orange-to-blue palette; the colours of Figures 4E, 4F and 5 replace it).
BARCOL = [("#f6d5a8", "#9ecae1", "spacer-1"), ("#cfe3d4", "#a6dbb4", "GA-rich"),
          ("#bcd6e8", "#d9d9d9", "spacer-2"), ("#9ec5e0", "#f2b705", "5'-flank"),
          ("#2a6faf", "#9e9e9e", "5S gene")]
for _o, _n, _nm in BARCOL:
    s1 = s1.replace("fill: %s" % _o, "fill: %s" % _n)

# e) the segment colour band and its names, by the same rule as Figure 4F: a wide segment is labelled in line, a narrow one on a second row with a leader
ELN = [("AluY", 1, 318), ("spacer-1", 319, 1178), ("GA-rich", 1179, 1551),
       ("spacer-2", 1552, 2031), ("5'-flank", 2032, 2121), ("5S gene", 2122, 2240)]
ELC = {"AluY": "#c7b2de", "spacer-1": "#9ecae1", "GA-rich": "#a6dbb4",
       "spacer-2": "#d9d9d9", "5'-flank": "#f2b705", "5S gene": "#9e9e9e"}
BPX = (XCUT - XL0) / 2240.0                    # how many pt one bp is
BY0, BH = FY1 + 19.0, 9.0                      # the top y of the band, and its height
LINE1, LINE2 = BY0 + BH + T * 1.05, BY0 + BH + T * 2.35
BAND, BTX = [], []
_last = {LINE1: -1e9, LINE2: -1e9}
for _nm, _a, _b in ELN:
    x0, x1 = XL0 + (_a - 1) * BPX, XL0 + _b * BPX
    BAND.append((x0, BY0, x1, BY0 + BH, ELC[_nm]))
    _w = FONT["ari"].text_length(_nm, T)
    _cx = (x0 + x1) / 2.0
    # where a name is wider than its band it collides with its neighbour, so the names alternate between two rows with a leader drawn
    _y = LINE1 if _cx - _w / 2.0 > _last[LINE1] + T * 0.4 else LINE2
    if _cx - _w / 2.0 <= _last[_y] + T * 0.4:
        _y = LINE2 if _y == LINE1 else LINE1
    _last[_y] = _cx + _w / 2.0
    BTX.append((_nm, _cx, _y, _y != LINE1))
NEWH = LINE2 + T * 0.45
tt = texts(s1)
ann = [t for t in tt if "sites" in t["text"] and "loci" in t["text"]]
lab = [t for t in tt if t["anchor"] == "end"]
assert len(ann) == 5 and len(lab) == 5, (len(ann), len(lab))
annw = max(FONT["ari"].text_length(t["text"], T) for t in ann)
NEWX = ann[0]["x"] + annw + T * 1.2
new1 = s1
for t in sorted(lab, key=lambda z: -z["span"][0]):
    r = t["raw"]
    r = re.sub(r"text-anchor: \w+", "text-anchor: start", r)
    r = re.sub(r'x="[-\d.]+"', 'x="%.3f"' % NEWX, r, count=1)
    r = re.sub(r"rotate\(-0 [-\d.]+ ", "rotate(-0 %.3f " % NEWX, r)
    new1 = new1[:t["span"][0]] + r + new1[t["span"][1]:]
new1 = TEXT.sub(_resize, new1)
RIGHT = NEWX + max(FONT["ari"].text_length(t["text"], T) for t in lab) + 2.0
XMIN = min(t["x"] for t in texts(new1) if t["anchor"] == "middle")
DX = -max(0.0, XMIN - T * 1.6)
NEWW = RIGHT + DX + 2.0
# the segment band and its names are written into the svg; the text is read back from this svg when the pdf text layer is rebuilt
_band = []
for x0, y0, x1, y1, c in BAND:
    _band.append('   <rect x="%.3f" y="%.3f" width="%.3f" height="%.3f" '
                 'style="fill: %s; stroke: none"/>' % (x0, y0, x1 - x0, y1 - y0, c))
for _nm, _cx, _y, _low in BTX:
    if _low:      # a name on the lower row is given a leader up to its band
        _band.append('   <path d="M %.3f %.3f L %.3f %.3f" style="stroke: #9aa3ad; '
                     'stroke-width: 0.5; fill: none"/>'
                     % (_cx, BY0 + BH, _cx, _y - T * 0.95))
    _band.append('   <text style="font-size: %gpx; font-family: \'Arial\', '
                 '\'DejaVu Sans\', sans-serif; text-anchor: middle" x="%.3f" y="%.3f" '
                 'transform="rotate(-0 %.3f %.3f)">%s</text>'
                 % (T, _cx, _y, _cx, _y, html.escape(_nm)))
# the band is attached outside figure_1, so the whole-panel translate is applied directly in its own transform
new1 = new1.replace("</svg>", NL.join(
    ['  <g id="element_band" transform="translate(%.3f,0)">' % DX] + _band + ["  </g>"])
    + NL + "</svg>", 1)
new1 = re.sub(r'(<svg[^>]*width=")[\d.]+(pt" height=")[\d.]+(pt")',
              r"\g<1>%.2f\g<2>%.2f\g<3>" % (NEWW, NEWH), new1, count=1)
new1 = re.sub(r'viewBox="0 0 [\d.]+ [\d.]+"',
              'viewBox="0 0 %.2f %.2f"' % (NEWW, NEWH), new1, count=1)
new1 = new1.replace('<g id="figure_1">',
                    '<g id="figure_1" transform="translate(%.3f,0)">' % DX, 1)
io.open(os.path.join(DST, D1 + ".svg"), "w", encoding="utf8", newline=NL).write(new1)
print("D1: row label x %.1f -> %.1f, whole shifted by %.1f and width %.1f -> %.1f"
      % (lab[0]["x"], NEWX, DX, W1, NEWW))


# ──────────────────────────────── the pdf: the text layer is rebuilt
# the lines and points of axes_33, the long plot at the bottom right, are not text and so are not stripped out
# and that area is therefore erased. The frame measures x 596.3-1104.8 and y 638.1-707.6, plus room for the ticks.
CUT33 = (585.0, 625.0, 1116.0, 748.8)


def rebuild_pdf(name, dx=0.0, neww=None, newh=None, cut=None, box=None,
                band=None, leads=(), recolor=None):
    shutil.copy2(os.path.join(SRC, name + ".pdf"), os.path.join(DST, name + ".pdf"))
    q = os.path.join(DST, name + ".pdf")
    d = pymupdf.open(q)
    pg = d[0]
    oldh = pg.rect.height
    # When the page is extended downwards the content must be raised by the same amount:
    # the pdf origin is at the bottom left, so without raising it the new space appears at the top.
    dy = (newh - oldh) if newh else 0.0
    pg.wrap_contents()
    cont = pg.read_contents()
    n0 = len(re.findall(rb"BT", cont))
    cont = re.sub(rb"BT.*?ET", b"", cont, flags=re.S)       # the text layer is stripped out
    # changing the bar colours requires editing the r g b rg operators of the pdf content stream as well as the svg
    # (the png is rendered from the pdf, so fixing only the svg leaves the two out of step)
    for _o, _n, _nm in (recolor or ()):
        _a = (" ".join("%.10g" % (int(_o[1 + 2 * i:3 + 2 * i], 16) / 255.0)
                       for i in range(3))).encode()
        _b = (" ".join("%.10g" % (int(_n[1 + 2 * i:3 + 2 * i], 16) / 255.0)
                       for i in range(3))).encode()
        if _a + b" rg" in cont:
            cont = cont.replace(_a + b" rg", _b + b" rg")
        else:
            print(" colour %s(%s) not found in the pdf" % (_nm, _o))
    cont = ("q 1 0 0 1 %.4f %.4f cm" % (dx, dy)).encode() + NL.encode() + cont \
        + (NL + "Q").encode()
    xref = pg.get_contents()[0]
    d.update_stream(xref, cont)
    pg.set_contents(xref)
    if neww or newh:
        pg.set_mediabox(pymupdf.Rect(0, 0, neww or pg.rect.width, newh or oldh))
        pg = d[0]
    if band:
        for x0, y0, x1, y1, c in band:
            pg.draw_rect(pymupdf.Rect(x0 + dx, y0, x1 + dx, y1),
                         color=None, fill=_rgb(c), width=0)
        print("   segment color band %d cells" % len(band))
    for _cx, _y0, _y1 in leads:
        pg.draw_line(pymupdf.Point(_cx + dx, _y0), pymupdf.Point(_cx + dx, _y1),
                     color=(0.604, 0.639, 0.694), width=0.5)
    if cut:
        # Redaction does not cut a line: IF_COVERED, the default, erases only what lies entirely inside.
        # That is enough for the long plot (D5), but the top and bottom frame lines of A only cross the box and would survive.
        # So A uses IF_TOUCHED, which erases a crossing line too, and the frame is redrawn as a box.
        pg.add_redact_annot(pymupdf.Rect(*cut))
        pg.apply_redactions(images=pymupdf.PDF_REDACT_IMAGE_NONE,
                            graphics=(pymupdf.PDF_REDACT_LINE_ART_REMOVE_IF_TOUCHED
                                      if box else
                                      pymupdf.PDF_REDACT_LINE_ART_REMOVE_IF_COVERED))
        print("   erased positions %s" % (cut,))
    if box:
        pg.draw_rect(pymupdf.Rect(*box), color=(0, 0, 0), width=0.8)
        print("   frame redrawn %s" % (box,))
    sv = io.open(os.path.join(DST, name + ".svg"), encoding="utf8").read()
    pg.insert_font(fontname="ari", fontfile=AR)
    pg.insert_font(fontname="arib", fontfile=ARB)
    put = 0
    for t in texts(sv):
        fn = "arib" if t["bold"] else "ari"
        w = FONT[fn].text_length(t["text"], t["size"])
        if t["rot"]:
            pg.insert_text(pymupdf.Point(t["x"] + dx, t["y"] + w / 2.0), t["text"],
                           fontname=fn, fontsize=t["size"], rotate=90)
        else:
            x = t["x"] - (w if t["anchor"] == "end" else
                          w / 2.0 if t["anchor"] == "middle" else 0.0)
            pg.insert_text(pymupdf.Point(x + dx, t["y"]), t["text"], fontname=fn,
                           fontsize=t["size"])
        put += 1
    d.save(os.path.join(DST, name + ".tmp.pdf"), deflate=True, garbage=3)
    d.close()
    os.replace(os.path.join(DST, name + ".tmp.pdf"), q)
    print("%-36s pdf: %d blocks of the existing text layer stripped, %d written again" % (name, n0, put))


rebuild_pdf(D5, cut=CUT33)
# A: the frame and tick artwork outside the window, to the right of 2240 bp, is erased and a new right line drawn.
# The coordinates are given after the whole-panel translate, redaction and drawing happening after it.
rebuild_pdf(D1, dx=DX, neww=NEWW, newh=NEWH,
            cut=(XCUT + DX + 0.6, FY0 - 6.0, XR0 + DX + 3.0, FY1 + 6.0),
            box=(XL0 + DX, FY0, XCUT + DX, FY1), band=BAND, recolor=BARCOL,
            leads=[(cx, BY0 + BH, _y - T * 0.95) for _n, cx, _y, low in BTX if low])
