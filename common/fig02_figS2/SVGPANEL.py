# -*- coding: utf-8 -*-
r"""Panel cutting module for the combined Fig 2 (placement, not recomputation).

   - In a published SVG (matplotlib), one `<g id="axes_N">` is one panel.
   - Each panel is located by the background rectangle of its axes, and regions are split at the midpoint to the neighbouring axes
     (so that axis labels and titles fall inside the region).
   - Panel and letter are matched by reading the `<text>` inside that axes, not guessed from position.
   Nothing is redrawn."""
import io, os, re, sys, copy
import xml.etree.ElementTree as ET

SVG = "http://www.w3.org/2000/svg"
XLINK = "http://www.w3.org/1999/xlink"
ET.register_namespace("", SVG)
ET.register_namespace("xlink", XLINK)
NUM = r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?"


def parse(path):
    t = io.open(path, encoding="utf8").read()
    t = re.sub(r"<!DOCTYPE[^>]*>", "", t, flags=re.S)
    return ET.fromstring(t)


def viewbox(root):
    return tuple(float(x) for x in root.get("viewBox").split())


def _axes_rect(g):
    """Axes background rectangle: uses the path of the first child <g id="patch_N">,
       otherwise the widest path (choosing by colour picks a white cell inside a table)."""
    for cp in g.iter("{%s}clipPath" % SVG):       # the most reliable position
        for r in cp.iter("{%s}rect" % SVG):
            try:
                x, y = float(r.get("x")), float(r.get("y"))
                w, h = float(r.get("width")), float(r.get("height"))
            except (TypeError, ValueError):
                continue
            if w > 1 and h > 1:
                return (x, y, x + w, y + h)

    def box_of(path_el):
        d = path_el.get("d") or ""
        if not d.strip().startswith("M"): return None
        nums = [float(x) for x in re.findall(NUM, d)]
        if len(nums) < 8: return None
        xs, ys = nums[0::2][:4], nums[1::2][:4]
        return (min(xs), min(ys), max(xs), max(ys))

    for ch in list(g):
        if (ch.get("id") or "").startswith("patch_"):
            for pe in ch.iter("{%s}path" % SVG):
                b = box_of(pe)
                if b and (b[2] - b[0]) > 1 and (b[3] - b[1]) > 1:
                    return b
            break
    best = None
    for pe in g.iter("{%s}path" % SVG):
        b = box_of(pe)
        if not b: continue
        a = (b[2] - b[0]) * (b[3] - b[1])
        if best is None or a > best[0]: best = (a, b)
    return best[1] if best else None


def _content_bbox(g):
    """bbox of text and shape coordinates inside the group (for restoring), only translate is considered."""
    xs, ys = [], []
    for e in g.iter():
        tag = e.tag.split("}")[-1]
        if tag == "text":
            try:
                x, y = e.get("x"), e.get("y")
                if x is not None and y is not None: xs.append(float(x)); ys.append(float(y))
            except ValueError: pass
            tr = e.get("transform") or ""
            m = re.search(r"translate\(\s*(%s)[\s,]+(%s)" % (NUM, NUM), tr)
            if m: xs.append(float(m.group(1))); ys.append(float(m.group(2)))
        elif tag in ("rect", "use", "image"):
            try:
                if e.get("x") is not None: xs.append(float(e.get("x")))
                if e.get("y") is not None: ys.append(float(e.get("y")))
            except ValueError: pass
        elif tag == "g":
            tr = e.get("transform") or ""
            m = re.search(r"translate\(\s*(%s)[\s,]+(%s)" % (NUM, NUM), tr)
            if m: xs.append(float(m.group(1))); ys.append(float(m.group(2)))
    if len(xs) < 2: return None
    return (min(xs), min(ys), max(xs), max(ys))


def texts(g):
    out = []
    for e in g.iter("{%s}text" % SVG):
        s = "".join(e.itertext()).strip()
        if s: out.append(s)
    return out


def _outside_text_rows(root):
    """y range of the <text> elements outside the axes groups, usually the suptitle."""
    inside = set()
    for g in root.iter("{%s}g" % SVG):
        if (g.get("id") or "").startswith("axes_"):
            for e in g.iter("{%s}text" % SVG): inside.add(id(e))
    ys = []
    for e in root.iter("{%s}text" % SVG):
        if id(e) in inside: continue
        s = "".join(e.itertext()).strip()
        if not s: continue
        y = e.get("y")
        if y is None:
            tr = e.get("transform") or ""
            m = re.search(r"translate\(\s*%s\s+(%s)" % (NUM, NUM), tr)
            if m: y = m.group(1)
        if y is not None:
            try: ys.append(float(y))
            except ValueError: pass
    return ys


DROP_TITLE = True        # remove the panel title line from the region (avoids a duplicate of the composite cell letter)
_LAB = re.compile(r"^[A-Za-z]\d?$")


def _text_y(e):
    y = e.get("y")
    if y is None:
        tr = e.get("transform") or ""
        m = re.search(r"translate\(\s*%s[\s,]+(%s)" % (NUM, NUM), tr)
        if m: y = m.group(1)
    try: return float(y)
    except (TypeError, ValueError): return None


def strip_panel_letter(g, ax_rect, label):
    """Removes the leading letter of an original panel (`A  ...`, `B1  ...`) from the title only.
       - title = <text> above the axes rectangle (data labels are inside)
       - only the letter is removed and the title text is kept, because the composite adds its own letter and there would be two.
       - when the label is not letter-shaped (such as `(single)`), nothing is done.
       Returns: the number of removals"""
    if not label or not _LAB.match(label.strip()): return 0
    lab = label.strip()
    top = ax_rect[1]
    pat = re.compile(r"^%s(?=[\s ])[\s ]*" % re.escape(lab))
    n = 0
    for e in g.iter("{%s}text" % SVG):
        y = _text_y(e)
        if y is None or y >= top: continue          # text inside the axes is not touched
        if not pat.match("".join(e.itertext()).strip()): continue
        for node in [e] + list(e.iter()):           # remove only from the first text node
            if node.text and node.text.strip():
                node.text = pat.sub("", node.text.lstrip(), count=1) or node.text.lstrip()
                n += 1
                break
    return n


INK_EXPAND = False # not used; REGION='mid+ink' is applied instead


def _tw(s, fs):
    """Estimated text width (pt): 1.0 em for Hangul and CJK characters, 0.56 em otherwise"""
    return sum((1.0 if ord(ch) >= 0x1100 else 0.56) * fs for ch in s)


def _ink_x(g, skip_above=None):
    """Horizontal ink range [x0, x1] of a group (text width is estimated). skip_above: text above this is excluded (title)"""
    xs = []
    for e in g.iter():
        tag = e.tag.split("}")[-1]
        if tag == "text":
            s = "".join(e.itertext()).strip()
            if not s: continue
            st = e.get("style") or ""
            m = re.search(r"font-size:\s*(%s)px" % NUM, st)
            fs = float(m.group(1)) if m else 10.0
            tr = e.get("transform") or ""
            x = y = None
            if e.get("x") is not None:
                try: x, y = float(e.get("x")), float(e.get("y"))
                except (TypeError, ValueError): pass
            mt = re.search(r"translate\(\s*(%s)[\s,]+(%s)" % (NUM, NUM), tr)
            if mt: x, y = float(mt.group(1)), float(mt.group(2))
            if x is None: continue
            if skip_above is not None and y is not None and y < skip_above: continue
            mr = re.search(r"rotate\(\s*(%s)" % NUM, tr)
            ang = abs(float(mr.group(1))) if mr else 0.0
            w = _tw(s, fs)
            if 60 < ang % 180 < 120:            # vertical text, horizontally only the text height
                xs += [x - fs, x + 0.3 * fs]; continue
            an = "middle" if "text-anchor: middle" in st else ("end" if "text-anchor: end" in st else "start")
            a = {"start": x, "middle": x - w / 2, "end": x - w}[an]
            xs += [a, a + w]
        elif tag == "path" and e.get("d"):
            v = [float(t) for t in re.findall(NUM, e.get("d"))]
            _d = re.sub(r"\d[eE][-+]?\d", "0", e.get("d"))
            if len(v) >= 2 and not re.search(r"[mlhvcsqta]", _d) and not re.search(r"[HV]", _d):   # absolute coordinate paths only (z close is allowed)
                xs += v[0::2]
        elif tag == "use" and e.get("x") is not None:
            try: xs.append(float(e.get("x")))
            except ValueError: pass
    return (min(xs), max(xs)) if xs else None


REGION = "mid+ink"        # midline region, widened up to EXT only on the side where its own ink spills over ('mid' = old rule)
EXT = 40.0
PAD = 3.0


def _rot(tr):
    m = re.search(r"rotate\(\s*(%s)" % NUM, tr or "")
    return float(m.group(1)) if m else 0.0


def _text_box(e):
    """Estimated ink rectangle of one text element (rotation included), (x0, y0, x1, y1) or None"""
    import math
    s = "".join(e.itertext()).strip()
    if not s: return None
    st = e.get("style") or ""
    m = re.search(r"font-size:\s*(%s)px" % NUM, st)
    fs = float(m.group(1)) if m else 10.0
    tr = e.get("transform") or ""
    x = y = None
    if e.get("x") is not None:
        try: x, y = float(e.get("x")), float(e.get("y"))
        except (TypeError, ValueError): pass
    mt = re.search(r"translate\(\s*(%s)[\s,]+(%s)" % (NUM, NUM), tr)
    if mt: x, y = float(mt.group(1)), float(mt.group(2))
    if x is None: return None
    w = sum((1.0 if ord(ch) >= 0x1100 else 0.55) * fs for ch in s)
    k = 0.5 if "text-anchor: middle" in st else (1.0 if "text-anchor: end" in st else 0.0)
    th = math.radians(_rot(tr))
    dx, dy = math.cos(th), math.sin(th)          # text direction
    nx, ny = -dy, dx                             # perpendicular to it (text height direction)
    sx, sy = x - k * w * dx, y - k * w * dy
    pts = []
    for a in (0.0, w):
        for b in (-0.85 * fs, 0.30 * fs):        # ascender and descender
            pts.append((sx + a * dx + b * nx, sy + a * dy + b * ny))
    xs, ys = [p[0] for p in pts], [p[1] for p in pts]
    return (min(xs), min(ys), max(xs), max(ys))


def _is_title(e, axr):
    y = _text_y(e)
    return (y is not None and y < axr[1] - 1
            and not (60 < abs(_rot(e.get("transform"))) % 180 < 120))


def _ink_box(g, axr, skip_title=False):
    """Ink rectangle of the whole group: text width is estimated, and shapes with a clip-path are bounded to the axes rectangle"""
    pm = {c: p for p in g.iter() for c in p}
    def in_defs(e):
        while e is not None and e is not g:
            if e.tag.split("}")[-1] in ("defs", "clipPath", "marker"): return True
            e = pm.get(e)
        return False
    def clipped(e):
        while e is not None:
            if e.get("clip-path"): return True
            if e is g: return False
            e = pm.get(e)
        return False
    X0, Y0, X1, Y1 = axr
    xs, ys = [], []
    for e in g.iter():
        tag = e.tag.split("}")[-1]
        if in_defs(e): continue                                    # marker and tick shapes defined at the origin
        if tag == "text":
            if skip_title and _is_title(e, axr): continue
            b = _text_box(e)
            if b: xs += [b[0], b[2]]; ys += [b[1], b[3]]
        elif tag in ("path", "use"):
            if tag == "path" and e.get("d"):
                d = re.sub(r"\d[eE][-+]?\d", "0", e.get("d"))
                if re.search(r"[mlhvcsqtaHV]", d): continue          # relative coordinates and single-axis moves are skipped
                v = [float(t) for t in re.findall(NUM, e.get("d"))]
                px, py = v[0::2], v[1::2]
            else:
                try: px, py = [float(e.get("x"))], [float(e.get("y"))]
                except (TypeError, ValueError): continue
            if not px or not py: continue
            if clipped(e):
                px = [min(max(p, X0), X1) for p in px]; py = [min(max(p, Y0), Y1) for p in py]
            xs += px; ys += py
    if not xs: return None
    return (min(xs), min(ys), max(xs), max(ys))


def _drop_title_texts(g, axr):
    """Deletes title text (above the axes, not rotated) from the group. The region is not cut, so vertical axis titles survive"""
    pm = {c: p for p in g.iter() for c in p}
    top, n = axr[1], 0
    for e in list(g.iter("{%s}text" % SVG)):
        if not _is_title(e, axr): continue          # vertical text (vertical axis titles) and text inside the axes are kept
        p = pm.get(e)
        if p is not None: p.remove(e); n += 1
    return n


_HEADLET = re.compile(r"^([A-Za-z]\d?)[\s ]{1,4}(?=\S)")


def drop_label_title(g, label):
    """Panels with the axes turned off, such as table figures: the 'above the axes' test does not work, so text starting with the original panel letter is deleted"""
    if not label or not _LAB.match(label.strip()): return 0
    lab = label.strip()
    pm = {c: p for p in g.iter() for c in p}
    n = 0
    for e in list(g.iter("{%s}text" % SVG)):
        s = "".join(e.itertext()).strip()
        m = _HEADLET.match(s)
        if m and m.group(1) == lab and pm.get(e) is not None:
            pm[e].remove(e); n += 1
    return n


def clean_figure_cell(root, fig, box):
    """Cell that uses a whole figure: (1) deletes the top large title (outside the axes, above the topmost axes),
       (2) deletes the inner panel title if there is only one, and if there are two or more removes only the leading letter (A, B) and keeps the description.
       (3) shrinks the top of the cell by the deleted amount.  Returns: the new box"""
    axes = []
    for g in fig.iter("{%s}g" % SVG):
        if (g.get("id") or "").startswith("axes_"):
            r = _axes_rect(g)
            if r: axes.append((g, r))
    if not axes: return box
    inside = set()
    for g, _ in axes:
        for e in g.iter("{%s}text" % SVG): inside.add(id(e))
    pm = {c: p for p in fig.iter() for c in p}
    top = min(r[1] for _, r in axes)
    for e in list(fig.iter("{%s}text" % SVG)):          # (1) large title
        if id(e) in inside: continue
        y = _text_y(e)
        if y is not None and y < top - 1 and pm.get(e) is not None:
            pm[e].remove(e)
    titled = []
    for g, r in axes:                                     # (2) panel titles
        ts = [e for e in g.iter("{%s}text" % SVG) if _is_title(e, r)]
        if ts: titled.append((g, r, ts))
    if len(titled) == 1:
        g, r, ts = titled[0]
        pg = {c: p for p in g.iter() for c in p}
        for e in ts:
            if pg.get(e) is not None: pg[e].remove(e)
    else:
        for g, r, ts in titled:
            for e in ts:
                for node in [e] + list(e.iter()):
                    if node.text and node.text.strip():
                        node.text = _HEADLET.sub("", node.text.lstrip(), count=1)
                        break
    ib = _ink_box(fig, (box[0], box[1], box[0] + box[2], box[1] + box[3]))   # (3) shrink the top
    if ib and ib[1] - PAD > box[1]:
        t = ib[1] - PAD
        return (box[0], t, box[2], box[3] - (t - box[1]))
    return box


def panels(path):
    """[(axes id, Element, region(x,y,w,h), axes rectangle, text list)]"""
    root = parse(path)
    _, _, W, H = viewbox(root)
    SUP = _outside_text_rows(root)
    A = []
    for g in root.iter("{%s}g" % SVG):
        i = g.get("id") or ""
        if i.startswith("axes_"):
            r = _axes_rect(g)
            if r and (r[2] - r[0]) * (r[3] - r[1]) < 0.02 * W * H:
                r2 = _content_bbox(g)          # table figure with the axes off, located again from its contents
                if r2: r = r2
            if r: A.append([i, g, r])
    out, INKC = [], {}
    for i, (nm, g, (x0, y0, x1, y1)) in enumerate(A):
        L, R, T, B = 0.0, W, 0.0, H
        for j, (_, _, (a0, b0, a1, b1)) in enumerate(A):
            if i == j: continue
            if min(y1, b1) - max(y0, b0) > 0.3 * min(y1 - y0, b1 - b0):
                if a1 <= x0: L = max(L, (a1 + x0) / 2.0)
                elif a0 >= x1: R = min(R, (x1 + a0) / 2.0)
            if min(x1, a1) - max(x0, a0) > 0.3 * min(x1 - x0, a1 - a0):
                if b1 <= y0: T = max(T, (b1 + y0) / 2.0)
                elif b0 >= y1: B = min(B, (y1 + b0) / 2.0)
        for sy in SUP:                 # the suptitle is removed from the regions (it came in pieces into each cell)
            if sy < y0 and sy + 4 > T: T = max(T, sy + 4)
        Tmid = T
        if DROP_TITLE and REGION == "mid":
            # Remove the title line of the original panel (= text above the axes) from the region.
            #   That way the SVG and the PNG both become the same 'panel without title',
            #   and there is no duplicate of the cell letter the composite adds. The title text is in the cell layout table and the legend.
            ty = [v for v in (_text_y(e) for e in g.iter("{%s}text" % SVG)) if v is not None]
            ty = [v for v in ty if v < y0 - 1]
            if ty: T = max(T, max(ty) + 1.5)
        if INK_EXPAND:
            # If its own ink (legend or labels outside the axes) goes past the midline, widen to that ink.
            #   But stop before the ink of the neighbouring panel, since the PNG is cut from the published PNG and crossing would bring a piece of the neighbour along.
            own = _ink_x(g, skip_above=(y0 - 1) if DROP_TITLE else None)
            if own:
                capL, capR = 0.0, W
                for j, (_, gj, (a0, b0, a1, b1)) in enumerate(A):
                    if i == j: continue
                    if min(B, b1) - max(T, b0) <= 0: continue          # no vertical overlap, so it does not matter
                    ink = INKC.setdefault(j, _ink_x(gj))
                    if not ink: continue
                    if a0 >= x1: capR = min(capR, ink[0] - 1.5)
                    elif a1 <= x0: capL = max(capL, ink[1] + 1.5)
                if own[1] + 1.5 > R: R = max(R, min(own[1] + 1.5, capR))
                if own[0] - 1.5 < L: L = min(L, max(own[0] - 1.5, capL))
        if REGION in ("ink", "mid+ink"):
            # Own ink measured without the title text. Deletion happens when the cell is cut out (the composite deletes it).
            ib = _ink_box(g, (x0, y0, x1, y1), skip_title=DROP_TITLE)
            if ib and REGION == "ink":
                L, T = max(0.0, ib[0] - PAD), max(0.0, ib[1] - PAD)
                R, B = min(W, ib[2] + PAD), min(H, ib[3] + PAD)
            elif ib:
                # Start from the midline region and widen only on the side where its own ink spills over (up to EXT pt on one side)
                L = max(min(L, ib[0] - PAD), L - EXT)
                R = min(max(R, ib[2] + PAD), R + EXT)
                B = min(max(B, ib[3] + PAD), B + EXT)
                if DROP_TITLE:      # the top is set to 'the top of the ink without the title' so vertical axis titles are not cut
                    T = min(max(ib[1] - PAD, Tmid - EXT), y0)
                L, T, R, B = max(0.0, L), max(0.0, T), min(W, R), min(H, B)
        out.append((nm, g, (L, T, R - L, B - T), (x0, y0, x1, y1), texts(g)))
    return root, out


def find_panel(path, title_head, panel_label=None):
    """Selects the axes by panel title (its leading part). If there is only one axes, that one is used."""
    root, ps = panels(path)
    if len(ps) == 1:                      # with only one panel, the whole figure is that panel
        return root, ps[0]
    if panel_label and ("전체" in panel_label or "단일" in panel_label):
        # the whole figure as one cell (the panel is not split)
        W_, H_ = viewbox(root)[2], viewbox(root)[3]
        fig = None
        for g in root.iter("{%s}g" % SVG):
            if (g.get("id") or "").startswith("figure_"): fig = g; break
        if fig is None: fig = root
        tx = []
        for e in fig.iter("{%s}text" % SVG):
            s = "".join(e.itertext()).strip()
            if s: tx.append(s)
        return root, ("figure", fig, (0.0, 0.0, W_, H_), (0.0, 0.0, W_, H_), tx)
    base = re.split(r"[(\u2605]", title_head)[0]
    norm = lambda s: re.sub(r"\s+", "", s)
    # (1) by the original panel letter first, the Korean and English versions are the same (only the titles are translated)
    _lm = re.match(r"^\s*([A-Za-z]\d?)", panel_label or "")
    if _lm:
        lab = _lm.group(1)
        hit = [q for q in ps
               if any(re.match(r"^%s\s" % re.escape(lab), s.strip()) for s in q[4])]   # letter followed by whitespace, to avoid category labels such as "A>G"
        if len(hit) == 1:
            return root, hit[0]
        # (2) panels with only a bare letter (P08 figures, no title, only 'a' via ax.text)
        # looked at only when the method above found nothing, so it cannot change the judgement for existing figures
        if not hit and not any(re.match(r"^[A-Za-z]\d?\s", s.strip()) for q in ps for s in q[4]):
            hit = [q for q in ps if any(s.strip() == lab for s in q[4])]
            if len(hit) == 1:
                return root, hit[0]
    for n in (16, 12, 9, 6):
        key = norm(base)[:n]
        if len(key) < 4: break
        hit = [q for q in ps if any(key in norm(s) for s in q[4])]
        if len(hit) == 1:
            return root, hit[0]
        if len(hit) > 1 and panel_label:
            lab = norm(panel_label)
            hit2 = [q for q in hit if any(norm(s).startswith(lab) for s in q[4])]
            if len(hit2) == 1:
                return root, hit2[0]
    raise LookupError("%s : '%s' 로 axes 를 못 가렸다 (후보 제목 %s)"
                      % (os.path.basename(path), title_head[:30],
                         [next((s for s in q[4] if re.match(r"^[A-Za-z]\d?\s", s)), (q[4][0] if q[4] else "?"))[:30]
                          for q in ps]))


def _ref_ids(el):
    """Collects the referenced ids used inside el (url(#x), href="#x")."""
    out = set()
    for e in el.iter():
        for k, v in e.attrib.items():
            if not v: continue
            if "url(#" in v:
                out |= set(re.findall(r"url\(#([^)]+)\)", v))
            if k.endswith("href") and v.startswith("#"):
                out.add(v[1:])
    return out


def needed_defs(root, g):
    """Fetches from the figure root the defs elements that g references (transitively).
       matplotlib puts clipPath and marker paths in the root <defs>, and without fetching them the references break."""
    byid = {}
    for e in root.iter():
        i = e.get("id")
        if i and i not in byid: byid[i] = e
    inside = {e.get("id") for e in g.iter() if e.get("id")}   # what is already inside g is not fetched
    # Deferred fix (2026-10-19 decision 'fit it in the next time you re-render for another reason'): need.pop() below takes from a set, so
    #   the clipPath order changes from run to run depending on PYTHONHASHSEED -> the composite md5 varies even with the same content (main_en 3da9c603c4 vs 47992e07ea).
    #   Fix = sort the order of taking (for example use a list with need = sorted(...) or sort order by id). Until then, compare versions by png.
    #   Until then, do not use the svg version hash either (, the shared svg version-hash tool): if ids are renumbered in order of appearance and the definition order changes, it breaks.
    #   The number of ids does not catch it either (swapping two definitions still gives 606,  test). After the fix, record the version hash and the png md5 together.
    need, got, order = _ref_ids(g) - inside, set(inside), []
    while need:
        i = need.pop()
        if i in got: continue
        got.add(i)
        e = byid.get(i)
        if e is None: continue
        order.append(copy.deepcopy(e))          # the original tree is not touched
        need |= (_ref_ids(e) - got)
    # fix (user 'redraw it now', the deferred fix above): export sorted by id instead of the order taken (set order, depends on PYTHONHASHSEED) -> deterministic output
    return sorted(order, key=lambda e: e.get("id") or "")


def prefix_ids(el, pre):
    for e in el.iter():
        if e.get("id"): e.set("id", "%s_%s" % (pre, e.get("id")))
        for k in ("{%s}href" % XLINK, "href"):
            v = e.get(k)
            if v and v.startswith("#"): e.set(k, "#%s_%s" % (pre, v[1:]))
        for k, v in list(e.attrib.items()):
            if isinstance(v, str) and "url(#" in v:
                e.set(k, re.sub(r"url\(#([^)]+)\)", lambda m: "url(#%s_%s)" % (pre, m.group(1)), v))
    return el


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    P = (__import__("os").environ["RDNA_DATA_ROOT"] + r"/for_paper/00_paper_figures_distal_layout/figures")
    for f in sys.argv[1:]:
        root, ps = panels(os.path.join(P, f))
        print("■ %s · viewBox %s · axes %d" % (f, viewbox(root), len(ps)))
        for nm, g, box, ax, tx in ps:
            head = next((s for s in tx if re.match(r"^[A-Za-z]\d?\s", s)), (tx[0] if tx else "(없음)"))
            print("   %-8s 구역(%.0f,%.0f,%.0fx%.0f) 글자 %d · %s" % (nm, box[0], box[1], box[2], box[3], len(tx), head[:52]))
