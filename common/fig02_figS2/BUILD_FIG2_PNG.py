# -*- coding: utf-8 -*-
r"""Writes the combined Fig 2 and Supplementary as 300 dpi PNG.
   - Positions and crops are computed exactly as for the SVG composite (ROWS, width and margins of FIG2_COMPOSE unchanged)
   - Cell contents are cut and pasted from the published 300 dpi PNG, and the cell letters are drawn in Malgun Gothic
   - Width 414.4 pt = 5.756 in -> 1,727 px (the dpi metadata is also set to 300)"""
import io, os, sys, shutil
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.stdout.reconfigure(encoding="utf-8")
from PIL import Image, ImageDraw, ImageFont
import SVGPANEL as S
import FIG2_COMPOSE as C
import FIG2_PREVIEW as PV

P = (__import__("os").environ["RDNA_DATA_ROOT"] + r"/for_paper/00_paper_figures_distal_layout")
TAB = os.path.join(P, "tables", "FIG2_layout_2026-10-09.tsv")
OUT = os.path.join(P, "merged_Fig2_2026-10-10")
# before overwriting, keep the current combined figure in the old-version archive folder (skipped if the content is the same)
from ARCHIVE_OUT import archive_before
archive_before(OUT)
DPI = 300
PX = DPI / 72.0
FONT = (__import__("os").environ.get("RDNA_FONT_DIR", __import__("os").path.join(__import__("os").environ.get("SystemRoot", ""), "Fonts")) + "/malgun.ttf")


CHROME = __import__("os").environ.get("CHROME", __import__("os").path.join(__import__("os").environ.get("ProgramFiles", ""), "Google", "Chrome", "Application", "chrome.exe"))
PROF = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_chrome_prof")


def build(figdir, which, out):
    """Renders the SVG composite (the .svg of the same name in the same folder) as is with Chrome, so the SVG and the PNG are the same figure."""
    import re, math, subprocess
    svg = out[:-4] + ".svg"
    head = io.open(svg, encoding="utf8").read(3000)
    wpt = float(re.search(r'width="([\d.]+)pt"', head).group(1))
    hpt = float(re.search(r'height="([\d.]+)pt"', head).group(1))
    wpx, hpx = wpt * 96 / 72, hpt * 96 / 72                  # CSS px
    tw, th = int(round(wpt * PX)), int(round(hpt * PX))      # target 300 dpi
    LOC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_png_local")
    os.makedirs(LOC, exist_ok=True)
    tmp = os.path.join(LOC, os.path.basename(out)[:-4] + "_raw.png")   # render in the working folder
    # --disable-lcd-text, without it vertical text (axis titles) gets LCD subpixel colour fringes (purple, yellow, sky blue), measured (87,0,87) and others
    r = subprocess.run([CHROME, "--headless=new", "--disable-gpu", "--disable-lcd-text", "--hide-scrollbars",
                        "--no-first-run", "--no-default-browser-check", "--user-data-dir=" + PROF,
                        "--force-device-scale-factor=%.6f" % (tw / wpx),
                        "--window-size=%d,%d" % (math.ceil(wpx), math.ceil(hpx)),
                        "--screenshot=" + tmp, "file:///" + svg.replace("\\", "/")],
                       capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=180)
    assert r.returncode == 0 and os.path.exists(tmp), "Chrome 렌더 실패 : %s" % (r.stderr or "")[-300:]
    im = Image.open(tmp).convert("RGB")
    bg = Image.new("RGB", (max(im.width, tw), max(im.height, th)), "white")
    bg.paste(im, (0, 0))
    img = bg.crop((0, 0, tw, th))
    loc = os.path.join(LOC, os.path.basename(out))
    img.save(loc, dpi=(DPI, DPI))                    # save to the working folder first
    os.remove(tmp)
    import hashlib, time, shutil
    h = hashlib.md5(open(loc, "rb").read()).hexdigest()
    for att in range(4):                              # copy finished files to F: one at a time and check md5
        try:
            shutil.copyfile(loc, out); time.sleep(0.3)
            if hashlib.md5(open(out, "rb").read()).hexdigest() == h: break
        except OSError as e:
            if getattr(e, 'errno', None) == 22:                 # absolute rule: Errno 22 = USB hub fault -> stop at once, no retry (2026-10-18 night)
                raise SystemExit('Errno 22 — 즉시 중단 (재시도 안 함) : %s' % e)
            print("      %s 복사 시도 %d : %s" % (os.path.basename(out), att + 1, e)); time.sleep(3.0)
    else:
        raise OSError("%s 를 F: 에 못 썼다 (네 번)" % os.path.basename(out))
    # request 'vector pdf for NAR submission, in the same call as the png': print-to-pdf of the same svg with the same Chrome (page size = svg size in pt, margin 0)
    #   the pdf contains its creation time, so its md5 changes every time -> judge versions by the png
    body = io.open(svg, encoding="utf8").read()
    body = re.sub(r"^\s*<\?xml[^>]*\?>", "", body)
    body = re.sub(r"<!DOCTYPE[^>]*>", "", body, count=1)
    htm = os.path.join(LOC, os.path.basename(out)[:-4] + "_pdf.html")
    io.open(htm, "w", encoding="utf8", newline="\n").write(
        "<!doctype html><html><head><meta charset='utf-8'><style>@page{size:%.4fpt %.4fpt;margin:0}"
        "html,body{margin:0;padding:0;background:#fff}svg{display:block;width:%.4fpt;height:%.4fpt}</style></head><body>%s</body></html>"
        % (wpt, hpt, wpt, hpt, body))
    lpdf = os.path.join(LOC, os.path.basename(out)[:-4] + ".pdf")
    if os.path.exists(lpdf): os.remove(lpdf)
    r = subprocess.run([CHROME, "--headless=new", "--disable-gpu", "--no-first-run", "--no-default-browser-check", "--user-data-dir=" + PROF,
                        "--no-pdf-header-footer", "--print-to-pdf-no-header", "--print-to-pdf=" + lpdf, "file:///" + htm.replace("\\", "/")],
                       capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=180)
    assert r.returncode == 0 and os.path.exists(lpdf), "Chrome pdf 실패 : %s" % (r.stderr or "")[-300:]
    os.remove(htm)
    # Chrome rounds the page size to 0.72 pt steps (measured 423.50 x 643.32 -> 424.08 x 642.96), and the figure starts at the top left of the page, so
    #   reset the page boxes (MediaBox, CropBox) to the svg size anchored at the top left (the drawing itself is not touched)
    from pypdf import PdfReader, PdfWriter
    from pypdf.generic import RectangleObject
    rd = PdfReader(lpdf); assert len(rd.pages) == 1, len(rd.pages)
    pg = rd.pages[0]; top = float(pg.mediabox.top); left = float(pg.mediabox.left)
    box = RectangleObject([left, top - hpt, left + wpt, top])
    pg.mediabox = box; pg.cropbox = box
    wr = PdfWriter(); wr.add_page(pg)
    with open(lpdf, "wb") as g: wr.write(g)
    opdf = out[:-4] + ".pdf"
    hp = hashlib.md5(open(lpdf, "rb").read()).hexdigest()
    try:
        shutil.copyfile(lpdf, opdf); time.sleep(0.3)
    except OSError as e:
        if getattr(e, 'errno', None) == 22:
            raise SystemExit('Errno 22 — 즉시 중단 (재시도 안 함) : %s' % e)
        raise
    assert hashlib.md5(open(opdf, "rb").read()).hexdigest() == hp, "pdf 복사 md5 다름"
    n = len(C.load_plan(TAB, which))
    return img.size, hpt, n


JOBS = [("Fig2", "fig2_main"), ("SFig2", "fig2_supp")]
print("■ 300 dpi PNG")
for which, stem in JOBS:
    for lang, d in (("ko", "그림"), ("en", "그림_영어판")):
        o = os.path.join(OUT, "%s_%s.png" % (stem, lang))
        size, H, n = build(os.path.join(P, d), which, o)
        print("   %-20s %4d x %5d px · 칸 %2d · %.1f x %.1f pt · %7.2f MB"
              % (os.path.basename(o), size[0], size[1], n, C.W, H, os.path.getsize(o) / 1e6))

old = os.path.join(OUT, "preview_for_checking")
if os.path.isdir(old):
    shutil.rmtree(old)
    print("\n   옛 저해상 미리보기 폴더는 지웠다 (위 300 dpi PNG 로 갈음)")
