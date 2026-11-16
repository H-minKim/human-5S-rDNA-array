# -*- coding: utf-8 -*-
# Render the publication png and svg from the pdf.
#
#   The submitted Fig1_pipeline_EN.png/.svg and the merged SFig1_*.png/.svg
#     were produced by hand and by no script; this fills that gap.
#     (The png and svg for SFig1_p1_* and SFig1_p2_* are written directly by matplotlib in the drawing script
#      and are not re-rendered here.)
#
#   The settings were fixed by comparison against the published files, not by assumption:
#     · svg: with text_as_path=False the output is byte-identical to the published file, confirmed for all three.
#     · png: at 400 dpi the pixels match the published file.
#       The bytes differ, because the published file was written by a different release of pymupdf
#       (it carries no Software chunk and records pHYs as 96 dpi).
#       So png files are compared by pixels; judging reproduction by a byte md5 would be wrong.
import os, sys, hashlib
import fitz
sys.stdout.reconfigure(encoding="utf-8")

WORK_ROOT = os.environ["WORK_ROOT"]
OUT_DIR = os.environ.get("OUT_DIR",
                         os.path.join(WORK_ROOT, "04_figures", "paper_figures_final_2026-10-07",
                                      "01_selected_for_paper"))
DPI = int(os.environ.get("PNG_DPI", "400"))
NAMES = ("Fig1_pipeline_EN", "SFig1_EN", "SFig1_KO")


def main():
    print("png·svg 내기 — %d dpi · %s" % (DPI, OUT_DIR))
    n_ok = 0
    for n in NAMES:
        src = os.path.join(OUT_DIR, n + ".pdf")
        if not os.path.exists(src):
            print("   없음 : %s.pdf" % n)
            continue
        d = fitz.open(src)
        pg = d[0]
        pm = pg.get_pixmap(dpi=DPI)
        pp = os.path.join(OUT_DIR, n + ".png")
        ps = os.path.join(OUT_DIR, n + ".svg")
        # Before overwriting, the pixel dimensions of the existing file are recorded, so that a change cannot pass silently.
        before = None
        if os.path.exists(pp):
            old = fitz.Pixmap(pp)
            before = ((old.width, old.height), hashlib.md5(old.samples).hexdigest()[:12])
        pm.save(pp)
        svg = pg.get_svg_image(text_as_path=False)
        with open(ps, "w", encoding="utf-8", newline="\n") as g:
            g.write(svg)
        # The page size is read before closing; touching the page object afterwards fails with "page is None".
        _w, _h = pg.rect.width, pg.rect.height
        d.close()
        now = ((pm.width, pm.height), hashlib.md5(pm.samples).hexdigest()[:12])
        tag = "처음" if before is None else ("화소 같음" if before == now else "화소 바뀜 %s -> %s" % (before[1], now[1]))
        print("   %-18s %.1f x %.1f pt -> %dx%d px · svg %d 자   %s"
              % (n, _w, _h, pm.width, pm.height, len(svg), tag))
        n_ok += 1
    print("낸 것 %d / %d" % (n_ok, len(NAMES)))


if __name__ == "__main__":
    main()
