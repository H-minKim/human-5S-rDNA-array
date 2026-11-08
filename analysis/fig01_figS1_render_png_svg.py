# -*- coding: utf-8 -*-
# 발행용 png·svg 를 pdf 에서 낸다 — 2026-10-21 에 새로 적은 글이다.
#
#   왜 새로 적나 : 제출본의 `Fig1_pipeline_EN.png/.svg` 와 합친 `SFig1_*.png/.svg` 는
#     그때그때 손으로 그렸고 어느 스크립트도 만들지 않았다. 그 구멍을 메운다.
#     (`SFig1_p1_*`·`SFig1_p2_*` 의 png·svg 는 그리는 스크립트가 matplotlib 으로 직접 낸다 —
#      여기서 다시 내지 않는다.)
#
#   어떻게 맞췄나 — 짐작이 아니라 발행본과 맞대어 정했다
#     · svg : `text_as_path=False` ⇒ 발행본과 바이트까지 같다 (셋 다 확인)
#     · png : 400 dpi ⇒ 발행본과 화소가 같다.
#       다만 바이트는 다르다 — 발행본은 다른 판의 pymupdf 가 썼다
#       (발행본 png 에 `Software` 덩어리가 없고 `pHYs` 가 96 dpi 로 적혀 있다).
#       ⇒ png 은 화소로 맞댄다. 바이트 md5 로 재현을 판정하면 틀린다.
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
        # 덮어쓰기 전에 전에 있던 것의 화소를 적어 둔다 — 바뀜을 조용히 넘기지 않으려고.
        before = None
        if os.path.exists(pp):
            old = fitz.Pixmap(pp)
            before = ((old.width, old.height), hashlib.md5(old.samples).hexdigest()[:12])
        pm.save(pp)
        svg = pg.get_svg_image(text_as_path=False)
        with open(ps, "w", encoding="utf-8", newline="\n") as g:
            g.write(svg)
        # 쪽 크기는 닫기 전에 꺼내 둔다 — 닫은 뒤 `pg` 를 만지면 `page is None` 으로 멈춘다.
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
