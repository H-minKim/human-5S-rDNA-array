# -*- coding: utf-8 -*-
# supplementary 를 한 쪽으로 합친다
#   위 : `SFig1_p1_*` (a·b·d / c·e)   아래 : `SFig1_p2_*` (f / g / h / i)
#   폭은 규격 414.4 pt 로 같다 — 세로로만 잇는다. 축소하지 않는다.
#   높이는 자유 — 겹침만 없으면 된다.
import os, sys, hashlib, fitz
sys.stdout.reconfigure(encoding="utf-8")
# -- (repo) paths come from the environment --------------------------------
#   기본값이 분석 때 쓴 경로와 같다 ⇒ 작업 폴더에서는 전과 똑같이 돈다.
#   그림 내용은 한 글자도 바뀜이 없다 — png 화소로 맞대어 확인한다.
WORK_ROOT = os.environ["WORK_ROOT"]
OUT_DIR = os.environ.get("OUT_DIR",
                         os.path.join(WORK_ROOT, "04_figures", "paper_figures_final_2026-10-07",
                                      "01_selected_for_paper"))
D = OUT_DIR
GAP = 6.0            # 두 쪽 사이 숨 (pt)

for lang in ("EN", "KO"):
    p1 = os.path.join(D, "SFig1_p1_%s.pdf" % lang)
    p2 = os.path.join(D, "SFig1_p2_%s.pdf" % lang)
    if not (os.path.exists(p1) and os.path.exists(p2)):
        print("%s : 쪽 파일이 없다" % lang)
        continue
    s1, s2 = fitz.open(p1), fitz.open(p2)
    r1, r2 = s1[0].rect, s2[0].rect
    W = max(r1.width, r2.width)
    H = r1.height + GAP + r2.height
    out = fitz.open()
    pg = out.new_page(width=W, height=H)
    pg.show_pdf_page(fitz.Rect((W - r1.width) / 2, 0,
                               (W + r1.width) / 2, r1.height), s1, 0)
    pg.show_pdf_page(fitz.Rect((W - r2.width) / 2, r1.height + GAP,
                               (W + r2.width) / 2, H), s2, 0)
    s1.close(); s2.close()
    O = os.path.join(D, "SFig1_%s.pdf" % lang)
    out.save(O + ".tmp", garbage=3, deflate=True)
    out.close()
    os.replace(O + ".tmp", O)
    d = fitz.open(O)
    ls = sorted({sp["text"].strip() for bl in d[0].get_text("dict")["blocks"]
                 for ln in bl.get("lines", []) for sp in ln.get("spans", [])
                 if len(sp["text"].strip()) == 1 and sp["text"].strip().isalpha()
                 and sp["text"].strip().islower() and sp["size"] > 10.5})
    small = sum(1 for bl in d[0].get_text("dict")["blocks"] for ln in bl.get("lines", [])
                for sp in ln.get("spans", []) if sp["text"].strip() and sp["size"] < 7.95)
    print("SFig1_%s.pdf  %.1f x %.1f pt · 쪽 %d · 패널 %s · 8pt미만 %d · md5 %s"
          % (lang, W, H, d.page_count, " ".join(ls), small,
             hashlib.md5(open(O, "rb").read()).hexdigest()[:12]))
    # 같은 줄 가로겹침
    sp = [(fitz.Rect(x["bbox"]), x["text"].strip(), x["size"])
          for bl in d[0].get_text("dict")["blocks"] for ln in bl.get("lines", [])
          for x in ln.get("spans", []) if x["text"].strip()]
    hit = 0
    for i in range(len(sp)):
        for j in range(i + 1, len(sp)):
            a, b = sp[i], sp[j]
            if abs((a[0].y0 + a[0].y1) / 2 - (b[0].y0 + b[0].y1) / 2) > 0.55 * max(a[2], b[2]):
                continue
            ov = min(a[0].x1, b[0].x1) - max(a[0].x0, b[0].x0)
            if ov > 0.5:
                print("   겹침 %.1f : %-26s <-> %s" % (ov, a[1][:26], b[1][:26]))
                hit += 1
    print("   글자 %d · 같은 줄 겹침 %d" % (len(sp), hit))
    d.close()
