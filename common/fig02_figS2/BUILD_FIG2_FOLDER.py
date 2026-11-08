# -*- coding: utf-8 -*-
r"""Writes the finished combined Fig 2 and Supplementary into one folder (user request 2026-10-10).
   - Publication = SVG (width 414.4 pt, text stays as text and is editable)
   - The preview PNG is for checking only (the published 300 dpi PNG cut and pasted into the same place)
   - Ledger: file md5, number of cells, height, and for each cell which panel of which figure it is"""
import io, os, sys, csv, shutil, hashlib, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.stdout.reconfigure(encoding="utf-8")
import FIG2_COMPOSE as C
import FIG2_PREVIEW as PV
HERE = os.path.dirname(os.path.abspath(__file__))
P = (__import__("os").environ["RDNA_DATA_ROOT"] + r"/for_paper/00_paper_figures_distal_layout")
TAB = os.path.join(P, "tables", "FIG2_layout_2026-10-09.tsv")
OUT = os.path.join(P, "merged_Fig2_2026-10-10")
# before overwriting, keep the current combined figure in the old-version archive folder (skipped if the content is the same)
from ARCHIVE_OUT import archive_before
archive_before(OUT)
PRE = os.path.join(OUT, "preview_for_checking")
SRC = os.path.join(HERE, "OUT_FIG2")
os.makedirs(PRE, exist_ok=True)
md5 = lambda q: hashlib.md5(open(q, "rb").read()).hexdigest()
JOBS = [("Fig2", "fig2_main", "Fig 2 (본문)"),
        ("SFig2", "fig2_supp", "Supplementary Fig 2 (한 장)")]

print("■ ⑴ 합본 SVG 를 다시 만들어 바로 폴더에 (정본 표·정본 그림에서 읽는다)")
LEDGER = []
for which, stem, title in JOBS:
    for lang, d in (("ko", "그림"), ("en", "그림_영어판")):
        o = os.path.join(OUT, "%s_%s.svg" % (stem, lang))
        built, miss, w, h = C.compose(os.path.join(P, d), TAB, o, which=which, title=None)
        bad = sum(1 for i in range(len(built)) for j in range(i + 1, len(built))
                  if built[i][1] < built[j][1] + built[j][3] and built[j][1] < built[i][1] + built[i][3]
                  and built[i][2] < built[j][2] + built[j][4] and built[j][2] < built[i][2] + built[i][4])
        m = md5(o); time.sleep(0.03)
        assert m == md5(o)
        assert not miss and bad == 0, (miss, bad)
        LEDGER.append((os.path.basename(o), which, lang, str(len(built)), "%.1f" % w, "%.1f" % h,
                       str(os.path.getsize(o)), m[:12]))
        print("   %-20s 칸 %-3d 폭 %.1f 높이 %7.1f pt · 겹침 %d · 빠진 칸 %s · md5 %s"
              % (os.path.basename(o), len(built), w, h, bad, miss or "없음", m[:10]))

print("\n■ ⑵ 미리보기 PNG (검증용 · 한글·영어)")
for which, stem, title in JOBS:
    for lang, d in (("ko", "그림"), ("en", "그림_영어판")):
        o = os.path.join(PRE, "%s_%s_preview.png" % (stem, lang))
        p, size, miss, h = PV.preview(os.path.join(P, d), TAB, o, which=which)
        print("   %-28s %dx%d px · 빠진 칸 %s" % (os.path.basename(o), size[0], size[1], miss or "없음"))

print("\n■ ⑶ 칸 구성표 (배치표에서 뽑은 사본)")
rows = []
for l in io.open(TAB, encoding="utf8"):
    if l.startswith("#") or "\t" not in l: continue
    f = l.rstrip("\n").split("\t")
    if f[0] == "그림": hdr = f; continue
    rows.append(f)
g = io.open(os.path.join(OUT, "00_panel_layout.tsv"), "w", encoding="utf8", newline="\n")
g.write("# 합본 Fig 2 · Supplementary 칸 구성 — 정본 배치표 FIG2_layout_2026-10-09.tsv (md5 %s) 에서 뽑았다\n" % md5(TAB)[:12])
g.write("\t".join(hdr[:5]) + "\n")
for r in rows: g.write("\t".join(r[:5]) + "\n")
g.close()
cnt = {}
for r in rows: cnt[r[0]] = cnt.get(r[0], 0) + 1
print("   00_panel_layout.tsv · %s" % cnt)

print("\n■ ⑷ 대장")
g = io.open(os.path.join(OUT, "00_merged_manifest.tsv"), "w", encoding="utf8", newline="\n")
g.write("# 합본 Fig 2 대장 ·  · 2026-10-10 · 폭 414.4 pt 고정 · 높이 자유 · 겹침 0 · 빠진 칸 0\n")
g.write("파일\t그림\t판\t칸수\t폭pt\t높이pt\t바이트\tmd5_12\n")
for r in LEDGER: g.write("\t".join(r) + "\n")
g.close()
print("   00_merged_manifest.tsv · %d 줄" % len(LEDGER))

RM = """# 합본 Fig 2 · Supplementary Fig 2 (완성본) — 2026-10-10

## 무엇이 들어 있나
| 파일 | 내용 | 칸 |
|---|---|---|
| `fig2_main_ko.svg` · `fig2_main_en.svg` | **본문 Fig 2** | 31 |
| `fig2_suppA_ko.svg` · `_en.svg` | Supplementary Fig 2A | 10 |
| `fig2_suppB_ko.svg` · `_en.svg` | Supplementary Fig 2B | 14 |
| `fig2_suppC_ko.svg` · `_en.svg` | Supplementary Fig 2C — array 내 동질화 (포크 ) | 10 |

- **발행물은 SVG 입니다.** 폭 **414.4 pt** 고정 · 높이 자유 · `svg.fonttype=none` 이라 글자가 글자로 남아 있어 편집·번역이 됩니다.
- `preview_for_checking/*.png` 는 **눈으로 확인하는 용도**입니다 — 발행 300 dpi PNG 를 합본과 똑같은 자리·똑같은 크롭으로 오려 붙이고 칸 글자를 새로 찍은 것이라, 칸 글자의 모양은 SVG 와 다릅니다.
- `00_panel_layout.tsv` : 칸마다 어느 그림의 어느 패널인지 (정본 배치표에서 뽑음)
- `00_merged_manifest.tsv` : 파일별 칸 수 · 높이 · md5

## 어떻게 만들었나
- 발행 SVG 에서 **패널 구역만 떼어** nested `<svg viewBox>` 로 넣었습니다. **다시 그리지 않습니다** — 그림 하나를 고치면 합본을 다시 조립하기만 하면 됩니다.
- assembly code : `FIG2_COMPOSE.py` · panel extraction : `SVGPANEL.py` · preview : `FIG2_PREVIEW.py`
- 칸 차례와 줄 묶음은 **정본 배치표**(`tables/FIG2_layout_2026-10-09.tsv`)가 정합니다.

## 알아 둘 것
- 칸 글자(a, b, c1 …)는 **이 배치용으로 새로 매긴 것**이고 원본 그림의 패널 글자(A, B …)와 다릅니다. 둘 다 보입니다.
- 칸 경계는 이웃 패널의 **가운데**에서 자릅니다. 그래서 옆 패널의 축 이름 끝이 조금 들어오는 칸이 있습니다(값은 가리지 않습니다).
- 180 mm 로 줄이는 것은 **사용자가 직접** 합니다.
- 본문 그림 번호는 붙이지 않았습니다.
"""
io.open(os.path.join(OUT, "00_README.md"), "w", encoding="utf8", newline="\n").write(RM)
print("\n■ 폴더 : %s" % OUT)
for f in sorted(os.listdir(OUT)):
    p = os.path.join(OUT, f)
    print("   %-26s %s" % (f, ("폴더 %d 개" % len(os.listdir(p))) if os.path.isdir(p)
                           else "%9d B · md5 %s" % (os.path.getsize(p), md5(p)[:10])))
