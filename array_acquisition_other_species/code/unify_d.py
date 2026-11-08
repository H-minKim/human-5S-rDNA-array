# -*- coding: utf-8 -*-
# 타 생물종 5S 산출물 형식 통일 (D 디스크 재편)
#   열쇠 : array_id = <종키>__<염색체>_<시작>_<끝>  (모든 표가 이 하나로 이어진다)
#   종키 : 조립본 단위 (human_CHM13v2 · zfinch_v03_dip · zebrafish_GRCz12tu · arabidopsis_ColCC ...)
#   모든 read/유닛/array 표 앞에 species · sample 열을 붙인다 (표만 보고도 출처를 알게)
#   폴더 : 00_catalog · 01_canonical(정본) · 02_consensus_map · 03_reads · 04_bam · 05_bedmethyl · 06_aux
#   옛 폴더(zebrafish · zfinch_ont · arabidopsis_ont 등)는 새 자리로 복사만 하고 지우지 않는다
# MIT License · Copyright (c) 2026 Hungmin Kim
import io
import os
import sys
import csv
import gzip
import shutil

NL = chr(10)
TAB = chr(9)
D = sys.argv[1] if len(sys.argv) > 1 else "$RDNA_OTHERSP_DIR"
OLD = os.path.join(D, "data")
NEW = os.path.join(D, "data_v2")
SUB = ["00_catalog", "01_canonical", "02_consensus_map", "03_reads", "04_bam", "05_bedmethyl", "06_aux", "cms_tbl"]
for s in SUB:
    os.makedirs(os.path.join(NEW, s), exist_ok=True)
# 종키 : 옛 폴더 이름 -> 새 종키, 그리고 표본 이름
SPK = {"zebrafish": "zebrafish_GRCz12tu", "zfinch_ont": "zfinch_v03_dip", "arabidopsis_ont": "arabidopsis_ColCC"}
# 제브라피시 read 표는 옛 스크립트(9 열)라 쓰지 않는다 — 같은 units_matrix.py 로 다시 계산한 원본을 쓴다
ZFU2 = "$SCRATCH_DIR/relay/zfu2"
LOG = []


def aid(sp, chrom, s, e):
    return "%s__%s_%s_%s" % (sp, chrom, s, e)


def read_tsv(p):
    with io.open(p, encoding="utf8") as f:
        rd = list(csv.reader(f, delimiter=TAB))
    return rd[0], rd[1:]


def write_tsv(p, hd, rows):
    with io.open(p, "w", encoding="utf8", newline=NL) as f:
        f.write(TAB.join(hd) + NL)
        for r in rows:
            f.write(TAB.join(r) + NL)
    LOG.append("%s  %d 행" % (os.path.relpath(p, D), len(rows)))


# 1) 카탈로그 : 두 파일 합치고 array_id 붙이기
SPMAP = {"Danio_rerio": "zebrafish_GRCz12tu", "Taeniopygia_guttata": "zfinch_v03_dip",
         "Arabidopsis_thaliana": "arabidopsis_ColCC", "Zea_mays": "maize_Mo17_T2T"}
rows, hd = [], None
for fn in ("5S_array_catalog.tsv", "5S_array_catalog_primates_mouse.tsv"):
    p = os.path.join(OLD, fn)
    if not os.path.exists(p):
        continue
    h, rs = read_tsv(p)
    hd = h
    for r in rs:
        sp = r[0]
        asm = r[1]
        key = SPMAP.get(sp)
        if key is None:
            key = sp if "_" in sp else sp
        # 영장류·마우스 카탈로그는 1 열이 이미 종키
        if sp in ("chimp_mPanTro3", "bonobo_mPanPan1", "gorilla_mGorGor1", "orang_sumatran_mPonAbe1",
                  "orang_bornean_mPonPyg2", "siamang_mSymSyn1", "rhesus_T2T_MMU8v2", "mouse_C57BL6J_T2T",
                  "human_CHM13v2"):
            key = sp
        elif sp == "Arabidopsis_thaliana" and "CEN" in asm:
            key = "arabidopsis_ColCEN12"
        elif sp == "Taeniopygia_guttata" and ("GCA_048771995" in asm or "GCA_048772025" in asm):
            # mat·pat 두 조립본을 합쳐 만든 dip 참조로 구조표·지도를 냈고 염색체 이름·좌표가 같으므로 종키를 맞춘다
            key = "zfinch_v03_dip"
        elif sp == "Danio_rerio" and "GRCz11" in asm:
            key = "zebrafish_GRCz11"
        rows.append([key, aid(key, r[2], r[3], r[4])] + r)
write_tsv(os.path.join(NEW, "00_catalog", "5S_array_catalog.tsv"),
          ["species_key", "array_id"] + hd, rows)
# 2) 구조표 : array_id 붙이기 (unit_identity_pct 열은 정렬 없이 잰 값이라 뺀다)
p = os.path.join(OLD, "array_structure", "5S_array_structure_all.tsv")
if os.path.exists(p):
    h, rs = read_tsv(p)
    drop = h.index("unit_identity_pct") if "unit_identity_pct" in h else -1
    nh = [x for i, x in enumerate(h) if i != drop]
    out = []
    for r in rs:
        out.append([aid(r[0], r[1], r[2], r[3])] + [x for i, x in enumerate(r) if i != drop])
    write_tsv(os.path.join(NEW, "00_catalog", "5S_array_structure.tsv"), ["array_id"] + nh, out)
# 3) 그대로 옮기는 표
for fn in ("dataset_access_status.tsv", "5S_methylation_cross_species_summary.tsv", "rhesus_5S_array_judgement.txt"):
    p = os.path.join(OLD, fn)
    if os.path.exists(p):
        shutil.copy2(p, os.path.join(NEW, "00_catalog", fn))
        LOG.append("00_catalog/" + fn)
# 4) read / 유닛 / array 표 : species · sample 열 앞에 붙여 복사
def add_cols(src, dst, sp, smp):
    op = gzip.open if src.endswith(".gz") else io.open
    wp = gzip.open if dst.endswith(".gz") else io.open
    with op(src, "rt", encoding="utf8") as f, wp(dst, "wt", encoding="utf8", newline=NL) as g:
        h = f.readline().rstrip(NL)
        g.write("species_key" + TAB + "sample" + TAB + h + NL)
        n = 0
        for L in f:
            L = L.rstrip(NL)
            if not L:
                continue
            s = smp if smp else L.split(TAB)[0]
            g.write(sp + TAB + s + TAB + L + NL)
            n += 1
    LOG.append("%s  %d 행" % (os.path.relpath(dst, D), n))


for old_dir, sp in SPK.items():
    src = os.path.join(OLD, old_dir)
    if not os.path.isdir(src):
        continue
    rd = os.path.join(NEW, "03_reads", sp)
    os.makedirs(rd, exist_ok=True)
    for base in (src, os.path.join(src, "read_unit")):
        if not os.path.isdir(base):
            continue
        for fn in sorted(os.listdir(base)):
            p = os.path.join(base, fn)
            if not os.path.isfile(p):
                continue
            if sp == "zebrafish_GRCz12tu" and os.path.isdir(ZFU2):
                continue
            if fn.endswith(".read_unit.tsv") or fn.endswith(".read_unit.tsv.gz"):
                smp = fn.split(".read_unit")[0]
                add_cols(p, os.path.join(rd, "%s__%s.read_unit.tsv%s" % (sp, smp, ".gz" if fn.endswith(".gz") else "")), sp, smp)
            elif fn in ("5S_unit_summary.tsv", "5S_array_summary.tsv"):
                add_cols(p, os.path.join(rd, fn), sp, None)
    # BAM · bedMethyl · 부가
    for sub, tgt in (("bam", "04_bam"), ("bedmethyl", "05_bedmethyl"), ("extract", "06_aux")):
        s2 = os.path.join(src, sub)
        if os.path.isdir(s2):
            t2 = os.path.join(NEW, tgt, sp)
            os.makedirs(t2, exist_ok=True)
            for fn in sorted(os.listdir(s2)):
                shutil.copy2(os.path.join(s2, fn), os.path.join(t2, fn))
            LOG.append("%s/%s  파일 %d" % (tgt, sp, len(os.listdir(s2))))
    a2 = os.path.join(NEW, "06_aux", sp)
    os.makedirs(a2, exist_ok=True)
    for fn in sorted(os.listdir(src)):
        p = os.path.join(src, fn)
        if os.path.isfile(p) and (fn.endswith(".txt") or fn in ("poscmp.tsv", "summary.tsv", "pileup_summary.txt",
                                                               "flowcell_basecall.tsv", "array_read_selection.tsv",
                                                               "batch_basecall.tsv")):
            shutil.copy2(p, os.path.join(a2, fn))
# 4b) 제브라피시 read 표 : 새 형식 원본에서
if os.path.isdir(ZFU2):
    sp = "zebrafish_GRCz12tu"
    rd = os.path.join(NEW, "03_reads", sp)
    os.makedirs(rd, exist_ok=True)
    for fn in sorted(os.listdir(ZFU2)):
        p2 = os.path.join(ZFU2, fn)
        if not os.path.isfile(p2):
            continue
        if fn.endswith(".read_unit.tsv"):
            smp = fn.split(".read_unit")[0]
            add_cols(p2, os.path.join(rd, "%s__%s.read_unit.tsv" % (sp, smp)), sp, smp)
        elif fn in ("5S_unit_summary.tsv", "5S_array_summary.tsv"):
            add_cols(p2, os.path.join(rd, fn), sp, None)
# 5) 금화조 PacBio 집계 · BED 등 부가
for fn in ("zfinch_5S_array_meth_summary.tsv", "zfinch_5S_units_meth.tsv", "grcz11_5S_regions.bed", "grcz12_5S_arrays.bed"):
    p = os.path.join(OLD, fn)
    if os.path.exists(p):
        os.makedirs(os.path.join(NEW, "06_aux"), exist_ok=True)
        shutil.copy2(p, os.path.join(NEW, "06_aux", fn))
        LOG.append("06_aux/" + fn)
p = os.path.join(OLD, "cms_tbl")
if os.path.isdir(p):
    for fn in sorted(os.listdir(p)):
        shutil.copy2(os.path.join(p, fn), os.path.join(NEW, "cms_tbl", fn))
    LOG.append("cms_tbl  파일 %d" % len(os.listdir(p)))
print(NL.join(LOG))
