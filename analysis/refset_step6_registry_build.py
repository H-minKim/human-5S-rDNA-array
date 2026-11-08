# -*- coding: utf-8 -*-
# → V11F4_DATASET_hap.tsv (286df1d6b860 · 616 행) · V11F4_DATASET_sample.tsv (48757fd4287a · 313 행)
#   단계 : ⑥ 자료 대장

# -*- coding: utf-8 -*-
# v11 자료 대장 — 하플로타입 한 줄 · 개체 한 줄 두 벌을 만든다
#   모든 칸은 실측이다. 모르는 것은 "?" 로 두고 추정하지 않는다.
#   입력 (전부 읽기만)
#     정본   08_versions\v11\02_canonical\hap_seq_node_v8.tsv              (출처·파일이름)
#     지도   08_versions\v11\03_maps\g8canonv10\CONSENSUS_SUM_G8_full.tsv (유닛수·bp·revcomp·n_break)
#     flank  08_versions\v11\12_flank\tests\FLANK_FINAL2.tsv            (양쪽 flank 실측)
#     read   08_versions\v11\08_methylation\raw\METH_READUNIT_*_{ao,jx}.tsv.gz (placed reads, methylation, data source)
#     QC     v11_readvar\02_tables\SLOT_QC.tsv                (gate·depth·eligible)
#     HiFi   03_result_tables\2026-09-07_map_test\JXHIFI_ALL.tsv.gz        (HiFi 합의 검정 자리)
#     HPRC   v11_readvar\09_inventory\hprc_r2_phasing_meta.csv (공개 조립본 release2)
import io
import os
import csv
import gzip
import glob
import collections

V = os.environ["RDNA_V11"]
F1 = os.environ["RDNA_V11_FIX1"]   # 정본 2 슬롯 수리 + 그 지도
F4 = os.environ["RDNA_V11_FIX4"]   # 메틸 raw 최신 (수정2 반영 + HG008-T CpG 전용)
RV = os.environ["RDNA_READVAR"]
OUT = os.path.dirname(os.path.abspath(__file__))
NL = "\n"
TAB = "\t"

CANCER = {"COLO829-BL", "COLO829-T", "COLO829-Tr", "COLO829-BLr",
          "HG008-T", "HG008-N-D", "HG008-N-P"}
GIAB = {"hg001", "hg002", "hg003", "hg004"}
SRCKO = {"arr": "HPRC public assembly", "hiarr": "our hifiasm assembly", "t1": "single spanning read",
         "c2seq": "placed read consensus", "pat": "T2T Q100"}

# ── 정본 출처
src = {}
cfile = {}
for x in io.open(os.path.join(V, r"02_canonical\hap_seq_node_v8.tsv"), encoding="utf8"):
    r = x.rstrip(NL).split(TAB)
    if len(r) >= 3 and r[0] != "sample":
        b = r[2].split("__")[-1]
        for t in ("hiarr_", "arr_", "t1_", "c2seq_", "pat"):
            if b.startswith(t):
                src[(r[0], r[1])] = t.rstrip("_")
                cfile[(r[0], r[1])] = r[2]
                break

# ── 지도 요약
SUM = {}
f = io.open(os.path.join(F1, r"03_maps\g8canonv11f1\CONSENSUS_SUM_G8_full.tsv"), encoding="utf8")
h = f.readline().rstrip(NL).split(TAB)
ix = {c: i for i, c in enumerate(h)}
for L in f:
    r = L.rstrip(NL).split(TAB)
    SUM[(r[0], r[1])] = (int(r[ix["N_units"]]), int(r[ix["bp"]]), r[ix["revcomp"]], int(r[ix["n_break"]]))

# ── flank
FL = {}
f = io.open(os.path.join(V, r"12_flank\tests\FLANK_FINAL2.tsv"), encoding="utf8")
h = f.readline().rstrip(NL).split(TAB)
jx = {c: i for i, c in enumerate(h)}
for L in f:
    r = L.rstrip(NL).split(TAB)
    FL[(r[0], r[1])] = r

# ── 배치 read / 메틸 / 자료원
reads = collections.defaultdict(lambda: collections.defaultdict(set))
runits = collections.defaultdict(lambda: collections.Counter())
dsrc = collections.defaultdict(collections.Counter)
for p in sorted(glob.glob(os.path.join(V, r"08_methylation\raw\METH_READUNIT_*.tsv.gz"))):
    tag = "JX" if "_jx" in os.path.basename(p) else "AO"
    g = gzip.open(p, "rt", encoding="utf8")
    hh = g.readline().rstrip(NL).split(TAB)
    kx = {c: i for i, c in enumerate(hh)}
    for L in g:
        r = L.rstrip(NL).split(TAB)
        k = (r[0], r[kx["slot"]])
        reads[k][tag].add(r[kx["read_id"]])
        runits[k][tag] += 1
        dsrc[k][r[kx["data_source"]]] += 1
    print("  읽음 %s" % os.path.basename(p), flush=True)

# ── SLOT_QC
QC = {}
f = io.open(os.path.join(RV, r"02_tables\SLOT_QC.tsv"), encoding="utf8")
h = f.readline().rstrip(NL).split(TAB)
qx = {c: i for i, c in enumerate(h)}
for L in f:
    r = L.rstrip(NL).split(TAB)
    QC[(r[0], r[1])] = r

# ── HiFi 검정
hifi = collections.Counter()
hifi5 = collections.Counter()
for L in gzip.open(os.path.join(os.environ["RDNA_RESULTS"], "2026-09-07_map_test",
                      "JXHIFI_ALL.tsv.gz"), "rt", encoding="utf8"):
    r = L.rstrip(NL).split(TAB)
    if len(r) < 8 or r[7] == "NA":
        continue
    hifi[(r[0], r[1])] += 1
    if int(r[4]) >= 5:
        hifi5[(r[0], r[1])] += 1

# ── HPRC release2 공개 조립본 (개체 단위로만 쓴다 — hap 번호는 우리 slot 과 대응하지 않는다)
HP = collections.defaultdict(list)
with io.open(os.path.join(RV, r"09_inventory\hprc_r2_phasing_meta.csv"), encoding="utf8") as fh:
    for row in csv.DictReader(fh):
        HP[row["sample_id"]].append(row)


def cohort(s):
    if s in CANCER:
        return "cancer"
    if s in GIAB:
        return "GIAB"
    if s.startswith("K562") or s.startswith("GM"):
        return "세포주"
    return "HPRC"


COLS = ["sample", "slot", "cohort", "canon_source", "canon_source_label", "canon_file",
        "N_units", "array_bp", "bp_per_unit", "revcomp", "n_break",
        "flank5_src", "flank5_len", "flank5_id", "flank3_src", "flank3_len", "flank3_id", "mid_exact",
        "ont_read_AO", "ont_read_JX", "readunit_AO", "readunit_JX", "coverage_x",
        "ont_data_source", "meth_rows",
        "qc_nread", "qc_gate", "qc_cov_frac", "qc_depth", "qc_nbreak", "qc_eligible",
        "hifi_pos", "hifi_pos_d5",
        "hprc_r2_sample", "hprc_r2_nhap", "hprc_r2_phasing", "hprc_r2_method", "hprc_r2_date"]

rows = []
for k in sorted(src):
    s, sl = k
    n, bp, rc, nb = SUM.get(k, (0, 0, "?", 0))
    fl = FL.get(k)
    a = len(reads[k]["AO"])
    j = len(reads[k]["JX"])
    ua = runits[k]["AO"]
    uj = runits[k]["JX"]
    q = QC.get(k)
    hp = HP.get(s, [])
    rows.append([
        s, sl, cohort(s), src[k], SRCKO.get(src[k], src[k]), cfile.get(k, "?"),
        n, bp, ("%.0f" % (bp / float(n))) if n else "?", rc, nb,
        fl[jx["src5"]] if fl else "?", fl[jx["len5"]] if fl else "?", fl[jx["id5"]] if fl else "?",
        fl[jx["src3"]] if fl else "?", fl[jx["len3"]] if fl else "?", fl[jx["id3"]] if fl else "?",
        fl[jx["mid_exact"]] if fl else "?",
        a, j, ua, uj, ("%.2f" % ((ua + uj) / float(n))) if n else "?",
        (dsrc[k].most_common(1)[0][0] if dsrc[k] else "none"), ua + uj,
        q[qx["nread"]] if q else "no_row", q[qx["gate"]] if q else "-",
        q[qx["cov_frac"]] if q else "-", q[qx["depth"]] if q else "-",
        q[qx["n_break"]] if q else "-", q[qx["eligible"]] if q else "-",
        hifi.get(k, 0), hifi5.get(k, 0),
        "Y" if hp else "N", len(hp),
        "/".join(sorted({x["phasing"] for x in hp})) if hp else "-",
        "/".join(sorted({x["assembly_method"] + " " + x["assembly_method_version"] for x in hp})) if hp else "-",
        "/".join(sorted({x["assembly_date"] for x in hp})) if hp else "-",
    ])

p1 = os.path.join(OUT, "V11F4_DATASET_hap.tsv")
with io.open(p1, "w", encoding="utf8", newline=NL) as g:
    g.write(TAB.join(COLS) + NL)
    for r in rows:
        g.write(TAB.join(str(x) for x in r) + NL)
print("하플로타입 표 %d 행 → %s" % (len(rows), os.path.basename(p1)))

# ── 개체 단위
bysmp = collections.defaultdict(list)
for r in rows:
    bysmp[r[0]].append(r)
C2 = ["sample", "cohort", "n_hap_canon", "slots", "canon_sources", "N_units_total", "array_bp_total",
      "ont_read_total", "coverage_min", "coverage_max", "meth_any",
      "qc_slots", "qc_eligible_slots", "hifi_pos_total",
      "hprc_r2_sample", "hprc_r2_nhap", "hprc_r2_phasing", "hprc_r2_method", "hprc_r2_date"]
srows = []
for s in sorted(bysmp):
    R = bysmp[s]
    cov = [float(x[COLS.index("coverage_x")]) for x in R if x[COLS.index("coverage_x")] != "?"]
    srows.append([
        s, R[0][2], len(R), "+".join(x[1] for x in R),
        "+".join(sorted({x[3] for x in R})),
        sum(int(x[6]) for x in R), sum(int(x[7]) for x in R),
        sum(int(x[18]) + int(x[19]) for x in R),
        ("%.2f" % min(cov)) if cov else "?", ("%.2f" % max(cov)) if cov else "?",
        "Y" if any(int(x[COLS.index("meth_rows")]) > 0 for x in R) else "N",
        sum(1 for x in R if x[COLS.index("qc_nread")] != "no_row"),
        sum(1 for x in R if x[COLS.index("qc_eligible")] == "1"),
        sum(int(x[COLS.index("hifi_pos")]) for x in R),
        R[0][COLS.index("hprc_r2_sample")], R[0][COLS.index("hprc_r2_nhap")],
        R[0][COLS.index("hprc_r2_phasing")], R[0][COLS.index("hprc_r2_method")],
        R[0][COLS.index("hprc_r2_date")],
    ])
p2 = os.path.join(OUT, "V11F4_DATASET_sample.tsv")
with io.open(p2, "w", encoding="utf8", newline=NL) as g:
    g.write(TAB.join(C2) + NL)
    for r in srows:
        g.write(TAB.join(str(x) for x in r) + NL)
print("개체 표 %d 행 → %s" % (len(srows), os.path.basename(p2)))
