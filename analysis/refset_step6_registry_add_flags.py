# -*- coding: utf-8 -*-
# → V11F4_DATASET_hap_flags.tsv (0b7c991fada5 · 616 행 · 42 열)
#   단계 : ⑥ 자료 대장

# -*- coding: utf-8 -*-
# Flags the haplotypes whose assembly and canonical sequence are not taken at face value
#   입력 : V11F4_DATASET_hap.tsv (수정1 지도 + 수정4 메틸) + JXHIFI (HiFi 재현율)
#   깃발은 "틀렸다" 가 아니라 "이 용도로는 쓰지 말라" 다. 쓰임마다 다르므로 열을 나눠 둔다.
import io
import os
import gzip
import collections

OUT = os.path.dirname(os.path.abspath(__file__))
NL = "\n"
TAB = "\t"
SRC = os.path.join(OUT, "V11F4_DATASET_hap.tsv")

# ── HiFi 합의 재현율 (다수형 기준) — 앞서 쓴 것과 같은 코드
SUB = {75: ("G", "A"), 369: ("A", "G"), 563: ("C", "G"), 568: ("G", "A"), 1096: ("T", "G"),
       1344: ("A", "C"), 1387: ("C", "T"), 1398: ("G", "C"), 1683: ("G", "T"), 1756: ("A", "C")}


def toks(s):
    d = {}
    if s in ("-", ""):
        return d
    for x in s.split(";"):
        p = x.rsplit(":", 1)
        if len(p) == 2:
            d[p[0]] = int(p[1])
    return d


def flip(S):
    o = set(S)
    for p, (cb, mb) in SUB.items():
        cur = [t for t in S if len(t) > 1 and t[:-1] == str(p) and t[-1] in "ACGT"]
        al = cur[0][-1] if cur else cb
        for t in cur:
            o.discard(t)
        if al != mb:
            o.add("%d%s" % (p, al))
    if "437D1" in S:
        o.discard("437D1")
    elif not any(t.startswith("437D") or t.startswith("437I") for t in S):
        o.add("437REF")
    return o


agg = collections.defaultdict(lambda: [0, 0])          # TP, FN
for L in gzip.open(os.path.join(os.environ["RDNA_RESULTS"], "2026-09-07_map_test",
                      "JXHIFI_ALL.tsv.gz"), "rt", encoding="utf8"):
    f = L.rstrip(NL).split(TAB)
    if len(f) < 8 or f[7] == "NA" or int(f[4]) < 5:
        continue
    J = flip(set(x for x, v in toks(f[6]).items() if v / float(int(f[4])) >= 0.5))
    M = flip(set() if f[7] == "-" else set(f[7].split(",")))
    a = agg[(f[0], f[1])]
    a[0] += len(J & M)
    a[1] += len(M - J)
REC = {k: 100.0 * a[0] / (a[0] + a[1]) for k, a in agg.items() if a[0] + a[1] >= 30}

# ── result of the flank-to-flank test with spanning reads (10 short arrays)
SPAN_OK = {("HG02792", "2"), ("HG03567", "2"), ("HG03704", "2"),
           ("NA19468", "2"), ("HG03816", "2"), ("HG01940", "1")}
SPAN_BAD = {("NA18565", "2"): "관통 read 5 개가 전부 37~38 유닛 (정본 35) — 3 유닛 차"}
SPAN_NONE = {("NA19035", "2"), ("HG04184", "2"), ("HG04153", "2")}
# ── 정본 fasta 가 지도보다 긴 3 슬롯 (원인 미규명)
LEN_ODD = {("HG01150", "2"): "2.11x", ("NA20503", "2"): "1.54x", ("NA20827", "1"): "1.36x"}

f = io.open(SRC, encoding="utf8")
h = f.readline().rstrip(NL).split(TAB)
ix = {c: i for i, c in enumerate(h)}
rows = [L.rstrip(NL).split(TAB) for L in f if L.strip()]

NEW = ["hifi_recall", "flag_n", "flags", "flag_note"]
FLAGDESC = {
    "T1_single_spanning_read": "the reference is a single ONT read itself, do not use for unit-level variant analysis (confirmed not replaceable)",
    "C2_read_consensus": "slot whose reference was replaced by a placed read consensus, not an assembly",
    "READ0_no_placement": "0 placed reads across both AO and JX, so there is no independent support",
    "COV_below_2x": "placed read-units divided by reference units < 2",
    "QC_no_row": "no row in SLOT_QC (0 placed reads), so there is no quality call at all",
    "QC_not_eligible": "SLOT_QC eligible=0, excluded from the zone and methylation analyses",
    "FLANK_both_missing": "neither flank recovered from sequence (mid_exact != Y), so the array boundary is unconfirmed",
    "FLANK_one_missing": "one flank not recovered (the middle span is confirmed)",
    "HIFI_untested": "no position with depth >= 5 to compare against the HiFi consensus",
    "HIFI_low_recall": "HiFi consensus recall < 95%",
    "SHORT_under_30_units": "reference unit count < 30, a length the eligible rule also excludes",
    "SPAN_mismatch_with_reference": "the spanning read has a different unit count from the reference, a candidate for reference repair",
    "LEN_longer_than_map": "the reference fasta is longer than the span the map covers, cause not established",
    "HAP_only_one": "only one reference haplotype for that individual",
}
nhap = collections.Counter(r[0] for r in rows)
out = []
for r in rows:
    k = (r[0], r[1])
    fl = []
    cs = r[ix["canon_source"]]
    cov = r[ix["coverage_x"]]
    covf = float(cov) if cov not in ("?", "") else -1.0
    if cs == "t1":
        fl.append("T1_single_spanning_read")
    if cs == "c2seq":
        fl.append("C2_read_consensus")
    if covf == 0:
        fl.append("READ0_no_placement")
    elif 0 <= covf < 2:
        fl.append("COV_below_2x")
    if r[ix["qc_nread"]] == "no_row":
        fl.append("QC_no_row")
    elif r[ix["qc_eligible"]] == "0":
        fl.append("QC_not_eligible")
    if r[ix["mid_exact"]] != "Y":
        fl.append("FLANK_both_missing")
    elif r[ix["flank5_src"]] in ("-", "?") or r[ix["flank3_src"]] in ("-", "?"):
        fl.append("FLANK_one_missing")
    if int(r[ix["hifi_pos_d5"]]) == 0:
        fl.append("HIFI_untested")
    rc = REC.get(k)
    if rc is not None and rc < 95:
        fl.append("HIFI_low_recall")
    if int(r[ix["N_units"]]) < 30:
        fl.append("SHORT_under_30_units")
    if k in SPAN_BAD:
        fl.append("SPAN_mismatch_with_reference")
    if k in LEN_ODD:
        fl.append("LEN_longer_than_map")
    if nhap[r[0]] == 1:
        fl.append("HAP_only_one")
    out.append(r + ["%.1f" % rc if rc is not None else "-", str(len(fl)),
                    ";".join(fl) if fl else "-",
                    " / ".join(FLAGDESC[x] for x in fl) if fl else "-"])

p = os.path.join(OUT, "V11F4_DATASET_hap_flags.tsv")
with io.open(p, "w", encoding="utf8", newline=NL) as g:
    g.write(TAB.join(h + NEW) + NL)
    for r in out:
        g.write(TAB.join(str(x) for x in r) + NL)
print("→ %s  (%d 행 · 열 %d)" % (os.path.basename(p), len(out), len(h) + len(NEW)))

HP = [r for r in out if r[ix["cohort"]] == "HPRC"]
print()
print("HPRC %d 하플로 · 깃발 분포" % len(HP))
c = collections.Counter(x for r in HP for x in (r[-2].split(";") if r[-2] != "-" else []))
for k, v in c.most_common():
    print("   %-20s %3d   %s" % (k, v, FLAGDESC[k][:46]))
print()
d = collections.Counter(int(r[-3]) for r in HP)
for k in sorted(d):
    print("   깃발 %d 개 : %3d 하플로" % (k, d[k]))
print("   깃발 0 (문제 안 잡힘) : %d / %d (%.1f%%)" % (d[0], len(HP), 100.0 * d[0] / len(HP)))
