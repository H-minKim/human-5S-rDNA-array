# -*- coding: utf-8 -*-
"""유닛 보조표 v2 : 유전자부 β · array_pos(배치기 정렬 좌표) · 고정 구역 (2026-10-12)

기존 UNIT_<label>.tsv 와 같은 unit_key("smp|read|k")를 쓴다 -> 이어 붙일 수 있다.
  유전자부 = 유닛 끝 119 bp = [TSS, TSS+118]  (TSS = 유닛 안 마지막 GTCTACGGC · 유닛 = ... + 118)
  array_pos : 배치기 BAM 의 정렬(참조 K562r_h1 등)에서 TSS 의 참조 좌표 -> MAP_G8 의 유닛 구간
              (1 = 원위 · N = 근위). 정렬이 방향·클립을 이미 갖고 있다.
usage: python cell_units2.py LABEL BAM SMP [MAPFILE MAPSAMPLE]   (여러 BAM 은 ; 로 이어 쓴다)
"""
import io
import os
import sys
import bisect
import collections
import pysam

NL = chr(10)
TB = chr(9)
ANCH = "GTCTACGGC"
MINGAP = 500
MLCUT = 128
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
LABEL = sys.argv[1]
BAMS = sys.argv[2].split(";")
SMPS = sys.argv[3].split(";")
MAPF = sys.argv[4] if len(sys.argv) > 4 else None
MAPS = sys.argv[5] if len(sys.argv) > 5 else None
_SELF = os.path.basename(os.path.abspath(__file__))


def _need(k):
    # No defaults are set: if a variable is missing the script stops with a readable message
    v = os.environ.get(k, "")
    if not v:
        raise SystemExit(
            "environment variable %s is not set; it is required by %s." % (k, _SELF)
            + chr(10)
            + "The deposited copy has no default; see README.md section 2")
    return v


OUT = (_need("RDNA_CELLS_DIR") if LABEL != "HEPG2" else _need("RDNA_HEPG2_DIR"))


def rc(s):
    return "".join(COMP.get(c, "N") for c in reversed(s))


def anchors(seq):
    hits = []
    i = seq.find(ANCH)
    while i >= 0:
        hits.append(i)
        i = seq.find(ANCH, i + 1)
    return [h for n, h in enumerate(hits) if not (n + 1 < len(hits) and hits[n + 1] - h < MINGAP)]


def parse_c(seq, mm, ml):
    out = {}
    k = 0
    for b in [b for b in mm.rstrip(";").split(";") if b]:
        head, _, body = b.partition(",")
        deltas = [int(v) for v in body.split(",") if v != ""] if body else []
        if head.startswith("C+m"):
            idx = [i for i, c in enumerate(seq) if c == "C"]
            j = 0
            for d in deltas:
                j += d
                if j >= len(idx):
                    break
                out[idx[j]] = ml[k] if k < len(ml) else 0
                k += 1
                j += 1
        else:
            k += len(deltas)
    return out


# 유닛 지도
UMAP = {}
NSLOT = {}
if MAPF:
    f = io.open(MAPF, encoding="utf8")
    ix = dict((n, i) for i, n in enumerate(f.readline().rstrip(NL).split(TB)))
    tmp = collections.defaultdict(list)
    for l in f:
        x = l.rstrip(NL).split(TB)
        if x[ix["sample"]] != MAPS or x[ix["seg_type"]] != "UNIT":
            continue
        tmp[x[ix["slot"]]].append((int(x[ix["unit_start"]]), int(x[ix["unit_end"]]),
                                   int(x[ix["array_pos"]])))
        NSLOT[x[ix["slot"]]] = int(x[ix["N"]])
    for s, v in tmp.items():
        v.sort()
        UMAP[s] = ([a for a, _, _ in v], v)
    print("지도 슬롯 :", dict((s, (len(v[1]), NSLOT[s])) for s, v in UMAP.items()))


def q2r(r, qpos):
    """저장된 query 좌표 -> 참조 좌표(0-based) · 클립/삽입이면 None."""
    q = 0
    rp = r.reference_start
    for op, ln in r.cigartuples or []:
        if op in (0, 7, 8):
            if q <= qpos < q + ln:
                return rp + (qpos - q)
            q += ln
            rp += ln
        elif op in (1, 4):
            if q <= qpos < q + ln:
                return None
            q += ln
        elif op in (2, 3):
            rp += ln
    return None


fo = io.open(os.path.join(OUT, "UNIT2_%s.tsv" % LABEL), "w", encoding="utf8", newline=NL)
fo.write(TB.join(["unit_key", "read_id", "k", "n_cpg_gene", "n_meth_gene", "beta_gene",
                  "slot", "array_pos", "N", "zone_fixed"]) + NL)
ST = collections.Counter()
SEEN = set()
for bp, smp in zip(BAMS, SMPS):
    bam = pysam.AlignmentFile(bp, "rb", check_sq=False)
    for r in bam.fetch(until_eof=True):
        if r.is_secondary or r.is_supplementary or not r.query_sequence:
            continue
        if (smp, r.query_name) in SEEN:
            continue
        SEEN.add((smp, r.query_name))
        seq = r.query_sequence
        orig = rc(seq) if r.is_reverse else seq
        try:
            mm = r.get_tag("MM")
            ml = list(r.get_tag("ML"))
        except KeyError:
            continue
        nf, nr = orig.count(ANCH), rc(orig).count(ANCH)
        if max(nf, nr) < 2:
            continue
        L = len(orig)
        cm = parse_c(orig, mm, ml)
        aflip = nr > nf
        if aflip:
            s2 = rc(orig)
            cm = dict((L - 2 - p, q) for p, q in cm.items() if L - 2 - p >= 0)
        else:
            s2 = orig
        cpg = set(i for i in range(L - 1) if s2[i] == "C" and s2[i + 1] == "G")
        AN = anchors(s2)
        for k, tss in enumerate(AN):
            if tss - 1120 < 0 or tss + 1119 > L - 1:
                continue
            ukey = "%s|%s|%d" % (smp, r.query_name, k)
            nc = nm = 0
            for i in range(tss, min(tss + 119, L)):
                if i in cpg and i in cm:
                    nc += 1
                    if cm[i] >= MLCUT:
                        nm += 1
            slot = ap = Nn = zone = ""
            if UMAP and r.reference_name:
                try:
                    slot = r.reference_name.rsplit("_h", 1)[1]
                except IndexError:
                    slot = ""
                if slot in UMAP:
                    pos_s = tss + 4                                   # 닻 9-mer 의 가운데
                    orig_pos = (L - 1 - pos_s) if aflip else pos_s
                    st_pos = (L - 1 - orig_pos) if r.is_reverse else orig_pos
                    rp = q2r(r, st_pos)
                    if rp is not None:
                        starts, v = UMAP[slot]
                        j = bisect.bisect_right(starts, rp + 1) - 1
                        if j >= 0 and v[j][0] <= rp + 1 <= v[j][1]:
                            ap = v[j][2]
                            Nn = NSLOT[slot]
                            zone = ("원위" if ap <= 20 else ("근위" if ap >= Nn - 19 else "가운데"))
            ST["unit"] += 1
            ST["ap"] += (ap != "")
            fo.write(TB.join([ukey, r.query_name, "%d" % k, "%d" % nc, "%d" % nm,
                              ("%.4f" % (nm / nc)) if nc else "", str(slot), str(ap), str(Nn),
                              zone]) + NL)
    bam.close()
fo.close()
print(LABEL, dict(ST))
