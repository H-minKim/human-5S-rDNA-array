# -*- coding: utf-8 -*-
# Per-unit m6A and CpG of the placed Fiber-seq reads (v12)
#   뼈대 = v11 mkbam_ao_v11.py (BAM·참조 만들기) + meth_units_ao_v11.py (좌표 규칙) 를 그대로 옮김
#   바뀐 것 (환경변수로만)
#     D     : 작업폴더 (~/aog8v12fs2 · ~/aog8v12pb)     — v11 은 ~/aog8v10
#     CFG   : AOJV_FS · AOJV_PB                          — 이 cfg 의 SAM 만 읽는다
#     LABEL : UCSC_ONT_FIBERSEQ · PACBIO_FIBERSEQ        — v11 은 srcmap_<host>.tsv
#     OUT   : 새 폴더 (v11 폴더에 쓰지 않는다)
#     FINAL : 최종 유닛표 (in_final 열 표시용 · 열 8 = array_pos G8)
#     SMPS  : (선택) 쉼표로 개체 제한 — 파일럿
#     SAMD  : (선택) SAM 폴더 — 기본 D/sam
#   MM 블록 끝 플래그가 없으면 SAM 명세대로 '.' 와 같게 본다 (ft 출력 A+a,/T-a, · jasmine·ONT 는 '.' 명시라 영향 없음)
#   CpG : v11 코드 그대로 (C+m · 역방향 -1 · ML>=128) -> METH_READUNIT / METH_SITE  + 적중률(정본 C·G)
#   m6A : 추가
#     · A+a (원본 read 의 A) / T-a (원본 read 의 T = 반대가닥 A · PacBio 에만)
#     · 플래그 '.' = 목록에 없는 염기는 0 -> 분모 = 해당 염기 전부 (C+m '?' 는 목록만)
#     · 한 염기짜리 변형이라 CpG 의 -1 보정 없음
#     · mod_strand = 정본 기준 m6A 가 놓인 가닥 : '+' 면 정본 A 자리 · '-' 면 정본 T 자리
#         A+a : 정방향 read '+' · 역방향 read '-'     T-a : 정방향 '-' · 역방향 '+'
#     · 적중률 = 정본 염기가 기대 염기(A/T)와 같은 비율 (좌표 한 칸 틀리면 ~0.3 로 떨어짐)
#   ARRAY_ONLY·junction 모두 ao_g8cut 로 정본 창에 다시 정렬됐으므로 offset 0
import io
import os
import re
import sys
import gzip
import glob
import bisect
import subprocess
import collections

NL = chr(10)
TAB = chr(9)
H = os.path.expanduser
D = H(os.environ["D"])
CFG = os.environ["CFG"]
LABEL = os.environ["LABEL"]
OUT = H(os.environ["OUT"])
FINAL = os.environ.get("FINAL", "")
SMPS = set(x for x in os.environ.get("SMPS", "").split(",") if x)
SAMD = H(os.environ.get("SAMD", os.path.join(D, "sam")))   # SAM 폴더 (기본 D/sam · ft 재호출판은 ~/pbft_v12/sam)
MAPF = H("~/MAP_V11.tsv")
HOST = os.uname()[1]
CIG = re.compile(r"(\d+)([MIDNSHP=X])")
CMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
MLCUT = 128
OFF = 0
os.makedirs(os.path.join(OUT, "bam", "ref"), exist_ok=True)


def rc(s):
    return "".join(CMP.get(c, "N") for c in reversed(s))


def rd(p):
    return "".join(x.strip() for x in io.open(p, encoding="utf8").read().split(NL)[1:])


# ── 정본 서열 (mkbam_ao_v11 과 같은 방식)
SEQF = {}
for L in io.open(os.path.join(D, "g8in", "hap_seq_node_v7.tsv"), encoding="utf8"):
    f = L.rstrip(NL).split(TAB)
    if len(f) >= 3:
        SEQF[(f[0], f[1])] = f[2]
FLIP = {}
f = io.open(os.path.join(D, "SUM_G8.tsv"), encoding="utf8")
f.readline()
for L in f:
    r = L.rstrip(NL).split(TAB)
    FLIP[(r[0], r[1])] = (r[4] == "Y")

# ── 유닛 경계 (meth_units_ao_v11 과 같음)
U = collections.defaultdict(list)
f = io.open(MAPF, encoding="utf8")
ix = dict((n, i) for i, n in enumerate(f.readline().rstrip(NL).split(TAB)))
for L in f:
    r = L.rstrip(NL).split(TAB)
    if r[ix["seg_type"]] != "UNIT":
        continue
    U[(r[0], r[1])].append((int(r[ix["unit_start"]]), int(r[ix["unit_end"]]), int(r[ix["array_pos"]])))
for k in U:
    U[k].sort()
STARTS = dict((k, [x[0] for x in v]) for k, v in U.items())

# ── 최종 유닛표 (sample, read, slot, array_pos) OK
FIN = set()
if FINAL:
    for fp in glob.glob(H(FINAL)):
        for L in io.open(fp, encoding="utf8"):
            r = L.rstrip(NL).split(TAB)
            if len(r) > 13 and r[4] == CFG and r[13] == "OK":
                FIN.add((r[0], r[1], r[3], r[7]))                  # r[7]=array_pos(G8) · r[8]=G6 번호 (fix_relabel_only 주석)


# ── 1. BAM + 참조 (mkbam_ao_v11 그대로 · 출력 폴더만 새로)
def mkbam(sp, smp):
    recs = [L.rstrip(NL) for L in io.open(sp, encoding="utf8") if L.strip() and not L.startswith("@")]
    nzr = sum(1 for x in recs if (TAB + "ZR:Z:") in x)
    recs = [TAB.join(g for g in x.split(TAB) if not g.startswith("ZR:Z:")) for x in recs]
    used = sorted(set(x.split(TAB, 3)[2] for x in recs))
    hdr = ["@HD" + TAB + "VN:1.6" + TAB + "SO:unsorted"]
    refs = {}
    fa = io.open(os.path.join(OUT, "bam", "ref", "%s.%s.G8v3ref.fa" % (smp, CFG)), "w", encoding="utf8", newline=NL)
    for nm in used:
        s_, sl = nm.rsplit("_h", 1)
        if (s_, sl) not in SEQF:
            return None, "NO_REFSEQ"
        a = rd(os.path.join(D, "cmp4d_v7", SEQF[(s_, sl)])).upper()
        a = rc(a) if FLIP.get((s_, sl)) else a
        refs[nm] = a
        hdr.append("@SQ" + TAB + "SN:" + nm + TAB + "LN:%d" % len(a))
        fa.write(">" + nm + NL + a + NL)
    fa.close()
    hdr.append("@CO" + TAB + "Fiber-seq placed reads (AOJV) re-aligned to haplotype canonical (G8v3 orientation) by ao_g8cut_v12; SEQ/QUAL/MM/ML untouched")
    tmp = os.path.join(OUT, "bam", "tmp_%s.%s.sam" % (smp, CFG))
    io.open(tmp, "w", encoding="utf8", newline=NL).write(NL.join(hdr + recs) + NL)
    out = os.path.join(OUT, "bam", "%s.%s.%s.G8v3.bam" % (smp, CFG, LABEL))
    r1 = subprocess.run(["samtools", "sort", "-@2", "-o", out, tmp], capture_output=True, text=True)
    r2 = subprocess.run(["samtools", "index", out], capture_output=True, text=True) if r1.returncode == 0 else r1
    os.remove(tmp)
    if r1.returncode or r2.returncode:
        return None, "SORT_FAIL"
    return (out, refs), "ZR%d rec%d" % (nzr, len(recs))


def blocks(seq, flag, mm, ml):
    """(블록이름, [(SEQ좌표, 값)]) — MM 은 원본 read 방향. '.' 블록은 목록 밖 염기를 0 으로 채움."""
    orig = rc(seq) if (flag & 16) else seq
    vals = [int(x) for x in ml.split(",") if x != ""]
    res = []
    base = 0
    for blk in [p for p in mm.split(";") if p]:
        f2 = blk.split(",")
        deltas = [int(x) for x in f2[1:]]
        head = f2[0]
        nm = head.rstrip(".?")
        b0 = nm[0]
        if nm in ("C+m", "A+a", "T-a"):
            idx = [i for i, c in enumerate(orig) if c == b0]
            got = {}
            j = -1
            for t, d in enumerate(deltas):
                j += d + 1
                if j >= len(idx):
                    break
                got[j] = vals[base + t] if base + t < len(vals) else 0
            if nm == "C+m":
                items = sorted(got.items())                       # v11 : 목록만
            elif head.endswith(".") or head[-1] not in ".?":
                items = [(j, got.get(j, 0)) for j in range(len(idx))]   # 분모 = 염기 전부
            else:
                items = sorted(got.items())
            out = []
            for j, v in items:
                op = idx[j]
                sp = (len(seq) - 1 - op) if (flag & 16) else op
                out.append((sp, v))
            res.append((nm, head[-1] if head[-1] in ".?" else "", out))
        base += len(deltas)
    return res


RU = gzip.open(os.path.join(OUT, "METH_READUNIT_%s.tsv.gz" % HOST), "wt", encoding="utf8", newline=NL)
RU.write(TAB.join(["sample", "data_source", "cfg", "kind", "slot", "array_pos", "read_id",
                   "strand", "n_cpg", "n_meth", "beta", "in_final", "calls"]) + NL)
SI = gzip.open(os.path.join(OUT, "METH_SITE_%s.tsv.gz" % HOST), "wt", encoding="utf8", newline=NL)
SI.write(TAB.join(["sample", "data_source", "slot", "array_pos", "unit_pos", "n", "n_meth", "beta"]) + NL)
AU = gzip.open(os.path.join(OUT, "M6A_READUNIT_%s.tsv.gz" % HOST), "wt", encoding="utf8", newline=NL)
AU.write(TAB.join(["sample", "data_source", "cfg", "kind", "slot", "array_pos", "read_id", "strand",
                   "mod_strand", "n_a", "n_m6a", "frac", "in_final", "calls"]) + NL)
AS = gzip.open(os.path.join(OUT, "M6A_SITE_%s.tsv.gz" % HOST), "wt", encoding="utf8", newline=NL)
AS.write(TAB.join(["sample", "data_source", "slot", "array_pos", "unit_pos", "mod_strand", "n", "n_m6a", "frac"]) + NL)
LOG = io.open(os.path.join(OUT, "LOG_%s.tsv" % HOST), "w", encoding="utf8", newline=NL)
LOG.write(TAB.join(["host", "cfg", "sample", "bam", "read", "no_mm", "mm_empty", "no_m6a_block", "read_unit", "in_final",
                    "cpg_calls", "cpg_hit", "cpg_hit_noshift",
                    "m6a_calls_plus", "m6a_hit_plus", "m6a_calls_minus", "m6a_hit_minus", "m6a_hit_shift1"]) + NL)

TOT = collections.Counter()
for sp_ in sorted(glob.glob(os.path.join(SAMD, "*.%s.sam" % CFG))):
    smp = os.path.basename(sp_).split(".")[0]
    if SMPS and smp not in SMPS:
        continue
    got, msg = mkbam(sp_, smp)
    if got is None:
        LOG.write(TAB.join([HOST, CFG, smp, msg]) + NL)
        continue
    b, refs = got
    T = collections.Counter()
    SITE = collections.Counter()
    SITEM = collections.Counter()
    ASITE = collections.Counter()
    ASITEM = collections.Counter()
    pr = subprocess.Popen(["samtools", "view", b], stdout=subprocess.PIPE, text=True)
    for L in pr.stdout:
        r = L.rstrip(NL).split(TAB)
        if len(r) < 11:
            continue
        flag = int(r[1])
        if flag & 0x900 or r[2] == "*":
            continue
        sn, pos, cig, seq = r[2], int(r[3]), r[5], r[9]
        sl = sn.split("_h")[-1]
        key = (smp, sl)
        if key not in U:
            T["no_map"] += 1
            continue
        mm = ml = None
        for x in r[11:]:
            if x.startswith("MM:Z:"):
                mm = x[5:]
            elif x.startswith("ML:B:C,"):
                ml = x[7:]
        T["read"] += 1
        if not mm:
            T["no_mm"] += 1
            continue
        if not ml:
            T["mm_empty"] += 1                    # MM 블록만 있고 호출 0 (ML 비어 있음) -> 0 으로 채우지 않고 뺀다
            continue
        if "A+a" not in mm and "T-a" not in mm:
            T["no_m6a_block"] += 1                # 원자료에 m6A 호출 없음 (CpG 는 계속)
        ref = refs[sn]
        m = {}
        rp = pos - 1
        q = 0
        for g in CIG.finditer(cig):
            k, op = int(g.group(1)), g.group(2)
            if op == "S":
                q += k
            elif op == "H":
                pass
            elif op in "M=X":
                for t in range(k):
                    m[q + t] = rp + t
                rp += k
                q += k
            elif op == "I":
                q += k
            elif op in "DN":
                rp += k
        st = "-" if (flag & 16) else "+"
        S = STARTS[key]
        per = collections.defaultdict(list)
        aper = collections.defaultdict(list)
        for nm, fl, calls in blocks(seq, flag, mm, ml):
            if nm == "C+m":
                for sp, v in calls:
                    if sp not in m:
                        continue
                    cp = m[sp] - OFF
                    # 검정 : 보정 전 자리
                    if 0 <= cp < len(ref) and ref[cp] == ("G" if flag & 16 else "C"):
                        T["cpg_hit_noshift"] += 1
                    if flag & 16:
                        cp -= 1
                    T["cpg_calls"] += 1
                    if 0 <= cp < len(ref) - 1 and ref[cp:cp + 2] == "CG":
                        T["cpg_hit"] += 1
                    p1 = cp + 1
                    if p1 < 1:
                        continue
                    i = bisect.bisect_right(S, p1) - 1
                    if i < 0:
                        continue
                    s0, e0, ap = U[key][i]
                    if p1 > e0:
                        continue
                    per[ap].append((p1 - s0 + 1, v))
            else:
                ms = "+" if ((nm == "A+a") != bool(flag & 16)) else "-"
                want = "A" if ms == "+" else "T"
                for sp, v in calls:
                    if sp not in m:
                        continue
                    cp = m[sp] - OFF
                    T["m6a_calls" + ms] += 1
                    if 0 <= cp < len(ref) and ref[cp] == want:
                        T["m6a_hit" + ms] += 1
                    if 0 <= cp + 1 < len(ref) and ref[cp + 1] == want:
                        T["m6a_hit_shift1"] += 1
                    p1 = cp + 1
                    if p1 < 1:
                        continue
                    i = bisect.bisect_right(S, p1) - 1
                    if i < 0:
                        continue
                    s0, e0, ap = U[key][i]
                    if p1 > e0:
                        continue
                    aper[(ap, ms)].append((p1 - s0 + 1, v))
        for ap in sorted(set(list(per.keys()) + [a for a, _ in aper.keys()])):
            T["read_unit"] += 1
            inf = "Y" if (smp, r[0], sl, str(ap)) in FIN else "N"
            if inf == "Y":
                T["in_final"] += 1
            lst = sorted(per.get(ap, []))
            if lst:
                nm_ = sum(1 for _, v in lst if v >= MLCUT)
                RU.write(TAB.join([smp, LABEL, CFG, "aojv", sl, str(ap), r[0], st, str(len(lst)), str(nm_),
                                   "%.4f" % (nm_ / float(len(lst))), inf,
                                   ",".join("%d:%d" % (a, v) for a, v in lst)]) + NL)
                for a, v in lst:
                    SITE[(sl, ap, a)] += 1
                    if v >= MLCUT:
                        SITEM[(sl, ap, a)] += 1
            for ms in "+-":
                al = sorted(aper.get((ap, ms), []))
                if not al:
                    continue
                na = sum(1 for _, v in al if v >= MLCUT)
                AU.write(TAB.join([smp, LABEL, CFG, "aojv", sl, str(ap), r[0], st, ms, str(len(al)), str(na),
                                   "%.4f" % (na / float(len(al))), inf,
                                   ",".join("%d:%d" % (a, v) for a, v in al)]) + NL)
                for a, v in al:
                    ASITE[(sl, ap, a, ms)] += 1
                    if v >= MLCUT:
                        ASITEM[(sl, ap, a, ms)] += 1
    pr.stdout.close()
    pr.wait()
    for k in sorted(SITE, key=lambda x: (int(x[0]), x[1], x[2])):
        n = SITE[k]
        nm_ = SITEM.get(k, 0)
        SI.write(TAB.join([smp, LABEL, k[0], str(k[1]), str(k[2]), str(n), str(nm_), "%.4f" % (nm_ / float(n))]) + NL)
    for k in sorted(ASITE, key=lambda x: (int(x[0]), x[1], x[2], x[3])):
        n = ASITE[k]
        na = ASITEM.get(k, 0)
        AS.write(TAB.join([smp, LABEL, k[0], str(k[1]), str(k[2]), k[3], str(n), str(na), "%.4f" % (na / float(n))]) + NL)

    def rt(a, b):
        return "%.4f" % (T[a] / float(T[b])) if T[b] else "NA"
    tot_m6a = T["m6a_calls+"] + T["m6a_calls-"]
    LOG.write(TAB.join([HOST, CFG, smp, os.path.basename(b), str(T["read"]), str(T["no_mm"]), str(T["mm_empty"]), str(T["no_m6a_block"]), str(T["read_unit"]),
                        str(T["in_final"]), str(T["cpg_calls"]), rt("cpg_hit", "cpg_calls"), rt("cpg_hit_noshift", "cpg_calls"),
                        str(T["m6a_calls+"]), rt("m6a_hit+", "m6a_calls+"), str(T["m6a_calls-"]), rt("m6a_hit-", "m6a_calls-"),
                        "%.4f" % (T["m6a_hit_shift1"] / float(tot_m6a)) if tot_m6a else "NA"]) + NL)
    LOG.flush()
    TOT.update(T)
    print("  %s %s read %d · read-unit %d · in_final %d" % (CFG, smp, T["read"], T["read_unit"], T["in_final"]), flush=True)

for fh in (RU, SI, AU, AS, LOG):
    fh.close()
print("DONE %s %s" % (HOST, CFG))
