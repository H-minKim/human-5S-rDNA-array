# -*- coding: utf-8 -*-
# read x 유닛 메틸을 coding(5S 유전자 cmsearch 구간) / 나머지(스페이서)로 나눠 다시 센다.
#
#   기존 03_reads/*.read_unit.tsv 는 유닛 전체 CpG 만 합쳐 두어 coding 만 따로 볼 수 없다.
#   그래서 배치 BAM(04_bam)의 MM/ML 을 직접 읽는다 (로컬에 samtools/pysam/modkit 없음 → 순수 파이썬 BAM 파서).
#   호출 규칙은 units_matrix_g8.py 와 같게 맞춘다 :
#     보조/이차 제외(0x900) · mapq 로 거르지 않음 · 참조 CpG 자리만(정방향 read = 참조 C 다음 G ·
#     역방향 read = 참조 G 앞 C, 자리는 G) · 5mC 호출 = mod_qual >= 0.5 (ML >= 128) ·
#     유닛 = g8map chrom_start..chrom_end (1기반 양끝 포함)
#   유전자 = cmsearch RF00001 포함(!) 히트 중 그 유닛 구간과 겹치는 것 (1기반 양끝 포함)
#   참조 서열 = 01_canonical/<species>/<array>.fa (flank 20kb · strand=- 이면 역상보 → 되돌려 씀)
#
#   검산 : 유닛 전체 (n_calls, n_5mC) 가 기존 read_unit.tsv 와 같아야 한다. 다르면 멈춘다(판정은 로그).
#   인자 : [종키 ...]  (없으면 전 종)  ·  --pilot 이면 검산만
#   산출 : 10_other_species/tables/OTHERSP_READUNIT_CODING.tsv.gz · 로그 X11_CODING_BUILD.txt
import io, os, sys, gzip, glob, struct, re, collections
import numpy as np
import pandas as pd

import os as _os_need, sys as _sys_need


def _need(k):
    """No defaults are set: if an environment variable is missing the script stops instead of silently reading the wrong location."""
    v = _os_need.environ.get(k)
    if not v:
        _sys_need.exit("environment variable %s is not set (the deposited copy has no default; "
                       "see README.md, section 'How to run it')" % k)
    return v


B = _need("RDNA_DATA_ROOT")
OUT = _need("RDNA_PREV_OUT_DIR")
PILOT = "--pilot" in sys.argv
SPS = [a for a in sys.argv[1:] if not a.startswith("--")] or [
    "gorilla_mGorGor1", "bonobo_mPanPan1", "chimp_mPanTro3", "orang_bornean_mPonPyg2",
    "orang_sumatran_mPonAbe1", "siamang_mSymSyn1", "zfinch_v03_dip", "arabidopsis_ColCC", "zebrafish_GRCz12tu"]
LOG = io.open(OUT + ("/X11_CODING_PILOT.txt" if PILOT else "/X11_CODING_BUILD.txt"), "a", encoding="utf8", newline="\n")


def P(s=""):
    try:
        print(s)
    except UnicodeEncodeError:
        print(s.encode("ascii", "replace").decode())
    LOG.write(s + "\n")
    LOG.flush()


COMP = bytes.maketrans(b"ACGTNacgtn", b"TGCANtgcan")
NIBCOMP = np.full(16, 15, dtype=np.uint8)
NIBCOMP[[1, 2, 4, 8]] = [8, 4, 2, 1]


def read_fa(p):
    L = io.open(p, encoding="utf8").read().split("\n")
    return L[0], "".join(x.strip() for x in L[1:]).upper()


def load_species(sp):
    """array 마다 : 정방향 참조(numpy uint8) · 0기반 시작 좌표 · 유닛 구간 · 유전자 구간"""
    A = {}
    tbl = collections.defaultdict(list)
    for L in io.open(B + "/cms_tbl/%s.tbl" % sp, encoding="utf8"):
        if L.startswith("#"):
            continue
        x = L.split()
        if len(x) < 17 or x[16] != "!":
            continue
        a, b = int(x[7]), int(x[8])
        tbl[x[0]].append((min(a, b), max(a, b)))
    for c in tbl:
        tbl[c].sort()
    for f in sorted(glob.glob(B + "/02_consensus_map/%s/*.g8map.tsv" % sp)):
        if not re.search(r"_[0-9]+_[0-9]+[.]g8map[.]tsv$", f):   # 대체 지도(.k2.g8map.tsv 등) 제외 — 본 지도만
            continue
        g = pd.read_csv(f, sep="\t")
        aid = g.array_id.iloc[0]
        fa = B + "/01_canonical/%s/%s.fa" % (sp, aid)
        if not os.path.exists(fa):
            continue
        hd, s = read_fa(fa)
        m = re.search(r"flank(\d+)_array_flank(\d+) strand=([+-]) chrom=(\S+) array=(\d+)-(\d+)", hd)
        f5, f3, strand, chrom, a0, a1 = int(m.group(1)), int(m.group(2)), m.group(3), m.group(4), int(m.group(5)), int(m.group(6))
        if strand == "-":
            s = s.encode().translate(COMP)[::-1].decode()
        # 정방향 참조 첫 글자의 1기반 좌표
        lo1 = a0 - (f5 if strand == "+" else f3)
        if len(s) != (a1 - a0 + 1) + f5 + f3:
            P("   ! 길이 불일치 %s len=%d 기대=%d" % (aid, len(s), (a1 - a0 + 1) + f5 + f3))
        us = np.minimum(g.chrom_start.values, g.chrom_end.values)
        ue = np.maximum(g.chrom_start.values, g.chrom_end.values)
        o = np.argsort(us)
        genes = []
        for i in o:
            genes.append([(max(a, us[i]), min(b, ue[i])) for a, b in tbl[chrom] if a <= ue[i] and b >= us[i]])
        A[aid] = dict(chrom=chrom, seq=np.frombuffer(s.encode(), dtype=np.uint8), lo1=lo1,
                      us=us[o], ue=ue[o], pos=g.array_pos.values[o], genes=genes)
    return A


def bam_records(path):
    f = gzip.open(path, "rb")
    assert f.read(4) == b"BAM\1"
    lt, = struct.unpack("<i", f.read(4)); f.read(lt)
    nref, = struct.unpack("<i", f.read(4))
    refs = []
    for _ in range(nref):
        ln, = struct.unpack("<i", f.read(4))
        refs.append(f.read(ln)[:-1].decode()); f.read(4)
    while True:
        h = f.read(4)
        if len(h) < 4:
            break
        bs, = struct.unpack("<i", h)
        yield refs, f.read(bs)


def parse_tags(buf, p):
    T = {}
    n = len(buf)
    while p < n:
        tag = buf[p:p + 2].decode(); ty = chr(buf[p + 2]); p += 3
        if ty in "AcC":
            p += 1
        elif ty in "sS":
            p += 2
        elif ty in "iIf":
            p += 4
        elif ty in "ZH":
            e = buf.index(b"\0", p)
            if tag == "MM" or tag == "Mm":
                T["MM"] = buf[p:e].decode()
            p = e + 1
        elif ty == "B":
            st = chr(buf[p]); cnt, = struct.unpack("<i", buf[p + 1:p + 5]); p += 5
            w = {"c": 1, "C": 1, "s": 2, "S": 2, "i": 4, "I": 4, "f": 4}[st]
            if tag == "ML" or tag == "Ml":
                T["ML"] = np.frombuffer(buf[p:p + cnt * w], dtype=np.uint8)
            p += cnt * w
        else:
            raise ValueError("tag type " + ty)
    return T


STAT = collections.Counter()


def process(sp, sample, bam, A, acc):
    bychrom = collections.defaultdict(list)
    for aid, a in A.items():
        bychrom[a["chrom"]].append(aid)
    for refs, r in bam_records(bam):
        refid, pos0, lrn, mapq, _bin, ncig, flag, lseq = struct.unpack("<iiBBHHHi", r[:20])
        STAT["rec"] += 1
        if refid < 0 or flag & 0x904:
            continue
        chrom = refs[refid]
        if chrom not in bychrom:
            continue
        rid = r[32:32 + lrn - 1].decode()
        p = 32 + lrn
        cig = np.frombuffer(r[p:p + 4 * ncig], dtype="<u4"); p += 4 * ncig
        sq = r[p:p + (lseq + 1) // 2]; p += (lseq + 1) // 2
        p += lseq
        T = parse_tags(r, p)
        if "MM" not in T or "ML" not in T:
            STAT["no_mm"] += 1
            continue
        ops, lens = cig & 15, (cig >> 4).astype(np.int64)
        if (ops == 5).any():
            STAT["hardclip"] += 1
            continue
        # 질의 -> 참조 (0기반)
        q2r = np.full(lseq, -1, dtype=np.int64)
        q, rr = 0, pos0
        for o, n in zip(ops.tolist(), lens.tolist()):
            if o in (0, 7, 8):
                q2r[q:q + n] = np.arange(rr, rr + n); q += n; rr += n
            elif o in (1, 4):
                q += n
            elif o in (2, 3):
                rr += n
        b = np.frombuffer(sq, dtype=np.uint8)
        nib = np.empty(2 * len(b), dtype=np.uint8); nib[0::2] = b >> 4; nib[1::2] = b & 15
        codes = nib[:lseq]                       # 1=A 2=C 4=G 8=T 15=N
        rev = bool(flag & 16)
        # 원래 read 방향에서의 C 위치
        ocodes = NIBCOMP[codes][::-1] if rev else codes
        cpos = np.flatnonzero(ocodes == 2)
        ml = T["ML"]; mi = 0
        m_idx, m_ml = None, None
        for ent in T["MM"].rstrip(";").split(";"):
            if not ent:
                continue
            parts = ent.split(",")
            head = parts[0]
            mode = "?"
            if head[-1] in "?.":
                mode = head[-1]; head = head[:-1]
            base, strand_, mc = head[0], head[1], head[2:]
            d = np.array([int(x) for x in parts[1:]], dtype=np.int64)
            nc = len(mc) if not mc.isdigit() else 1
            vals = ml[mi:mi + nc * len(d)].reshape(-1, nc) if len(d) else np.zeros((0, nc), np.uint8)
            mi += nc * len(d)
            if base != "C" or strand_ != "+" or "m" not in mc:
                continue
            STAT["mode" + mode] += 1
            if len(d) == 0:
                continue
            k = np.cumsum(d + 1) - 1
            if k[-1] >= len(cpos):
                STAT["mm_overflow"] += 1
                break
            m_idx, m_ml = cpos[k], vals[:, mc.index("m")]
            if mode == ".":                          # 생략된 C = 비메틸로 간주(modkit inferred 와 같게)
                full = np.zeros(len(cpos), dtype=np.uint8)
                full[k] = m_ml
                m_idx, m_ml = cpos, full
        if m_idx is None:
            continue
        sidx = (lseq - 1 - m_idx) if rev else m_idx
        ref0 = q2r[sidx]
        ok = ref0 >= 0
        ref0, mlv = ref0[ok], m_ml[ok]
        for aid in bychrom[chrom]:
            a = A[aid]
            k0 = ref0 - (a["lo1"] - 1)            # 정방향 참조 배열 안 0기반 자리
            inr = (k0 >= 1) & (k0 < len(a["seq"]) - 1)
            if not inr.any():
                continue
            kk, vv, rp = k0[inr], mlv[inr], ref0[inr] + 1
            s = a["seq"]
            if rev:
                cg = (s[kk] == 71) & (s[kk - 1] == 67)
            else:
                cg = (s[kk] == 67) & (s[kk + 1] == 71)
            kk, vv, rp = kk[cg], vv[cg], rp[cg]
            ui = np.searchsorted(a["us"], rp, side="right") - 1
            good = (ui >= 0) & (rp <= a["ue"][np.clip(ui, 0, None)])
            ui, vv, rp = ui[good], vv[good], rp[good]
            if len(ui) == 0:
                continue
            ingene = np.zeros(len(rp), dtype=bool)
            for u in np.unique(ui):
                for g0, g1 in a["genes"][u]:
                    ingene |= (ui == u) & (rp >= g0) & (rp <= g1)
            call = vv >= 128
            for u in np.unique(ui):
                w = ui == u
                wg = w & ingene
                acc.append((sp, sample, rid, mapq, aid, int(a["pos"][u]), len(a["genes"][u]),
                            int(w.sum()), int((call & w).sum()), int(wg.sum()), int((call & wg).sum())))


COLS = ["species_key", "sample", "read_id", "mapq", "array_id", "array_pos", "n_gene_hits",
        "n_calls", "n_5mC", "n_calls_gene", "n_5mC_gene"]
allrows = []
for sp in SPS:
    A = load_species(sp)
    bams = sorted(glob.glob(B + "/04_bam/%s/*.5S.bam" % sp))
    if PILOT:
        bams = bams[:1]
    ref = pd.concat([pd.read_csv(f, sep="\t") for f in glob.glob(B + "/03_reads/%s/*.read_unit.tsv*" % sp)])
    ref = ref[ref.context == "CG"]
    for bam in bams:
        smp = os.path.basename(bam).replace(".5S.bam", "")
        acc = []
        STAT.clear()
        process(sp, smp, bam, A, acc)
        D = pd.DataFrame(acc, columns=COLS)
        # 같은 read 가 array 두 번 이상 걸리면 합친다 (units_matrix 도 read·유닛 키로 합산)
        D = D.groupby(["species_key", "sample", "read_id", "mapq", "array_id", "array_pos", "n_gene_hits"],
                      as_index=False).sum()
        rs = ref[ref["sample"].astype(str).isin([smp, smp.replace(".g12", "")])]
        M = rs.merge(D, on=["read_id", "array_id", "array_pos"], how="outer", suffixes=("_old", ""), indicator=True)
        both = M[M._merge == "both"]
        same = ((both.n_calls_old == both.n_calls) & (both.n_5mC_old == both.n_5mC)).mean() if len(both) else np.nan
        P("%s %s : 기록 %d · 행 새 %d / 옛 %d · 공통 %d · 옛만 %d · 새만 %d · (n_calls,n_5mC) 완전일치 %.4f · 통계 %s"
          % (sp, smp, STAT["rec"], len(D), len(rs), len(both), (M._merge == "left_only").sum(),
             (M._merge == "right_only").sum(), same, dict(STAT)))
        if len(both):
            dd = both[(both.n_calls_old != both.n_calls) | (both.n_5mC_old != both.n_5mC)]
            if len(dd):
                P("   불일치 예 : " + dd[["array_pos", "n_calls_old", "n_calls", "n_5mC_old", "n_5mC"]].head(5).to_string().replace("\n", " | "))
            P("   유전자 CpG 수 중앙 %.0f · 유전자 없는 유닛 비율 %.3f"
              % (D[D.n_gene_hits > 0].n_calls_gene.median(), (D.n_gene_hits == 0).mean()))
        allrows.append(D)
if not PILOT:
    R = pd.concat(allrows, ignore_index=True)
    R.to_csv(OUT + "/tables/OTHERSP_READUNIT_CODING.tsv.gz", sep="\t", index=False, compression="gzip")
    P("저장 : tables/OTHERSP_READUNIT_CODING.tsv.gz (%d 행)" % len(R))
LOG.close()
