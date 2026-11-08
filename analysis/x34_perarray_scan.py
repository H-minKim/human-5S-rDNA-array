# -*- coding: utf-8 -*-
# array 마다 두 가지를 BAM 한 번 훑기로 같이 모은다 (시아망 제외 53 array)
#   (1) 유닛 안 위치별 메틸 프로파일 — 분자 등급(ALL / HIGH β>=0.70 / LOW β<0.35) 별   → 대주제 G (사람 예시 패널 e 형식)
#   (2) array 경계에서의 거리별 메틸 — 안/밖 · 게놈 좌표 쪽(낮은 좌표 / 높은 좌표)   → 대주제 E 지형 지도
# Calling rules : 5mC = ML>=128, flag&0x904 excluded, mapq unfiltered, reference CpG positions only
#       잘 덮인 read-유닛만(n_calls >= array 중앙값 x 0.5 · 중앙값은 모든 read 에서) · 보노보 guppy 6.3.8 만
#   방향 : 유닛 안 좌표는 g8 틀(유전자가 유닛 끝) · array 지형은 게놈 좌표 로 모은다
#           (열 이름이 아니라 실제 좌표로 — 낮은 좌표 쪽 / 높은 좌표 쪽)
#   산출 : tables/PERARRAY_WITHIN_UNIT.tsv.gz · tables/PERARRAY_LANDSCAPE.tsv.gz · tables/PERARRAY_META.tsv
#   ※ 그림은 x35(G) · x36(E) 가 이 표만 읽어 그린다 — 그림을 고칠 때 BAM 을 다시 읽지 않게.
import io, os, sys, glob, re, struct, collections
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


NEW = _need("RDNA_OUT_DIR")
B = _need("RDNA_DATA_ROOT")
OLD = _need("RDNA_PREV_OUT_DIR")
sys.path.insert(0, NEW + "/code")
import arch  # noqa
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
LOG = io.open(NEW + "/logs/X34_PERARRAY.txt", "w", encoding="utf8", newline="\n")
REASON = "array 마다 유닛 안 프로파일 + 경계 거리별 메틸 (2026-10-01)"
CUT, LOWC, HIGHC, FLMAX = 128, 0.35, 0.70, 20000


def P(s=""):
    print(s)
    LOG.write(s + "\n")
    LOG.flush()


src = io.open(OLD + "/code/x11_coding_readunit.py", encoding="utf8").read()
L11 = 'LOG = io.open(OUT + ("/X11_CODING_PILOT.txt" if PILOT else "/X11_CODING_BUILD.txt"), "a", encoding="utf8", newline="\\n")'
head = src.split("STAT = collections.Counter()")[0].replace(L11, "LOG = io.StringIO()")
G = {"__name__": "x11hdr"}
exec(compile(head, OLD + "/code/x11_coding_readunit.py", "exec"), G)
NIBCOMP, load_species, bam_records, parse_tags, read_fa = (
    G["NIBCOMP"], G["load_species"], G["bam_records"], G["parse_tags"], G["read_fa"])

A = pd.read_csv(NEW + "/tables/FIGDATA_A_array_summary.tsv", sep="\t")
LAB = dict(zip(A.array_id, A.이름표))
KO = dict(zip(A.array_id, A.종))
RU = pd.read_csv(NEW + "/tables/OTHERSP_READUNIT_CODING.tsv.gz", sep="\t", usecols=["array_id", "n_calls"])
MED = RU.groupby("array_id").n_calls.median().to_dict()
del RU
BC = pd.read_csv(NEW + "/tables/BONOBO_READ_BASECALLER.tsv", sep="\t")
B637 = set(BC.loc[BC.basecaller.astype(str).str.contains("6.3.7"), "read_id"])
SPK = [("고릴라", "gorilla_mGorGor1"), ("보노보", "bonobo_mPanPan1"), ("침팬지", "chimp_mPanTro3"),
       ("보르네오오랑우탄", "orang_bornean_mPonPyg2"), ("수마트라오랑우탄", "orang_sumatran_mPonAbe1"),
       ("애기장대", "arabidopsis_ColCC"), ("금화조", "zfinch_v03_dip"), ("제브라피시", "zebrafish_GRCz12tu")]


def strand_of(sp):
    out = {}
    for f in sorted(glob.glob(B + "/02_consensus_map/%s/*.g8map.tsv" % sp)):
        if not re.search(r"_[0-9]+_[0-9]+[.]g8map[.]tsv$", f.replace("\\", "/")):
            continue
        g = pd.read_csv(f, sep="\t", usecols=["array_id", "strand"])
        out[g.array_id.iloc[0]] = g.strand.iloc[0]
    return out


def win_edges(inside, imax):
    if inside:
        e = list(np.arange(0, 2000, 50)) + list(np.arange(2000, 10000, 200)) + list(np.arange(10000, max(imax, 10001), 1000))
    else:
        e = list(np.arange(0, 2000, 100)) + list(np.arange(2000, 10000, 500)) + list(np.arange(10000, FLMAX + 1, 5000))
    return np.array(e, dtype=float)


WU, LS, MET = [], [], []
for ko, sp in SPK:
    Asp = load_species(sp)
    SD = strand_of(sp)
    keep = [a for a in Asp if a in LAB and KO.get(a) == ko]
    bychrom = collections.defaultdict(list)
    for a in keep:
        bychrom[Asp[a]["chrom"]].append(a)
    ST, LW = {}, {}
    for a in keep:
        v = Asp[a]
        hd, _ = read_fa(B + "/01_canonical/%s/%s.fa" % (sp, a))
        m = re.search(r"array=(\d+)-(\d+)", hd)
        a0, a1 = int(m.group(1)), int(m.group(2))
        ul = int(np.median(v["ue"] - v["us"] + 1))
        imax = int((21) * ul)
        ST[a] = dict(a0=a0, a1=a1, ul=ul, imax=imax, strand=SD[a],
                     acc={k: np.zeros((2, ul + 2)) for k in ("ALL", "HIGH", "LOW")},
                     cnt=collections.Counter())
        LW[a] = {(side, ins): [np.zeros(len(win_edges(ins, imax))), np.zeros(len(win_edges(ins, imax)))]
                 for side in ("낮은 좌표", "높은 좌표") for ins in (True, False)}
    for bam in sorted(glob.glob(B + "/04_bam/%s/*.5S.bam" % sp)):
        for refs, r in bam_records(bam):
            refid, pos0, lrn, mapq, _bin, ncig, flag, lseq = struct.unpack("<iiBBHHHi", r[:20])
            if refid < 0 or (flag & 0x904) or refs[refid] not in bychrom:
                continue
            rid = r[32:32 + lrn - 1].decode()
            if sp == "bonobo_mPanPan1" and rid in B637:
                continue
            p = 32 + lrn
            cig = np.frombuffer(r[p:p + 4 * ncig], dtype="<u4"); p += 4 * ncig
            sq = r[p:p + (lseq + 1) // 2]; p += (lseq + 1) // 2 + lseq
            T = parse_tags(r, p)
            if "MM" not in T or "ML" not in T:
                continue
            ops, lens = cig & 15, (cig >> 4).astype(np.int64)
            if (ops == 5).any():
                continue
            q2r = np.full(lseq, -1, dtype=np.int64)
            q, rr = 0, pos0
            for o, n in zip(ops.tolist(), lens.tolist()):
                if o in (0, 7, 8):
                    q2r[q:q + n] = np.arange(rr, rr + n); q += n; rr += n
                elif o in (1, 4):
                    q += n
                elif o in (2, 3):
                    rr += n
            bb = np.frombuffer(sq, dtype=np.uint8)
            nib = np.empty(2 * len(bb), dtype=np.uint8); nib[0::2] = bb >> 4; nib[1::2] = bb & 15
            codes = nib[:lseq]
            rev = bool(flag & 16)
            ocodes = NIBCOMP[codes][::-1] if rev else codes
            cpos = np.flatnonzero(ocodes == 2)
            ml = T["ML"]; mi = 0
            m_idx = m_ml = None
            for ent in T["MM"].rstrip(";").split(";"):
                if not ent:
                    continue
                parts = ent.split(","); hd2 = parts[0]; mode = "?"
                if hd2[-1] in "?.":
                    mode = hd2[-1]; hd2 = hd2[:-1]
                base, st_, mc = hd2[0], hd2[1], hd2[2:]
                d = np.array([int(x) for x in parts[1:]], dtype=np.int64)
                nc = len(mc) if not mc.isdigit() else 1
                vals = ml[mi:mi + nc * len(d)].reshape(-1, nc) if len(d) else np.zeros((0, nc), np.uint8)
                mi += nc * len(d)
                if base != "C" or st_ != "+" or "m" not in mc or len(d) == 0:
                    continue
                k = np.cumsum(d + 1) - 1
                if k[-1] >= len(cpos):
                    break
                m_idx, m_ml = cpos[k], vals[:, mc.index("m")]
                if mode == ".":
                    full = np.zeros(len(cpos), dtype=np.uint8); full[k] = m_ml
                    m_idx, m_ml = cpos, full
            if m_idx is None:
                continue
            sidx = (lseq - 1 - m_idx) if rev else m_idx
            ref0 = q2r[sidx]
            ok = ref0 >= 0
            ref0, mlv = ref0[ok], m_ml[ok]
            for aid in bychrom[refs[refid]]:
                v, s_ = Asp[aid], Asp[aid]["seq"]
                k0 = ref0 - (v["lo1"] - 1)
                inr = (k0 >= 1) & (k0 < len(s_) - 1)
                if not inr.any():
                    continue
                kk, vv, rp = k0[inr], mlv[inr], ref0[inr] + 1
                cg = (s_[kk] == 71) & (s_[kk - 1] == 67) if rev else (s_[kk] == 67) & (s_[kk + 1] == 71)
                vv, rp = vv[cg], rp[cg]
                if not len(rp):
                    continue
                call = (vv >= CUT).astype(float)
                st = ST[aid]
                # ---------- (2) 경계 거리별 (게놈 좌표 쪽)
                a0, a1 = st["a0"], st["a1"]
                ins_mask = (rp >= a0) & (rp <= a1)
                for side, d_in, d_out in (("낮은 좌표", rp - a0, a0 - rp), ("높은 좌표", a1 - rp, rp - a1)):
                    for inside, dd in ((True, d_in), (False, d_out)):
                        lim = st["imax"] if inside else FLMAX
                        m = (dd >= 0) & (dd <= lim) & (ins_mask if inside else ~ins_mask)
                        if not m.any():
                            continue
                        e = win_edges(inside, st["imax"])
                        bi = np.clip(np.searchsorted(e, dd[m], side="right") - 1, 0, len(e) - 1)
                        np.add.at(LW[aid][(side, inside)][0], bi, 1.0)
                        np.add.at(LW[aid][(side, inside)][1], bi, call[m])
                # ---------- (1) 유닛 안 프로파일
                ui = np.searchsorted(v["us"], rp, side="right") - 1
                good = (ui >= 0) & (rp <= v["ue"][np.clip(ui, 0, None)])
                if not good.any():
                    continue
                ui2, call2, rp2 = ui[good], call[good], rp[good]
                inu = (rp2 - v["us"][ui2] + 1) if st["strand"] == "+" else (v["ue"][ui2] - rp2 + 1)
                gate = 0.5 * MED.get(aid, 0)
                L = st["ul"]
                for u in np.unique(ui2):
                    w = ui2 == u
                    n_ = int(w.sum())
                    if n_ < gate:
                        st["cnt"]["덮임미달"] += 1
                        continue
                    beta = float(call2[w].mean())
                    cls = "LOW" if beta < LOWC else ("HIGH" if beta >= HIGHC else "MID")
                    pos = np.clip(inu[w].astype(np.int64), 1, L)
                    np.add.at(st["acc"]["ALL"][0], pos, 1.0)
                    np.add.at(st["acc"]["ALL"][1], pos, call2[w])
                    st["cnt"]["ALL"] += 1
                    if cls == "MID":
                        st["cnt"]["MID"] += 1
                    else:
                        np.add.at(st["acc"][cls][0], pos, 1.0)
                        np.add.at(st["acc"][cls][1], pos, call2[w])
                        st["cnt"][cls] += 1
    for aid in keep:
        st = ST[aid]
        v = Asp[aid]
        L = st["ul"]
        bw = max(1, int(round(L / 110.0)))
        edges = np.arange(1, L + 1 + bw, bw)
        idx = np.digitize(np.arange(1, L + 1), edges) - 1
        ctr = np.array([np.mean(np.arange(1, L + 1)[idx == j]) if (idx == j).any() else np.nan
                        for j in range(len(edges) - 1)])
        for cls in ("ALL", "HIGH", "LOW"):
            n = st["acc"][cls][0][1:L + 1]
            m = st["acc"][cls][1][1:L + 1]
            nb = np.bincount(idx, weights=n, minlength=len(edges) - 1)[:len(edges) - 1]
            mb = np.bincount(idx, weights=m, minlength=len(edges) - 1)[:len(edges) - 1]
            for c, nn, mm in zip(ctr, nb, mb):
                if nn > 0:
                    WU.append((ko, aid, LAB[aid], cls, L, bw, round(float(c), 1), int(nn), int(mm), round(mm / nn, 6)))
        for (side, inside), (n, mth) in LW[aid].items():
            e = win_edges(inside, st["imax"])
            for j in range(len(e)):
                if n[j] > 0:
                    LS.append((ko, aid, LAB[aid], side, "array 안" if inside else "flank", float(e[j]),
                               int(n[j]), int(mth[j]), round(mth[j] / n[j], 6)))
        # 유전자 구간(유닛 좌표) · 부분 유닛 · 유닛 수
        gw = []
        for i in range(len(v["us"])):
            for (x, y) in v["genes"][i]:
                gw.append((x - v["us"][i] + 1, y - v["us"][i] + 1) if st["strand"] == "+"
                          else (v["ue"][i] - y + 1, v["ue"][i] - x + 1))
        MET.append(dict(종=ko, array_id=aid, 이름표=LAB[aid], chrom=v["chrom"], strand=st["strand"],
                        array시작=st["a0"], array끝=st["a1"], 유닛길이=L, 유닛수=int(v["pos"].max()),
                        유전자시작=(float(np.median([x[0] for x in gw])) if gw else np.nan),
                        유전자끝=(float(np.median([x[1] for x in gw])) if gw else np.nan),
                        부분유닛_낮은쪽=float(v["us"].min() - st["a0"]), 부분유닛_높은쪽=float(st["a1"] - v["ue"].max()),
                        유닛1쪽=("낮은 좌표" if (v["us"][int(np.argmin(v["pos"]))] - st["a0"]) < (st["a1"] - v["ue"][int(np.argmin(v["pos"]))]) else "높은 좌표"),
                        read유닛_ALL=st["cnt"]["ALL"], read유닛_HIGH=st["cnt"]["HIGH"],
                        read유닛_LOW=st["cnt"]["LOW"], read유닛_MID=st["cnt"]["MID"]))
    P("%-10s array %2d · read-유닛 %s" % (ko, len(keep), format(sum(ST[a]["cnt"]["ALL"] for a in keep), ",")))

WUD = pd.DataFrame(WU, columns=["종", "array_id", "이름표", "등급", "유닛길이", "칸폭", "위치", "호출수", "메틸호출수", "beta"])
LSD = pd.DataFrame(LS, columns=["종", "array_id", "이름표", "쪽", "위치", "d", "호출수", "메틸호출수", "beta"])
MTD = pd.DataFrame(MET)
TMP = os.environ.get("TEMP", ".") + "/x34"
os.makedirs(TMP, exist_ok=True)
WUD.to_csv(TMP + "/PERARRAY_WITHIN_UNIT.tsv.gz", sep="\t", index=False, encoding="utf8", lineterminator="\n")
LSD.to_csv(TMP + "/PERARRAY_LANDSCAPE.tsv.gz", sep="\t", index=False, encoding="utf8", lineterminator="\n")
MTD.to_csv(TMP + "/PERARRAY_META.tsv", sep="\t", index=False, encoding="utf8", lineterminator="\n")
for f, rel in ((TMP + "/PERARRAY_WITHIN_UNIT.tsv.gz", "tables/PERARRAY_WITHIN_UNIT.tsv.gz"),
               (TMP + "/PERARRAY_LANDSCAPE.tsv.gz", "tables/PERARRAY_LANDSCAPE.tsv.gz"),
               (TMP + "/PERARRAY_META.tsv", "tables/PERARRAY_META.tsv")):
    arch.put(f, rel, REASON)
P("")
P("유닛안 %s 행 · 지형 %s 행 · array %d" % (format(len(WUD), ","), format(len(LSD), ","), len(MTD)))
LOG.close()
