# -*- coding: utf-8 -*-
# array 마다 다시 모으기 v2 — 반영
#   ① 등급을 그 array 의 상위 25% / 하위 25% 로 (고정 문턱 0.35/0.70 대신 · 사분위는 그 array 의 read-유닛 β 에서)
#   ② mat / pat 를 합친다 (금화조 하플로 쌍 — 같은 좌위로 묶음)
#   ③ 모든 종 : 시아망 포함 (단 시아망은 호출기가 dorado 로 달라 표에 적어 둔다)
#   나머지 호출 규칙은 그대로 : 5mC = ML>=128 · flag&0x904 제외 · mapq 무여과 · n_calls >= 그 묶음 중앙값 x 0.5 · 보노보 6.3.8
#   산출 : tables/PERARRAY_WITHIN_UNIT_v3.tsv.gz · tables/PERARRAY_META_v3.tsv · tables/PERARRAY_UNIT_BETA_v3.tsv.gz(유닛 β 분포)
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
LOG = io.open(NEW + "/logs/X42_PERARRAY_V2.txt", "w", encoding="utf8", newline="\n")
REASON = "array 마다 v3 — 순위 기준 상/하위 25%(동점 처리) · mat/pat 합침 · 시아망 포함(2026-10-01)"
CUT = 128


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
BC = pd.read_csv(NEW + "/tables/BONOBO_READ_BASECALLER.tsv", sep="\t")
B637 = set(BC.loc[BC.basecaller.astype(str).str.contains("6.3.7"), "read_id"])
ZF = pd.read_csv(B + "/00_catalog/ZEBRAFISH_ALL_ARRAYS.tsv", sep="\t")
ZNAME = dict(zip(ZF.array_id, ZF.chrom_name))
FN = pd.read_csv(B + "/06_aux/zfinch_5S_array_meth_summary.tsv", sep="\t")
FNAME = {(t.chrom_genbank, int(t.array_start), int(t.array_end)): t.chrom_bigwig for t in FN.itertuples(index=False)}
ORI = pd.read_csv(NEW + "/tables/ORIENTATION.tsv", sep="\t").set_index("array_id")
SPK = [("고릴라", "gorilla_mGorGor1"), ("보노보", "bonobo_mPanPan1"), ("침팬지", "chimp_mPanTro3"),
       ("보르네오오랑우탄", "orang_bornean_mPonPyg2"), ("수마트라오랑우탄", "orang_sumatran_mPonAbe1"),
       ("시아망", "siamang_mSymSyn1"), ("애기장대", "arabidopsis_ColCC"),
       ("금화조", "zfinch_v03_dip"), ("제브라피시", "zebrafish_GRCz12tu")]


def label(aid, chrom, a0, a1):
    if aid in ZNAME and isinstance(ZNAME[aid], str):
        return "%s · %.1f Mb" % (ZNAME[aid], a0 / 1e6)
    k = (chrom, int(a0), int(a1))
    if k in FNAME:
        return str(FNAME[k])
    return "%s · %.1f Mb" % (chrom, a0 / 1e6)


WU, MET, UB = [], [], []
for ko, sp in SPK:
    Asp = load_species(sp)
    STRAND = {}                                   # g8map 의 array strand (유닛 안 좌표를 g8 틀로 맞출 때 씀)
    for f_ in sorted(glob.glob(B + "/02_consensus_map/%s/*.g8map.tsv" % sp)):
        if not re.search(r"_[0-9]+_[0-9]+[.]g8map[.]tsv$", f_.replace("\\", "/")):
            continue
        g_ = pd.read_csv(f_, sep="\t", usecols=["array_id", "strand"])
        STRAND[g_.array_id.iloc[0]] = g_.strand.iloc[0]
    # ---- 유닛 β (read x 유닛) : 사분위·여과 기준
    fs = sorted(glob.glob(B + "/03_reads/%s/*.read_unit.tsv" % sp) + glob.glob(B + "/03_reads/%s/*.read_unit.tsv.gz" % sp))
    D = pd.concat([pd.read_csv(f, sep="\t") for f in fs], ignore_index=True)
    D = D[D.context == "CG"].copy()
    med = D.groupby("array_id").n_calls.transform("median")      # 중앙값은 모든 read 에서
    D = D[D.n_calls >= 0.5 * med].copy()
    if sp == "bonobo_mPanPan1":
        D = D[~D.read_id.isin(B637)].copy()
    D["b"] = D.n_5mC / D.n_calls
    # ---- 묶음(좌위) 정의 : 금화조는 mat/pat 를 합친다
    grp_of = {}
    for aid in Asp:
        if aid not in set(D.array_id):
            continue
        v = Asp[aid]
        hd, _ = read_fa(B + "/01_canonical/%s/%s.fa" % (sp, aid))
        m = re.search(r"array=(\d+)-(\d+)", hd)
        a0, a1 = int(m.group(1)), int(m.group(2))
        lab = label(aid, v["chrom"], a0, a1)
        key = lab
        grp_of[aid] = (key, a0, a1, lab)
    if sp == "zfinch_v03_dip":
        # mat / pat 를 합친다 : 염색체 이름(_mat/_pat 뗀 것)마다, 시작 좌표 순서 로 짝을 맞춘다
        #   (유닛 수는 하플로마다 달라서 짝짓기 기준이 될 수 없다 — chr4 132/77 · chr9 70/41)
        byhap = collections.defaultdict(list)
        for aid, (key, a0, a1, lab) in grp_of.items():
            base = re.sub(r"_(mat|pat)$", "", lab)
            hap = "mat" if lab.endswith("_mat") else ("pat" if lab.endswith("_pat") else "?")
            byhap[(base, hap)].append((a0, aid))
        rank = {}
        for (base, hap), lst in byhap.items():
            for i, (a0, aid) in enumerate(sorted(lst)):
                rank[aid] = (base, i)
        for aid, (key, a0, a1, lab) in list(grp_of.items()):
            base, i = rank[aid]
            nsame = len([1 for (b, h) in byhap if b == base])
            newkey = "%s%s (mat+pat)" % (base, ("" if len([1 for x in byhap.get((base, "mat"), [])]) <= 1 else " #%d" % (i + 1)))
            grp_of[aid] = (newkey, a0, a1, lab)
    # 이름표 유일화 : 서로 다른 array 가 같은 이름표로 묶이면 안 된다
    #   (제브라 36,223,097 과 36,238,873 이 'chr4 · 36.2 Mb' 로 같아져 두 array 가 한 묶음이 되는 버그가 있었다)
    cnt_lab = collections.Counter(k for k, _a0, _a1, _l in grp_of.values())
    dup = {k for k, c in cnt_lab.items() if c > 1}
    if dup and sp != "zfinch_v03_dip":            # 금화조는 mat/pat 를 일부러 합치므로 건드리지 않는다
        for aid, (key, a0, a1, lab) in list(grp_of.items()):
            if key in dup:
                for nd in (2, 3, 4):
                    cand = "%s · %.*f Mb" % (key.split(" · ")[0], nd, a0 / 1e6)
                    if sum(1 for a2, (k2, b0, _b1, _l2) in grp_of.items()
                           if ("%s · %.*f Mb" % (k2.split(" · ")[0], nd, b0 / 1e6)) == cand) == 1:
                        grp_of[aid] = (cand, a0, a1, lab)
                        break
    groups = collections.defaultdict(list)
    for aid, (key, a0, a1, lab) in grp_of.items():
        groups[key].append(aid)
    if sp != "zfinch_v03_dip":
        assert all(len(v) == 1 for v in groups.values()), "이름표가 아직 겹친다 : %s" % [k for k, v in groups.items() if len(v) > 1]
    # ---- 묶음마다 사분위
    cut, CLS, TIE = {}, {}, {}
    rngq = np.random.default_rng(0)
    for key, aids in groups.items():
        sub = D[D.array_id.isin(aids)]
        b = sub.b.values
        q1, q3 = float(np.quantile(b, 0.25)), float(np.quantile(b, 0.75))
        cut[key] = (q1, q3)
        TIE[key] = float((b >= 0.999).mean())
        # 순위 기준 : 동점은 씨앗 0 의 무작위 순서로 섞어 아래 25% / 위 25% 를 정확히 가른다
        order = np.lexsort((rngq.random(len(b)), b))
        n = len(b)
        lo_n, hi_n = int(round(0.25 * n)), int(round(0.25 * n))
        lab = np.full(n, "MID", dtype=object)
        lab[order[:lo_n]] = "LOW"
        lab[order[n - hi_n:]] = "HIGH"
        for (aid_, rid_, ap_), lb in zip(zip(sub.array_id.values, sub.read_id.values, sub.array_pos.values), lab):
            CLS[(aid_, rid_, int(ap_))] = lb
        for a in aids:
            UB.append((ko, key, a, len(b), round(q1, 5), round(q3, 5)))
    # ---- BAM 한 번 훑어 유닛 안 프로파일 (묶음 단위로 누적)
    bychrom = collections.defaultdict(list)
    for aid in grp_of:
        bychrom[Asp[aid]["chrom"]].append(aid)
    MEDC = D.groupby("array_id").n_calls.median().to_dict()
    ACC = {}
    for key, aids in groups.items():
        L = int(np.median([np.median(Asp[a]["ue"] - Asp[a]["us"] + 1) for a in aids]))
        ACC[key] = dict(L=L, acc={k: np.zeros((2, L + 2)) for k in ("ALL", "HIGH", "LOW")}, cnt=collections.Counter())
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
                key = grp_of[aid][0]
                st = ACC[key]
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
                ui = np.searchsorted(v["us"], rp, side="right") - 1
                good = (ui >= 0) & (rp <= v["ue"][np.clip(ui, 0, None)])
                if not good.any():
                    continue
                ui2, call2, rp2 = ui[good], call[good], rp[good]
                # 유닛 안 좌표는 g8 틀 : array strand 를 맞춰 뒤집어야 유전자가 유닛 끝에 온다
                inu = (rp2 - v["us"][ui2] + 1) if STRAND[aid] == "+" else (v["ue"][ui2] - rp2 + 1)
                gate = 0.5 * MEDC.get(aid, 0)
                L = st["L"]
                q1, q3 = cut[key]
                for u in np.unique(ui2):
                    w = ui2 == u
                    n_ = int(w.sum())
                    if n_ < gate:
                        continue
                    cls = CLS.get((aid, rid, int(v["pos"][u])), None)
                    if cls is None:          # read_unit 표에 없는 read-유닛(드묾) — 등급 없이 ALL 에만
                        cls = "MID"
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
    for key, aids in groups.items():
        st = ACC[key]
        L = st["L"]
        # 칸으로 묶지 않는다 — 값이 있는 CpG 자리 그대로 내보낸다(빈 칸 때문에 선이 끊기던 문제)
        nall = st["acc"]["ALL"][0][1:L + 1]
        sites = np.flatnonzero(nall > 0) + 1
        for cls in ("ALL", "HIGH", "LOW"):
            n = st["acc"][cls][0][1:L + 1]; m = st["acc"][cls][1][1:L + 1]
            for p_ in sites:
                nn, mm = n[p_ - 1], m[p_ - 1]
                if nn > 0:
                    WU.append((ko, key, cls, L, 1, int(p_), int(nn), int(mm), round(mm / nn, 6)))
        v0 = Asp[aids[0]]
        gw = []
        for i in range(len(v0["us"])):
            for (x, y) in v0["genes"][i]:
                gw.append((x - v0["us"][i] + 1, y - v0["us"][i] + 1))
        q1, q3 = cut[key]
        MET.append(dict(종=ko, 묶음=key, array수=len(aids), array목록="|".join(sorted(aids)),
                        유닛길이=L, 유닛수합=int(sum(Asp[a]["pos"].max() for a in aids)),
                        read유닛_ALL=st["cnt"]["ALL"], read유닛_HIGH=st["cnt"]["HIGH"], read유닛_LOW=st["cnt"]["LOW"],
                        read유닛_MID=st["cnt"]["MID"], 하위25컷=round(q1, 4), 상위25컷=round(q3, 4), 동점_β1비율=round(TIE[key], 4),
                        유전자시작=(float(np.median([x[0] for x in gw])) if gw else np.nan),
                        유전자끝=(float(np.median([x[1] for x in gw])) if gw else np.nan),
                        호출기=("dorado 0.9.1(다름)" if sp == "siamang_mSymSyn1" else "")))
    P("%-10s 묶음 %2d · read-유닛 %s" % (ko, len(groups), format(sum(ACC[k]["cnt"]["ALL"] for k in groups), ",")))

W = pd.DataFrame(WU, columns=["종", "묶음", "등급", "유닛길이", "칸폭", "위치", "호출수", "메틸호출수", "beta"])
MT = pd.DataFrame(MET)
UBD = pd.DataFrame(UB, columns=["종", "묶음", "array_id", "read유닛", "하위25컷", "상위25컷"])
TMP = os.environ.get("TEMP", ".") + "/x42"
os.makedirs(TMP, exist_ok=True)
W.to_csv(TMP + "/PERARRAY_WITHIN_UNIT_v3.tsv.gz", sep="\t", index=False, encoding="utf8", lineterminator="\n")
MT.to_csv(TMP + "/PERARRAY_META_v3.tsv", sep="\t", index=False, encoding="utf8", lineterminator="\n")
UBD.to_csv(TMP + "/PERARRAY_UNIT_BETA_v3.tsv.gz", sep="\t", index=False, encoding="utf8", lineterminator="\n")
for f, rel in ((TMP + "/PERARRAY_WITHIN_UNIT_v3.tsv.gz", "tables/PERARRAY_WITHIN_UNIT_v3.tsv.gz"),
               (TMP + "/PERARRAY_META_v3.tsv", "tables/PERARRAY_META_v3.tsv"),
               (TMP + "/PERARRAY_UNIT_BETA_v3.tsv.gz", "tables/PERARRAY_UNIT_BETA_v3.tsv.gz")):
    arch.put(f, rel, REASON)
P("")
P("묶음 %d · 유닛안 %s 행" % (len(MT), format(len(W), ",")))
LOG.close()
