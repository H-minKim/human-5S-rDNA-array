# -*- coding: utf-8 -*-
# Collecting per array, with these rules:
#   1) the classes are that array's top and bottom quartile, not fixed thresholds at 0.35 and 0.70, with the quartiles taken from that array's read-unit beta;
#   2) the maternal and paternal copies are merged, for the zebra finch haplotype pairs, as one locus;
#   3) every species is included, the siamang among them, although the siamang is noted in the table as having been basecalled with dorado instead;
#   The remaining calling rules: 5mC at ML >= 128, flag & 0x904 excluded, no mapq filter, n_calls >= 0.5 times that locus's median, and the bonobo on guppy 6.3.8.
#   Outputs: tables/PERARRAY_WITHIN_UNIT_v3.tsv.gz, tables/PERARRAY_META_v3.tsv and tables/PERARRAY_UNIT_BETA_v3.tsv.gz for the unit beta distribution.
import io, os, sys, glob, re, struct, collections
import numpy as np
import pandas as pd

import os as _os_need, sys as _sys_need


def _need(k):
    """No defaults are set: if an environment variable is missing the script stops instead of silently reading the wrong location."""
    v = _os_need.environ.get(k)
    if not v:
        _sys_need.exit("environment variable %s is not set (the deposited copy has no default; "
                       "see README.md, section 'Environment variables')" % k)
    return v


NEW = _need("RDNA_OUT_DIR")
B = _need("RDNA_DATA_ROOT")
OLD = _need("RDNA_PREV_OUT_DIR")
sys.path.insert(0, NEW + "/code")
import arch, tab_en  # noqa
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
LOG = io.open(NEW + "/logs/X42_PERARRAY_V2.txt", "w", encoding="utf8", newline="\n")
REASON = "array 마다 v3 — 순위 기준 상/하위 25%(동점 처리) · mat/pat 합침 · 시아망 포함"
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
ORI = pd.read_csv(tab_en.open_ko(NEW + "/tables/ORIENTATION.tsv"), sep="\t").set_index("array_id")
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
    STRAND = {}                                   # the array strand from the g8 map, used to bring the within-unit coordinate into the g8 frame
    for f_ in sorted(glob.glob(B + "/02_consensus_map/%s/*.g8map.tsv" % sp)):
        if not re.search(r"_[0-9]+_[0-9]+[.]g8map[.]tsv$", f_.replace("\\", "/")):
            continue
        g_ = pd.read_csv(f_, sep="\t", usecols=["array_id", "strand"])
        STRAND[g_.array_id.iloc[0]] = g_.strand.iloc[0]
    # ---- unit beta, per read x unit: the basis of the quartiles and the filter
    fs = sorted(glob.glob(B + "/03_reads/%s/*.read_unit.tsv" % sp) + glob.glob(B + "/03_reads/%s/*.read_unit.tsv.gz" % sp))
    D = pd.concat([pd.read_csv(f, sep="\t") for f in fs], ignore_index=True)
    D = D[D.context == "CG"].copy()
    med = D.groupby("array_id").n_calls.transform("median")      # the median is taken over all the reads
    D = D[D.n_calls >= 0.5 * med].copy()
    if sp == "bonobo_mPanPan1":
        D = D[~D.read_id.isin(B637)].copy()
    D["b"] = D.n_5mC / D.n_calls
    # ---- defining the locus: for the zebra finch the maternal and paternal copies are merged
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
        # The maternal and paternal copies are merged by chromosome name, with the _mat and _pat suffix stripped, pairing them in order of start coordinate
        #   (the unit count differs between haplotypes and so cannot serve as the pairing key: chr4 has 132 against 77 and chr9 70 against 41).
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
    # the labels are made unique, so that two different arrays cannot fall under one label
    #   (zebrafish 36,223,097 and 36,238,873 would both render as "chr4 · 36.2 Mb" and be merged into one locus).
    cnt_lab = collections.Counter(k for k, _a0, _a1, _l in grp_of.values())
    dup = {k for k, c in cnt_lab.items() if c > 1}
    if dup and sp != "zfinch_v03_dip":            # the zebra finch maternal and paternal copies are merged deliberately and are left alone here
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
        assert all(len(v) == 1 for v in groups.values()), "labels still overlap: %s" % [k for k, v in groups.items() if len(v) > 1]
    # ---- the quartiles, per locus
    cut, CLS, TIE = {}, {}, {}
    rngq = np.random.default_rng(0)
    for key, aids in groups.items():
        sub = D[D.array_id.isin(aids)]
        b = sub.b.values
        q1, q3 = float(np.quantile(b, 0.25)), float(np.quantile(b, 0.75))
        cut[key] = (q1, q3)
        TIE[key] = float((b >= 0.999).mean())
        # Ranking: ties are shuffled in a random order from seed 0, so that the bottom and top 25% are cut exactly.
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
    # ---- one pass over the BAM for the within-unit profile, accumulated per locus
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
                # the within-unit coordinate is in the g8 frame: the array strand must be applied so that the gene falls at the end of the unit
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
                    if cls is None:          # a read-unit absent from the read_unit table, which is rare, is counted only in ALL and given no class
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
        # Nothing is binned: the CpG positions that carry a value are emitted as they are, so the lines do not break at empty bins.
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
    P("%-10s groups %2d, read-units %s" % (ko, len(groups), format(sum(ACC[k]["cnt"]["ALL"] for k in groups), ",")))

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
P("groups %d, within-unit %s rows" % (len(MT), format(len(W), ",")))
LOG.close()
