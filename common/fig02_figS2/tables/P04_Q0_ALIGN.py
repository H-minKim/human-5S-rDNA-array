# -*- coding: utf-8 -*-
r"""Protocol 04 section 2 Q0, alignment (hard gate for Q4)
   Aligner: edlib installation failed (Python 3.13.14, no wheel, no C compiler)
             => parasail 1.3.4 is used. Citable:
               Daily J. Parasail: SIMD C library for global, semi-global, and local pairwise
               sequence alignments. BMC Bioinformatics 2016;17:81. DOI 10.1186/s12859-016-0930-z
   Method: section 2.3 consensus profile projection, every unit is globally aligned (NW) to CAN (2,240) individually
   Scores: match +2, mismatch -3, gap open 10, gap extend 1 (recorded in the block)
   Left normalization: indel placement inside repeat blocks is ambiguous => shift to the leftmost position
   QC set (section 2.4): (a) back-projection reconstruction (b) token reproduction (c) orientation
   Output: only the list of differences per unit is saved (whole sequences are not saved, to protect array sequences)
"""
import os, sys, io, re, time, collections, tarfile, pickle
import numpy as np
import parasail
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, __import__("os").environ.get("V11_CODE", __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), __import__("os").pardir)))
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "OUT_P04")
os.makedirs(OUT, exist_ok=True)
assert not OUT.startswith(("F:", "f:"))
TAR = (__import__("os").environ["RDNA_DATA_ROOT"] + r"/for_paper/08_versions/v11/12_flank/final/flank_final2.tar")
SEED = 20261010
MATCH, MISM, GO, GE = 2, -3, 10, 1
import v11_io as G

CAN = G.CAN
LU = len(CAN)
MAT = parasail.matrix_create("ACGT", MATCH, MISM)
COMP = str.maketrans("ACGTN", "TGCAN")
rc = lambda s: s.translate(COMP)[::-1]
CIG = re.compile(r"(\d+)([=XIDMS])")


def tsv(p):
    L = [l.rstrip("\n").split("\t") for l in io.open(p, encoding="utf8") if l.strip() and not l.startswith("#")]
    return {c: i for i, c in enumerate(L[0])}, L[1:]


# ---------- ruler ----------
U, _, _ = G.load("h1")
VPU = sum(len(r[4]) for v in U.values() for r in v) / sum(len(v) for v in U.values())
print("RULER var_per_unit = %.6f  (정본 6.866)" % VPU)
if abs(VPU - 6.868174) >= 0.000001: raise AssertionError("자가 다르다 (541 값 6.868174) — 훅 확인")   # 541: the value was 6.866 at 557
H1 = {("%s" % k[0], "%s" % k[1]) for k in U}
TOKS = {}
for k, rs in U.items():
    for r in rs:
        TOKS[("%s" % k[0], "%s" % k[1], r[6]["array_pos"])] = set(r[4])
print("   H1 슬롯 %d · 개체 %d · 유닛 %d" % (len(H1), len({k[0] for k in H1}), len(TOKS)))

MU = collections.defaultdict(list)
for f, Hx in G._rows(G.MAPF):
    key = (f[0], f[Hx["slot"]])
    if key in H1 and f[Hx["seg_type"]] == "UNIT":
        MU[key].append((int(f[Hx["array_pos"]]), int(f[Hx["unit_start"]]), int(f[Hx["unit_end"]]),
                        f[Hx["varlist"]], f[Hx["verify"]]))
for k in MU:
    MU[k].sort()
Hs, Rs = tsv((__import__("os").environ["RDNA_DATA_ROOT"] + r"/for_paper/08_versions/v11-fix1/03_maps/g8canonv11f1/CONSENSUS_SUM_G8_full.tsv"))
RCF = {(r[0], r[1]): r[Hs["revcomp"]] for r in Rs}
Hb, Rb = tsv((__import__("os").environ["RDNA_DATA_ROOT"] + r"/v11_readvar/02_tables/FASTA_REF_MAP_2026-10-07b.tsv"))
fa2slot = {}
for r in Rb:
    if "길이" in r[Hb["method"]]:
        rn = r[Hb["ref_name"]]
        if "_h" in rn:
            s, h = rn.rsplit("_h", 1)
            fa2slot[r[Hb["fasta"]]] = (s, h)
# analysis slots that the 2026-10-07 FASTA map does not list (their reference sequence was replaced later in v11 curation).
#   Their file name comes from the reference node table. A file that is not in the tar is read from the curated reference folder
#   (its whole sequence is the array, so it enters with flank5=0 and flank3=0) and is taken after the tar members, in slot order.
_NODE = (__import__("os").environ["RDNA_DATA_ROOT"] + r"/for_paper/08_versions/v11-add2_trio/02_canonical/HAP_SEQ_NODE_v11+add1+add2.tsv")
_node = {(r[0], r[1]): (r[2], r[4]) for r in (l.rstrip("\n").split("\t") for l in io.open(_NODE, encoding="utf8")) if len(r) > 4}
_known = set(fa2slot.values())
_EXRAW = []
for _k in sorted(MU):
    if _k in _known: continue
    _fa, _d = _node[_k]; _stem = _fa[:-3]; fa2slot[_stem] = _k
    with io.open((__import__("os").environ["RDNA_DATA_ROOT"] + "/for_paper/08_versions/" + _d + "/" + _fa), "rb") as _f:
        _b = _f.read()
    _seq = _b[_b.index(b"\n") + 1:].replace(b"\n", b"").replace(b"\r", b"").upper()
    _EXRAW.append((_stem, b">" + _stem.encode("ascii") + b"|flank5=0|flank3=0\n" + _seq + b"\n"))
print("   slots added from the reference node table %d : %s" % (len(_EXRAW), " ".join("%s|%s" % fa2slot[x[0]] for x in _EXRAW)))
HDR = re.compile(r"flank5=(\d+)\|flank3=(\d+)")


def align_map(u):
    """Globally aligns unit u to CAN -> (sub, dele, ins, score)
       sub  : {can_pos(1-based) : unit base}            only positions that differ from CAN
       dele : {can_pos : True}                          CAN positions missing from the unit
       ins  : [(can_pos_after, inserted sequence)]      between CAN positions p and p+1
       left normalization applied"""
    r = parasail.nw_trace_striped_16(u, CAN, GO, GE, MAT)
    cg = r.cigar.decode
    if isinstance(cg, bytes):
        cg = cg.decode("ascii")
    qi = ri = 0
    sub, dele, ins = {}, {}, []
    for n, op in CIG.findall(cg):
        n = int(n)
        if op in "=X":
            for t in range(n):
                if u[qi + t] != CAN[ri + t]:
                    sub[ri + t + 1] = u[qi + t]
            qi += n
            ri += n
        elif op == "M":
            for t in range(n):
                if u[qi + t] != CAN[ri + t]:
                    sub[ri + t + 1] = u[qi + t]
            qi += n
            ri += n
        elif op == "D":                     # only in CAN = deletion in the unit
            for t in range(n):
                dele[ri + t + 1] = True
            ri += n
        elif op == "I":                     # only in the unit = insertion
            ins.append((ri, u[qi:qi + n]))
            qi += n
        else:
            raise ValueError("미지의 CIGAR op %s" % op)
    # left normalization, deletions
    if dele:
        runs = []
        ps = sorted(dele)
        a = ps[0]
        for p in ps[1:] + [None]:
            if p is None or p != a + 1:
                runs.append((ps[0] if not runs else None, a))
                a = p
                ps = [p] if p else []
            else:
                a = p
        # simple implementation: regroup consecutive runs and shift them left
        segs = []
        cur = []
        for p in sorted(dele):
            if cur and p == cur[-1] + 1:
                cur.append(p)
            else:
                if cur:
                    segs.append(cur)
                cur = [p]
        if cur:
            segs.append(cur)
        dele = {}
        for sg in segs:
            s0, L = sg[0], len(sg)
            while s0 > 1 and CAN[s0 - 2] == CAN[s0 - 2 + L]:
                s0 -= 1
            for t in range(L):
                dele[s0 + t] = True
    # left normalization, insertions
    ins2 = []
    for p, s in ins:
        while p > 0 and s and CAN[p - 1] == s[-1]:
            s = CAN[p - 1] + s[:-1]
            p -= 1
        ins2.append((p, s))
    return sub, dele, ins2, int(r.score)


def rebuild(sub, dele, ins):
    """Back-projection: reconstructs the unit sequence from the mapping"""
    ia = collections.defaultdict(str)
    for p, s in ins:
        ia[p] += s
    out = [ia.get(0, "")]
    for p in range(1, LU + 1):
        if p not in dele:
            out.append(sub.get(p, CAN[p - 1]))
        if p in ia:
            out.append(ia[p])
    return "".join(out)


def to_tokens(sub, dele, ins):
    """Builds v4-style tokens from the alignment (to compare with the H1 tokens)"""
    t = set()
    for p, b in sub.items():
        t.add("%d%s%s" % (p, CAN[p - 1], b))
    segs = []
    cur = []
    for p in sorted(dele):
        if cur and p == cur[-1] + 1:
            cur.append(p)
        else:
            if cur:
                segs.append(cur)
            cur = [p]
    if cur:
        segs.append(cur)
    for sg in segs:
        t.add("D%d:%d" % (sg[0], len(sg)))
    for p, s in ins:
        t.add("I%d:%d" % (p + 1, len(s)))
    return t


print("\n■ tar 한 패스 — 정렬 (parasail nw_trace · match %d mism %d open %d ext %d)" % (MATCH, MISM, GO, GE))
DIFF = {}
ori_ok = ori_n = 0
qc_back = qc_back_n = 0
tok_exact = tok_n = 0
tok_sub_only = 0
sub_mismatch_detail = collections.Counter()
nlen_bad = 0
t0 = time.time()
tf = tarfile.open(TAR, "r:")
done = 0
def _items():
    seen = set()
    for m in tf.getmembers():
        if not m.isfile() or not m.name.endswith(".fa"):
            continue
        stem = os.path.basename(m.name)[:-3]; seen.add(stem)
        yield stem, (lambda m=m: tf.extractfile(m).read())
    for stem, raw in _EXRAW:          # reference files that are not in the tar
        if stem not in seen: yield stem, (lambda raw=raw: raw)
for stem, _get in _items():
    slot = fa2slot.get(stem)
    if slot is None or slot not in MU:
        continue
    raw = _get()
    nl = raw.index(b"\n")
    mm = HDR.search(raw[:nl].decode("ascii", "replace"))
    if not mm:
        continue
    f5, f3 = int(mm.group(1)), int(mm.group(2))
    seq = raw[nl + 1:].replace(b"\n", b"").replace(b"\r", b"").decode("ascii", "replace").upper()
    arr0 = seq[f5:len(seq) - f3]
    arr = rc(arr0) if RCF.get(slot, "N") == "Y" else arr0
    for ap, us, ue, vl, vf in MU[slot]:
        if vf != "OK":
            continue
        u = arr[us - 1:ue]
        if len(u) != ue - us + 1 or "N" in u:
            nlen_bad += 1
            continue
        sub, dele, ins, sc = align_map(u)
        # (c) orientation, on a sample
        if ori_n < 2000:
            ori_n += 1
            r2 = parasail.nw_striped_16(rc(u), CAN, GO, GE, MAT)
            if sc >= r2.score:
                ori_ok += 1
        # (a) back-projection
        if qc_back_n < 51801:
            qc_back_n += 1
            if rebuild(sub, dele, ins) == u:
                qc_back += 1
        # (b) token reproduction
        key = (slot[0], slot[1], ap)
        if key in TOKS:
            tok_n += 1
            got = to_tokens(sub, dele, ins)
            want = TOKS[key]
            gs = {x for x in got if re.match(r"^\d+[ACGT][ACGT]$", x)}
            ws = {x for x in want if re.match(r"^R?\d+[ACGT][ACGT]$", x)}
            ws = {x[1:] if x[0] == "R" else x for x in ws}
            if got == want:
                tok_exact += 1
            if gs == ws:
                tok_sub_only += 1
            else:
                sub_mismatch_detail["got_only"] += len(gs - ws)
                sub_mismatch_detail["want_only"] += len(ws - gs)
        DIFF[(slot[0], slot[1], ap)] = (sub, dele, ins, len(u))
        done += 1
    if done and done % 10000 < len(MU[slot]):
        print("   ... 유닛 %d · %.0fs" % (done, time.time() - t0), flush=True)
tf.close()
print("\n■ 정렬 끝 : 유닛 %d · %.0fs · 건너뜀(길이/N) %d" % (len(DIFF), time.time() - t0, nlen_bad))
print("\n■ QC 셋 (§2.4)")
print("   ⒜ 역투영 복원 : %d / %d (%.6f)" % (qc_back, qc_back_n, qc_back / max(qc_back_n, 1)))
print("   ⒝ 토큰 재현 : 전체 일치 %d / %d (%.5f) · 치환만 일치 %d / %d (%.5f)"
      % (tok_exact, tok_n, tok_exact / max(tok_n, 1), tok_sub_only, tok_n, tok_sub_only / max(tok_n, 1)))
print("      불일치 세부 : 정렬에만 %d · 토큰에만 %d" % (sub_mismatch_detail["got_only"], sub_mismatch_detail["want_only"]))
print("   ⒞ 방향 : 정방향 우세 %d / %d (%.5f)" % (ori_ok, ori_n, ori_ok / max(ori_n, 1)))
nsub = np.array([len(v[0]) for v in DIFF.values()])
ndel = np.array([len(v[1]) for v in DIFF.values()])
nins = np.array([sum(len(s) for _, s in v[2]) for v in DIFF.values()])
ulen = np.array([v[3] for v in DIFF.values()])
print("\n■ 정렬 요약 : 치환 median %d (평균 %.3f) · 결실 bp median %d (평균 %.3f) · 삽입 bp median %d (평균 %.3f)"
      % (np.median(nsub), nsub.mean(), np.median(ndel), ndel.mean(), np.median(nins), nins.mean()))
print("   unit 길이 : median %d · 10/90 %d/%d · min %d max %d · 2,240 인 몫 %.4f"
      % (np.median(ulen), *np.percentile(ulen, [10, 90]).astype(int), ulen.min(), ulen.max(), (ulen == LU).mean()))
print("   검산 : 길이 = 2240 − 결실 + 삽입 ? 일치 %d / %d"
      % (int((ulen == LU - ndel + nins).sum()), len(ulen)))
with open(os.path.join(OUT, "Q0_alignment.pkl"), "wb") as f:
    pickle.dump({"diff": DIFF, "params": dict(match=MATCH, mism=MISM, go=GO, ge=GE, seed=SEED,
                 aligner="parasail 1.3.4 nw_trace_striped_16", vpu=VPU)}, f, protocol=4)
g = io.open(os.path.join(OUT, "Q0_align_qc.tsv"), "w", encoding="utf8", newline="\n")
g.write("# P04 Q0 정렬 QC · parasail 1.3.4 nw_trace_striped_16 · match %d mism %d open %d ext %d · seed %d · 2026-10-10\n"
        % (MATCH, MISM, GO, GE, SEED))
g.write("# edlib 설치 실패 (Python 3.13.14 wheel 없음 · C 컴파일러 없음) ⇒ parasail (Daily 2016 BMC Bioinformatics 17:81)\n")
g.write("# 좌파 정규화 적용 (결실·삽입 둘 다)\n")
g.write("item\tvalue\n")
for k_, v_ in (("var_per_unit", VPU), ("n_unit_aligned", len(DIFF)), ("n_unit_skipped", nlen_bad),
               ("qc_backproject_ok", qc_back), ("qc_backproject_n", qc_back_n),
               ("qc_token_exact", tok_exact), ("qc_token_sub_only", tok_sub_only), ("qc_token_n", tok_n),
               ("qc_orient_fwd_wins", ori_ok), ("qc_orient_n", ori_n),
               ("mean_sub", nsub.mean()), ("mean_del_bp", ndel.mean()), ("mean_ins_bp", nins.mean()),
               ("frac_len_2240", (ulen == LU).mean()),
               ("len_identity_check", int((ulen == LU - ndel + nins).sum()))):
    g.write("%s\t%s\n" % (k_, v_))
g.close()
print("\n산출 :", os.path.join(OUT, "Q0_alignment.pkl"), "·", os.path.join(OUT, "Q0_align_qc.tsv"))
