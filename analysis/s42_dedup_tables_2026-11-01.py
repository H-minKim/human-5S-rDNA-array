# -*- coding: utf-8 -*-
#  S5 입력표 세 개를 NA21133 한 벌로 다시 낸다 (2026-11-01 사용자 ①)
#    규칙 = 꾸러미 disp_17_prob_all.py (PROB_A4_SITE) · disp_19_a4_persite_tone.py (둘) 의 식 그대로
#    ① 중복 그대로 돌려 꾸러미 입력표 셋을 바이트까지 재현 → <OUT>/dup/
#    ② 한 벌만 → <OUT>/dedup/  (원본 표는 건드리지 않는다)
#    ③ 좌위 원소 대조 : 개체×자리마다 (모든 CpG 수 · 메틸 수) 다중집합 ↔ 58 판 STRUCT_LOCI (CpG ≥ 10)
#  쓰기 : python -B s42_dedup_tables_2026-11-01.py <OUT>
import collections
import glob
import gzip
import hashlib
import io
import os
import statistics
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
def _need(_v):
    """The deposited copy has no default paths. Stop when a root variable is missing."""
    _x = os.environ.get(_v)
    if not _x:
        sys.exit("environment variable %s is not set; see README.md section Figures after 30 October 2026" % _v)
    return _x


OUT = os.path.abspath(sys.argv[1])
R = _need("DISPERSED_ROOT")
N = os.path.join(R, "02_output_tables_node_CMD03")
PK = os.path.join(R, "05_paper_bundle_2026-10-09", "input_tables")
R0 = os.path.join(_need("RDNA_VERSION_ROOT"), "05_checks", "2026-09-17_TUC_dropped_reads_5S_content_methylation")
BW, NB, MINCPG, MINN = 0.02, 50, 10, 150
DUP, KEEPNODE = "NA21133", _need("S42_KEEP_NODE")
JA = (("하드콜", "beta_all_hard"), ("확률평균", "beta_all_prob"))
GRPS = ("전체", "색조 낮음", "중간", "색조 높음")


def rd(p):
    L = io.open(p, encoding="utf8").read().rstrip("\n").split("\n")
    H = L[0].split("\t")
    return [dict(zip(H, x.split("\t"))) for x in L[1:]]


INFO = {r["site"]: (r["chr"], int(r["start"])) for r in rd(os.path.join(R, "02_output_tables_coding_region", "DISP_SITES_TABLE.tsv"))}
TONE = {r["sample"]: float(r["tone"]) for r in rd(os.path.join(N, "TONE_PERSAMPLE.tsv")) if r["tone"] != "NA"}
q1, q2 = statistics.quantiles(sorted(TONE.values()), n=3)
print("자리 %d · 개체 색조 %d · 경계 %.4f %.4f" % (len(INFO), len(TONE), q1, q2))


def grp_of(s):
    t = TONE.get(s)
    return None if t is None else ("색조 낮음" if t <= q1 else ("색조 높음" if t > q2 else "중간"))


def modes(v):
    h = [0] * NB
    for x in v:
        h[min(NB - 1, max(0, int(x / BW)))] += 1
    lo, hi = int(0.2 / BW), int(0.8 / BW)
    vi = min(range(lo, hi), key=lambda i: h[i])
    a, b = h[:vi], h[vi:]
    lm = (a.index(max(a)) + 0.5) * BW if a and max(a) else float("nan")
    hm = (vi + b.index(max(b)) + 0.5) * BW if b and max(b) else float("nan")
    pk = min(max(a) if a else 0, max(b) if b else 0)
    dep = (h[vi] / pk) if pk else float("nan")
    return (vi + 0.5) * BW, lm, hm, sum(b) / float(len(v)), statistics.mean(v), dep, h


SB = {"dup": collections.defaultdict(list), "dedup": collections.defaultdict(list)}
ST = {"dup": collections.defaultdict(list), "dedup": collections.defaultdict(list)}
EL = collections.defaultdict(collections.Counter)     # (sample, site) -> (n_cpg_all, 메틸 수) 다중집합 (한 벌)
for p in sorted(glob.glob(os.path.join(N, "*", "LOCUS_BETA_V2.tsv"))):
    node = os.path.basename(os.path.dirname(p))
    for r in rd(p):
        n = int(r["n_cpg_all"])
        if n < MINCPG or r["site"] not in INFO:
            continue
        g = grp_of(r["sample"])
        for basis in ("dup", "dedup"):
            if basis == "dedup" and r["sample"] == DUP and node != KEEPNODE:
                continue
            for ja, key in JA:
                b = float(r[key])
                SB[basis][(r["site"], ja)].append(b)
                ST[basis][(r["site"], ja, "전체")].append(b)
                if g:
                    ST[basis][(r["site"], ja, g)].append(b)
        if not (r["sample"] == DUP and node != KEEPNODE):
            EL[(r["sample"], r["site"])][(n, int(round(float(r["beta_all_hard"]) * n)))] += 1


def dump(p, head, rows):
    io.open(p, "w", encoding="utf8", newline="\n").write(head + "\n" + "\n".join("\t".join(str(x) for x in r) for r in rows) + "\n")


NAMES = ("PROB_A4_SITE.tsv", "PROB_A4_SITE_TONE.tsv", "PROB_A4_DIST_SITE_TONE.tsv")
for basis in ("dup", "dedup"):
    d = os.path.join(OUT, basis)
    os.makedirs(d, exist_ok=True)
    rs, rt, rd_ = [], [], []
    for s in sorted(INFO):
        ch, st = INFO[s]
        for ja, _k in JA:
            v = SB[basis].get((s, ja), [])
            if len(v) >= 200:
                val, lm, hm, fh, mn, dep, h = modes(v)
                ok = (lm < 0.35 and hm > 0.65 and 0.05 <= fh <= 0.95 and dep <= 0.5)
                rs.append([s, ch, "%.2f" % (st / 1e6), ja, len(v), "%.3f" % lm, "%.3f" % hm, "%.3f" % val,
                           "%.4f" % fh, "%.3f" % dep, "이봉" if ok else "단봉"])
    for s in sorted(INFO):
        ch, st = INFO[s]
        for ja, _k in JA:
            for grp in GRPS:
                v = ST[basis].get((s, ja, grp), [])
                if len(v) < MINN:
                    continue
                val, lm, hm, fh, mn, dep, h = modes(v)
                rt.append([s, ch, "%.2f" % (st / 1e6), ja, grp, "%.3f" % lm, "%.3f" % hm, "%.3f" % val,
                           "%.4f" % fh, "%.4f" % mn, "%.3f" % dep, len(v)])
                for i, c in enumerate(h):
                    rd_.append([s, ja, grp, "%.2f" % (i * BW), "%.5f" % (c / len(v) / BW), c])
    dump(os.path.join(d, NAMES[0]), "site\tchr\tMb\truler\tn_loci\tlow_mode\thigh_mode\tvalley\tfrac_high\tvalley_depth\tcall_A4", rs)
    dump(os.path.join(d, NAMES[1]), "site\tchr\tMb\truler\tgrp\tlow_mode\thigh_mode\tvalley\tfrac_high\tmean_beta\tvalley_depth\tn", rt)
    dump(os.path.join(d, NAMES[2]), "site\truler\tgrp\tbeta_lo\tdensity\tn", rd_)


def md5(p):
    return hashlib.md5(io.open(p, "rb").read()).hexdigest()


print("① 재현 (dup ↔ 꾸러미)")
okall = True
for f in NAMES:
    a, b, c = md5(os.path.join(OUT, "dup", f)), md5(os.path.join(PK, f)), md5(os.path.join(OUT, "dedup", f))
    okall &= a == b
    print("   %-30s dup %s · 꾸러미 %s · 같음 %s · dedup %s" % (f, a[:12], b[:12], a == b, c[:12]))
assert okall, "꾸러미 재현 실패"

#  dedup ↔ dup 줄마다 무엇이 바뀌었나
for f in NAMES[:2]:
    A = rd(os.path.join(OUT, "dup", f))
    B = rd(os.path.join(OUT, "dedup", f))
    assert len(A) == len(B)
    cnt = collections.Counter()
    for x, y in zip(A, B):
        for k in x:
            if x[k] != y[k]:
                cnt[k] += 1
    print("   %s 바뀐 칸 (열별) %s" % (f, dict(cnt)))
A = [r for r in rd(os.path.join(OUT, "dup", NAMES[0])) if r["ruler"] == "하드콜"]
B = [r for r in rd(os.path.join(OUT, "dedup", NAMES[0])) if r["ruler"] == "하드콜"]
print("   판정 바뀐 자리 :", [x["site"] for x, y in zip(A, B) if x["call_A4"] != y["call_A4"]],
      "· 빠지는 세 자리 n :", [(y["site"], x["n_loci"], y["n_loci"]) for x, y in zip(A, B) if y["site"] in ("S0285", "S0107", "S0108")])

#  ③ 원소 대조
SITE = {}
with gzip.open(os.path.join(R0, "CHM13_positions", "LOCI_SITE.tsv.gz"), "rt", encoding="utf8") as F:
    F.readline()
    for L in F:
        f = L.rstrip("\n").split("\t")
        if f[4] in INFO:
            SITE[(f[0], f[1], f[2])] = f[4]
E2 = collections.defaultdict(collections.Counter)
for p in sorted(glob.glob(os.path.join(R0, "structure_exact_cells", "per_individual", "*.STRUCT_LOCI.tsv.gz"))):
    smp = os.path.basename(p).split(".")[0]
    with gzip.open(p, "rt", encoding="utf8") as F:
        H = F.readline().rstrip("\n").split("\t")
        ic, il, ird, iub = H.index("class"), H.index("locus"), H.index("read_id"), H.index("ubins50_exact")
        for L in F:
            f = L.rstrip("\n").split("\t")
            if f[ic] == "ARRAY":
                continue
            s = SITE.get((smp, f[ird], f[il]))
            if s is None or f[iub] in ("", "NA"):
                continue
            n = k = 0
            for t in f[iub].split(";"):
                _b, a, c = t.split(":")
                n += int(a)
                k += int(c)
            if n >= MINCPG:
                E2[(smp, s)][(n, k)] += 1
S29 = set(x.split("\t")[0] for x in io.open(_need("S42SET29_SITES_TSV"), encoding="utf8").read().split("\n")[1:] if x)
for nm, ss in (("29 자리", S29), ("표 전체 자리", set(INFO))):
    keys = [k for k in set(EL) | set(E2) if k[1] in ss]
    dk = [k for k in keys if EL[k] != E2[k]]
    only_a = sum(sum((EL[k] - E2[k]).values()) for k in keys)
    only_b = sum(sum((E2[k] - EL[k]).values()) for k in keys)
    print("③ %s : 키 %d · 다른 키 %d · 좌위 %s ↔ %s · S5 쪽에만 %d · 58 판 쪽에만 %d" % (
        nm, len(keys), len(dk), format(sum(sum(EL[k].values()) for k in keys), ","),
        format(sum(sum(E2[k].values()) for k in keys), ","), only_a, only_b))
    for k in dk[:5]:
        print("    ", k, list((EL[k] - E2[k]).items())[:4], list((E2[k] - EL[k]).items())[:4])
