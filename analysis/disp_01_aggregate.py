# -*- coding: utf-8 -*-
# 흩어진 5S 좌위 → 자리(site) 단위 집계 (로컬 자료만, 읽기 전용)
#
# 입력 : $WORK_ROOT\05_checks\2026-09-17_TUC_dropped_reads_5S_content_methylation\
#   structure_exact_cells\per_individual\<S>.STRUCT_LOCI.tsv.gz  (좌위 1 행 · ubins50_exact = 50bp 칸 CpG:메틸)
#   CHM13_positions\LOCI_SITE.tsv.gz                (좌위 → 자리, core·mapq>=20)
# 출력 : dispersed_vs_array\DISP_SITE_*.tsv
#   · 자리 × 상태(low/mid/high/NA) × 칸 (n, meth) — 짝/홀 개체 반쪽 따로(반쪽 재현용)
#   · 자리 × 칸 덮음(좌위 수)
#   · 자리 × 개체 β (자리 안 이봉·개체 일관성용)
#   · ARRAY 대조(같은 코드) 상태 × 칸
# 인자 : 표본목록(쉼표) | ALL   태그
import io, os, sys, gzip, glob, collections

def _need(_v):
    """No defaults are set: the script stops if a variable is missing.
    이 글이 읽는 자료의 뿌리는 README.md 의
    "Running this code" 절에 적혀 있다."""
    import os as _o, sys as _s
    _x = _o.environ.get(_v)
    if not _x:
        _s.exit(
            "environment variable %s is not set; it is required by %s"
            % (_v, _o.path.basename(__file__))
            + "\n(a module stops the script that imports it)."
            + "\nThe deposited copy has no default; see README.md"
            + " section \"Running this code\"."
        )
    return _x


R = os.path.join(_need("RDNA_VERSION_ROOT"), "05_checks", "2026-09-17_TUC_dropped_reads_5S_content_methylation")
OUT = os.path.join(_need("RDNA_WORK_ROOT"), "dispersed_vs_array")
want = sys.argv[1] if len(sys.argv) > 1 else "ALL"
tag = sys.argv[2] if len(sys.argv) > 2 else "ALL"
wset = None if want == "ALL" else set(want.split(","))

SITE = {}
for i, L in enumerate(gzip.open(os.path.join(R, "CHM13_positions", "LOCI_SITE.tsv.gz"), "rt", encoding="utf8")):
    if i == 0:
        continue
    f = L.rstrip("\n").split("\t")
    if wset and f[0] not in wset:
        continue
    SITE[(f[0], f[1], f[2])] = f[4]
print("LOCI_SITE 좌위 %d" % len(SITE), flush=True)

files = sorted(glob.glob(os.path.join(R, "structure_exact_cells", "per_individual", "*.STRUCT_LOCI.tsv.gz")))
samples = [os.path.basename(p).split(".")[0] for p in files]
PAR = {s: i % 2 for i, s in enumerate(samples)}       # 반쪽 = 개체 순번 짝/홀

BINS = collections.defaultdict(lambda: [0, 0])     # (site, state, half, bin) -> [n, meth]
COV = collections.Counter()                         # (site, bin) -> 좌위 수
NLOC = collections.Counter()                        # site -> 좌위 수
GRP = collections.defaultdict(collections.Counter)  # site -> 구조 그룹
STC = collections.defaultdict(collections.Counter)  # site -> 상태 수
SB = collections.defaultdict(lambda: [0, 0])        # (site, sample) -> [cpg_spec, meth_spec]
LB = collections.defaultdict(list)                  # site -> 좌위 beta_spec (CpG>=10)
PROFCHK = collections.defaultdict(lambda: [0, 0])   # (group, bin) -> 파일럿 대조용 (자리 무관)
nsite = nnos = 0
for p in files:
    smp = os.path.basename(p).split(".")[0]
    if wset and smp not in wset:
        continue
    for i, L in enumerate(gzip.open(p, "rt", encoding="utf8")):
        if i == 0:
            continue
        f = L.rstrip("\n").split("\t")
        cls, grp, state, ub = f[2], f[19], f[18], f[20]
        bins = []
        if ub and ub != "NA":
            for t in ub.split(";"):
                b, n, k = t.split(":")
                bins.append((int(b), int(n), int(k)))
        for b, n, k in bins:
            v = PROFCHK[(grp, b)]
            v[0] += n
            v[1] += k
        if cls == "ARRAY":
            site = "ARRAY"
        else:
            site = SITE.get((smp, f[1], f[4]))
            if site is None:
                nnos += 1
                continue
        nsite += 1
        NLOC[site] += 1
        GRP[site][grp] += 1
        STC[site][state] += 1
        cs, ms = int(f[15]), int(f[16])
        v = SB[(site, smp)]
        v[0] += cs
        v[1] += ms
        if cs >= 10:
            LB[site].append(ms / float(cs))
        h = PAR[smp]
        for b, n, k in bins:
            w = BINS[(site, state, h, b)]
            w[0] += n
            w[1] += k
            COV[(site, b)] += 1
    print("  %s 끝 (자리 붙은 좌위 누계 %d · 자리 없음 %d)" % (smp, nsite, nnos), flush=True)

with io.open(os.path.join(OUT, "DISP_SITE_BINS_%s.tsv" % tag), "w", encoding="utf8", newline="\n") as O:
    O.write("site\tstate\thalf\tbin\tn\tmeth\n")
    for (s, st, h, b), (n, k) in sorted(BINS.items()):
        O.write("%s\t%s\t%d\t%d\t%d\t%d\n" % (s, st, h, b, n, k))
with io.open(os.path.join(OUT, "DISP_SITE_INFO_%s.tsv" % tag), "w", encoding="utf8", newline="\n") as O:
    O.write("site\tn_loci\tgroups\tstates\tcov_bins\n")
    for s in sorted(NLOC):
        O.write("%s\t%d\t%s\t%s\t%s\n" % (
            s, NLOC[s], ";".join("%s:%d" % x for x in GRP[s].most_common()),
            ";".join("%s:%d" % x for x in STC[s].most_common()),
            ";".join("%d:%d" % (b, COV[(s, b)]) for b in range(45) if COV[(s, b)])))
with io.open(os.path.join(OUT, "DISP_SITE_SAMPLE_%s.tsv" % tag), "w", encoding="utf8", newline="\n") as O:
    O.write("site\tsample\tcpg_spec\tmeth_spec\n")
    for (s, smp), (n, k) in sorted(SB.items()):
        O.write("%s\t%s\t%d\t%d\n" % (s, smp, n, k))
with io.open(os.path.join(OUT, "DISP_SITE_LOCUSBETA_%s.tsv" % tag), "w", encoding="utf8", newline="\n") as O:
    O.write("site\tbetas\n")
    for s in sorted(LB):
        O.write("%s\t%s\n" % (s, ",".join("%.3f" % x for x in LB[s])))
with io.open(os.path.join(OUT, "DISP_PROFCHK_%s.tsv" % tag), "w", encoding="utf8", newline="\n") as O:
    O.write("group\tbin\tn\tmeth\n")
    for (g, b), (n, k) in sorted(PROFCHK.items()):
        O.write("%s\t%d\t%d\t%d\n" % (g, b, n, k))
print("saved DISP_*_%s.tsv" % tag)
