# -*- coding: utf-8 -*-
# A test of the claim:
#   5S makes its transition by changing the proportion of bimodally separated weakly and highly methylated units,
#    whereas in DXZ4 the methylation level itself rises in steps of about 0.3 at the level of the unit.
#
# Why this is a separate script:
#   p10_19_unit_steps.py uses the per-base P10_CPG_ALL.tsv.gz,
#   which holds only the three families 5S, D4Z4_4q and D4Z4_10q (6.10, 3.52 and 2.37 million rows; the rest are 0).
#   DXZ4, RNU2 and SST1 are therefore absent from its step table. The z version, which standardises by the noise, needs per-base data and cannot be extended.
#   The absolute-size version (P10_UNIT_STEPS_ABS) is defined on the unit beta alone, so the same definition
#   can be applied to all six arrays from the per-unit table P10_FRAME_ALL.tsv.
#   Control: for 5S and D4Z4 the result is printed alongside the canonical ABS table to see whether it reproduces. The canonical 5S
#            figure is restricted to a sub-element and so cannot match exactly; that difference is reported as well.
#   In a female, pooling the two X chromosomes misleads for DXZ4, so the result is also reported per chromosome, keyed on (sample, hap).
#   The canonical tables are not touched; output goes only to 02_output_tables/array_steps_extended/ .
import io, os, sys, collections, statistics


def _need(_v):
    """The deposited copy has no default paths. Stop when a root variable is missing."""
    _x = os.environ.get(_v)
    if not _x:
        sys.exit("environment variable %s is not set; see README.md, section 'Environment variables'" % _v)
    return _x



H = _need("ARRAYMETH_HANDOFF_ROOT")
SRC = os.path.join(H, "02_output_tables", "P10")
OUT = os.path.join(H, "02_output_tables", "array_steps_extended")
os.makedirs(OUT, exist_ok=True)
FAMS = ("5S", "D4Z4_4q", "D4Z4_10q", "DXZ4", "RNU2", "SST1")
ZONE = {"5S": 25, "D4Z4_4q": 10, "D4Z4_10q": 10, "DXZ4": 10, "RNU2": 10, "SST1": 10}
DIRN = {"5S": ("distal",), "D4Z4_4q": ("proximal",), "D4Z4_10q": ("proximal",),
        "DXZ4": ("distal", "proximal"), "RNU2": ("distal", "proximal"), "SST1": ("distal", "proximal")}
MINCPG = 10

MINC, COVMIN = 20, 0.9        # complete units, by the same rule as fig39, fig43 and fig47
_TOTB = collections.defaultdict(collections.Counter)
with io.open(os.path.join(SRC, "P10_FRAME_ALL.tsv"), encoding="utf8") as _f:
    _h = _f.readline().rstrip(chr(10)).split(chr(9))
    _j = {k: _h.index(k) for k in _h}
    for _L in _f:
        _g = _L.rstrip(chr(10)).split(chr(9))
        if _g[_j['family']] not in FAMS:
            continue
        for _t in _g[_j['bins']].split(';'):
            if _t:
                _TOTB[_g[_j['family']]][int(_t.split(':')[0])] += 1
CPGB = {fam: {b for b, v in c.items() if v >= 0.10 * max(c.values())}
        for fam, c in _TOTB.items()}


def _complete(g, ix):
    fam = g[ix['family']]
    if int(g[ix['n_cpg']]) < MINC:
        return False
    bs = {int(t.split(':')[0]) for t in g[ix['bins']].split(';') if t}
    cp = CPGB.get(fam) or set()
    return bool(cp) and len(bs & cp) / float(len(cp)) >= COVMIN

# ---------- unit beta, per read, which is the unit of the original definition: pairs of neighbouring units within a read
SER = collections.defaultdict(dict)       # (fam, sample, hap, read) -> {(dirn, p): beta}
CH = {}
with io.open(os.path.join(SRC, "P10_FRAME_ALL.tsv"), encoding="utf8") as f:
    hdr = f.readline().rstrip("\n").split("\t")
    ix = {k: hdr.index(k) for k in hdr}
    for L in f:
        g = L.rstrip("\n").split("\t")
        fam = g[ix["family"]]
        if fam not in FAMS:
            continue
        u, nu = g[ix["unit_from_distal"]], g[ix["n_units"]]
        if not u.lstrip("-").isdigit() or not nu.isdigit():
            continue
        nc = int(g[ix["n_cpg"]])
        if not _complete(g, ix):      # complete units only
            continue
        d, N = int(u), int(nu)
        if d > N:
            continue
        k = (fam, g[ix["sample"]], g[ix["hap"]], g[ix["read"]])
        b = int(g[ix["n_meth"]]) / float(nc)
        SER[k][("distal", d)] = b
        SER[k][("proximal", N + 1 - d)] = b
        CH[k] = (g[ix["sample"]], g[ix["hap"]])


def per_read(fam, dirn, hi, chrom=None):
    out = []
    for k, d in SER.items():
        if k[0] != fam or (chrom is not None and CH[k] != chrom):
            continue
        ps = sorted(p for (dd, p) in d if dd == dirn and 1 <= p <= hi)
        rs = [d[(dirn, b)] - d[(dirn, a)] for a, b in zip(ps, ps[1:]) if b == a + 1]
        if rs:
            out.append(rs)
    return out


def absrow(pr):
    """The original definition : only reads with a total rise >0.15, and the rise is split into three step-size bins"""
    tot, mx, nup, a, b, c = [], [], [], 0.0, 0.0, 0.0
    for rs in pr:
        t = sum(rs)
        if t <= 0.15:
            continue
        tot.append(t); mx.append(max(rs)); nup.append(sum(1 for x in rs if x >= 0.1))
        a += sum(x for x in rs if abs(x) < 0.1)
        b += sum(x for x in rs if 0.1 <= abs(x) < 0.3)
        c += sum(x for x in rs if abs(x) >= 0.3)
    s = a + b + c
    if not tot or s == 0:
        return None
    return dict(n=len(tot), n_all=len(pr), med_total=statistics.median(tot), med_max=statistics.median(mx),
                f_lt=a / s, f_mid=b / s, f_ge=c / s, med_nup=statistics.median(nup))


CANON = {("5S", "distal"): (0.735, 0.661, 0.910, "핵 하위원소 제한판"),
         ("D4Z4_4q", "proximal"): (0.323, 0.193, 0.152, "유닛 전체"),
         ("D4Z4_10q", "proximal"): (0.278, 0.161, 0.058, "유닛 전체")}
LINES = ["family\tdirection\tchromosome\tn_rising_reads\tn_reads_all\tmed_total_rise\tmed_max_step"
         "\tfrac_rise_lt0.1\tfrac_rise_0.1_0.3\tfrac_rise_ge0.3\tmed_n_up_steps_ge0.1"]
print("read 수 %s · 여섯 array 공통 정의(유닛 beta · 완전유닛 n_cpg>=%d + 흔한칸 %.0f%%)"
      % (format(len(SER), ","), MINC, 100 * COVMIN))
print()
print("%-9s %-9s %-13s %6s %8s %9s %8s %8s %8s %6s  %s" %
      ("family", "방향", "염색체", "오름", "전체", "총상승", "최대걸음", "<0.1", ".1-.3", ">=.3", "정본 대조"))


def emit(fam, dirn, lab, chrom=None):
    R = absrow(per_read(fam, dirn, ZONE[fam], chrom))
    if not R:
        print("%-9s %-9s %-13s   (오르는 read 없음)" % (fam, dirn, lab))
        return
    cmp_ = ""
    if lab == "ALL" and (fam, dirn) in CANON:
        t, m, g, note = CANON[(fam, dirn)]
        cmp_ = "정본 %.3f/%.3f/%.3f (%s)" % (t, m, g, note)
    print("%-9s %-9s %-13s %6d %8d %9.3f %9.3f %8.3f %8.3f %8.3f %6.0f  %s" %
          (fam, dirn, lab, R["n"], R["n_all"], R["med_total"], R["med_max"],
           R["f_lt"], R["f_mid"], R["f_ge"], R["med_nup"], cmp_))
    LINES.append("%s\t%s\t%s\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.1f"
                 % (fam, dirn, lab, R["n"], R["n_all"], R["med_total"], R["med_max"],
                    R["f_lt"], R["f_mid"], R["f_ge"], R["med_nup"]))


for fam in FAMS:
    for dirn in DIRN[fam]:
        emit(fam, dirn, "ALL")
for ch in sorted({CH[k] for k in CH if k[0] == "DXZ4"}):
    for dirn in DIRN["DXZ4"]:
        emit("DXZ4", dirn, "%s_%s" % ch, ch)
io.open(os.path.join(OUT, "STEPS_UNITLEVEL_ALL6.tsv"), "w", encoding="utf8",
        newline="\n").write("\n".join(LINES) + "\n")

# ---------- 2) the form of the transition: a changing proportion of two states, or a shift in level?
#   Within one read, over the rising stretch, does the unit beta take only the two values 0.15 and 0.86, or does it take intermediate ones?
print()
print("전이 구간 유닛 beta 가 두 값만 쓰는가 (중간칸 0.35-0.65 의 몫)")
print("%-9s %-9s %8s %8s %8s %8s %8s" % ("family", "방향", "유닛", "<=0.25", "중간", ">=0.75", "골(.45-.55)"))
L3 = ["family\tdirection\tn_unit\tfrac_le025\tfrac_mid\tfrac_ge075\tfrac_valley\thist20"]
for fam in FAMS:
    for dirn in DIRN[fam]:
        V = []
        for k, d in SER.items():
            if k[0] != fam:
                continue
            V += [d[(dd, p)] for (dd, p) in d if dd == dirn and 1 <= p <= ZONE[fam]]
        if not V:
            continue
        n = float(len(V))
        lo = sum(1 for b in V if b <= 0.25) / n
        mid = sum(1 for b in V if 0.35 < b < 0.65) / n
        hi = sum(1 for b in V if b >= 0.75) / n
        vl = sum(1 for b in V if 0.45 <= b <= 0.55) / n
        hist = [0] * 20
        for b in V:
            hist[min(int(b * 20), 19)] += 1
        print("%-9s %-9s %8s %8.3f %8.3f %8.3f %8.3f" % (fam, dirn, format(int(n), ","), lo, mid, hi, vl))
        L3.append("%s\t%s\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%s"
                  % (fam, dirn, int(n), lo, mid, hi, vl, ",".join("%.4f" % (x / n) for x in hist)))
io.open(os.path.join(OUT, "TRANSITION_ZONE_UNITDIST.tsv"), "w", encoding="utf8",
        newline="\n").write("\n".join(L3) + "\n")

# ---------- 3) DXZ4 per chromosome, without pooling a female's two X chromosomes
UB = collections.defaultdict(list)
for k, d in SER.items():
    if k[0] != "DXZ4":
        continue
    UB[CH[k]] += [d[(dd, p)] for (dd, p) in d if dd == "distal"]
L2 = ["sample\thap\tn_unit\tmean_beta\tfrac_le03\tfrac_ge07\tfrac_mid\thist20"]
print()
print("DXZ4 염색체별 — 합치면 가짜 봉우리가 생긴다")
print("%-10s %4s %7s %9s %8s %8s %8s" % ("sample", "hap", "n_unit", "평균beta", "<=0.3", "중간", ">=0.7"))
for ch in sorted(UB):
    v = UB[ch]; n = float(len(v))
    hist = [0] * 20
    for b in v:
        hist[min(int(b * 20), 19)] += 1
    lo = sum(1 for b in v if b <= 0.3) / n
    mid = sum(1 for b in v if 0.35 < b < 0.65) / n
    hi = sum(1 for b in v if b >= 0.7) / n
    print("%-10s %4s %7d %9.3f %8.3f %8.3f %8.3f" % (ch[0], ch[1], int(n), statistics.mean(v), lo, mid, hi))
    L2.append("%s\t%s\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%s"
              % (ch[0], ch[1], int(n), statistics.mean(v), lo, mid, hi,
                 ",".join("%.4f" % (x / n) for x in hist)))
io.open(os.path.join(OUT, "DXZ4_BY_CHROMOSOME.tsv"), "w", encoding="utf8",
        newline="\n").write("\n".join(L2) + "\n")
print()
print("쓴 표 3 개 : %s" % OUT)
