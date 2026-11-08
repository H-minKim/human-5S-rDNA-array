# -*- coding: utf-8 -*-
"""항목 1+3 — 미처리 대조를 **같은 추정량**으로 돌리고, 비의 비(ratio-of-ratios)를 낸다.

왜 비의 비인가
  미처리 시료는 같은 5S 서열이므로 창별 심문 A 수가 **정확히 같다**.
  따라서 (ICR/flank)_처리 ÷ (ICR/flank)_미처리 는 서열 조성을 **완전히** 상쇄한다.
  A 수를 모형으로 보정할 필요가 없다(모형 (1-q)^n 이 부정확한 것은 이미 확인됨).

창 (전부 45 bp 로 길이 통일 · 기존 53 bp 창도 연속성 위해 병기)
  W_prom   tx  -50..  -6   TFIIIB 대역
  W_gene5  tx   +0.. +44   유전자 5' (ICR 바깥이나 Pol III 복합체 안일 수 있음)
  W_ICR    tx  +50.. +94   A-box/IE/Box C
  W_flank  tx +250..+294   스페이서
  (구형) ICR53 +45..+97 / FLK53 +250..+302
"""
import io, os, sys, subprocess, collections, bisect, math
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace", line_buffering=True)
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


W = _need("RDNA_TF_TMP_DIR")
GENE = _need("RDNA_GENE_QUERY_FA")
IND = _need("RDNA_FS_OUT_DIR")
SUF = os.environ.get("SUF", ".hit.fq")
TH = os.environ.get("TH", "16")
THR = int(os.environ.get("THR", "128"))
MINA = int(os.environ.get("MINA", "6"))

WIN = [("W_prom", -50, -5), ("W_gene5", 0, 45), ("W_ICR", 50, 95),
       ("W_flank", 250, 295), ("ICR53", 45, 98), ("FLK53", 250, 303)]


def parse_A(seq, mm, ml):
    if not mm or not ml:
        return None
    probs = [int(x) for x in ml.split(",")[1:]]
    i0 = 0
    for b in [x for x in mm.split(";") if x.strip()]:
        f = b.split(",")
        head = f[0]
        base = head[0]
        code = head.split("+")[1].rstrip("?.") if "+" in head else ""
        d = [int(x) for x in f[1:]]
        blk = probs[i0:i0 + len(d)]
        i0 += len(d)
        if base == "A" and code == "a":
            cp = [i for i, c in enumerate(seq) if c == "A"]
            pos = {}
            kk = -1
            for x, p in zip(d, blk):
                kk += x + 1
                if kk < len(cp):
                    pos[cp[kk]] = p
            return cp, pos
    return None


def wc(allA, mods, lo, hi):
    i = bisect.bisect_left(allA, lo)
    j = bisect.bisect_left(allA, hi)
    return j - i, sum(1 for k in range(i, j) if allA[k] in mods)


if not os.path.isdir(W):
    os.makedirs(W)
OUT = {}
for smp in sys.argv[1:]:
    fq = os.path.join(IND, smp + SUF)
    if not os.path.exists(fq):
        sys.stderr.write("없음 %s\n" % fq)
        continue
    seqs = {}
    with io.open(fq, encoding="utf8", errors="replace") as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            s = fh.readline().strip()
            fh.readline(); fh.readline()
            p = h.rstrip("\n").split("\t")
            mm = next((x[5:] for x in p if x.startswith("MM:Z:")), None)
            ml = next((x[5:] for x in p if x.startswith("ML:B:")), None)
            seqs[p[0][1:].split()[0]] = (s, mm, ml)
    fa = os.path.join(W, "un_%s.fa" % smp)
    with io.open(fa, "w", newline="\n") as g:
        for r, v in seqs.items():
            g.write(">%s\n%s\n" % (r, v[0]))
    db = os.path.join(W, "undb_%s" % smp)
    subprocess.run(["makeblastdb", "-in", fa, "-dbtype", "nucl", "-out", db],
                   check=True, capture_output=True)
    r1 = subprocess.run(["blastn", "-task", "blastn", "-query", GENE, "-db", db,
                         "-outfmt", "6 sseqid sstart send", "-dust", "no", "-word_size", "9",
                         "-evalue", "1e-5", "-max_target_seqs", "100000000",
                         "-num_threads", TH], capture_output=True, text=True, check=True)
    hits = collections.defaultdict(list)
    for L in r1.stdout.splitlines():
        f = L.split("\t")
        if len(f) >= 3:
            hits[f[0]].append((int(f[1]), int(f[2])))
    rows = []
    tot_A = tot_m = 0
    for r, (sq, mm, ml) in seqs.items():
        g = hits.get(r, [])
        if len(g) < 3:
            continue
        sgn = -1 if sum(1 for a, b in g if b < a) > len(g) / 2 else 1
        five = sorted((min(a, b) - 1) if sgn == 1 else (max(a, b) - 1) for a, b in g)
        if sgn == -1:
            five = sorted(five, reverse=True)
        u = []
        for i in range(len(five) - 1):
            if 2100 <= abs(five[i + 1] - five[i]) <= 2400:
                u.append((five[i], five[i + 1]))
        if not u:
            continue
        P = parse_A(sq, mm, ml)
        if P is None:
            continue
        allA, pos = P
        mods = set(p for p, v in pos.items() if v >= THR)
        RL = len(sq)
        tot_A += len(allA); tot_m += len(mods)
        for f0, f1 in u:
            rec = {}
            ok = True
            for nm, t0, t1 in WIN:
                a = f0 + sgn * t0
                b = f0 + sgn * t1
                a, b = min(a, b), max(a, b)
                if a < 0 or b > RL:
                    ok = False
                    break
                n, m = wc(allA, mods, a, b)
                rec[nm] = (n, 1 if m == 0 else 0)
            if ok:
                rows.append(rec)
    OUT[smp] = (rows, tot_m / float(tot_A) if tot_A else 0.0)
    for e in (".ndb", ".nhr", ".nin", ".not", ".nsq", ".ntf", ".nto"):
        try:
            os.remove(db + e)
        except OSError:
            pass
    try:
        os.remove(fa)
    except OSError:
        pass
    sys.stderr.write("  %s : 유닛 %d\n" % (smp, len(rows)))
    sys.stderr.flush()


def prot(rows, nm, minA=MINA):
    v = [r for r in rows if r[nm][0] >= minA]
    if not v:
        return float("nan"), 0
    return sum(r[nm][1] for r in v) / float(len(v)), len(v)


def pair_ratio(rows, a, b, minA=MINA):
    v = [r for r in rows if r[a][0] >= minA and r[b][0] >= minA]
    if len(v) < 30:
        return float("nan"), 0, float("nan")
    pa = sum(r[a][1] for r in v) / float(len(v))
    pb = sum(r[b][1] for r in v) / float(len(v))
    bb = sum(1 for r in v if r[a][1] and not r[b][1])
    cc = sum(1 for r in v if r[b][1] and not r[a][1])
    x = ((abs(bb - cc) - 1) ** 2 / float(bb + cc)) if (bb + cc) else 0.0
    p = math.erfc(math.sqrt(x / 2.0)) if x > 0 else 1.0
    return (pa / pb) if pb else float("nan"), len(v), p


print("THR=%d · MINA=%d" % (THR, MINA))
print("\n" + "=" * 100)
print("① 시료별 기본값")
print("=" * 100)
print("  시료          유닛     A당 표지율 q    창별 평균 심문 A 수")
for smp in OUT:
    rows, q = OUT[smp]
    if not rows:
        continue
    s = " · ".join("%s %.1f" % (nm, sum(r[nm][0] for r in rows) / float(len(rows)))
                   for nm, _, _ in WIN)
    print("  %-12s %5d    %.4f      %s" % (smp, len(rows), q, s))

print("\n" + "=" * 100)
print("② 창별 보호율 (창 안 m6A 0 개인 비율)")
print("=" * 100)
hdr = "  시료         " + "".join("%-11s" % nm for nm, _, _ in WIN)
print(hdr)
for smp in OUT:
    rows, q = OUT[smp]
    if not rows:
        continue
    line = "  %-12s" % smp
    for nm, _, _ in WIN:
        p, n = prot(rows, nm)
        line += "%-11s" % ("%.3f(%d)" % (p, n) if n else "-")
    print(line)

print("\n" + "=" * 100)
print("③ 같은 유닛 짝비교 — 각 창 / W_flank")
print("=" * 100)
print("  시료         " + "".join("%-20s" % nm for nm, _, _ in WIN[:3]))
for smp in OUT:
    rows, q = OUT[smp]
    if not rows:
        continue
    line = "  %-12s" % smp
    for nm, _, _ in WIN[:3]:
        r, n, p = pair_ratio(rows, nm, "W_flank")
        line += "%-20s" % ("%.3f (n=%d,p=%.2g)" % (r, n, p) if n else "-")
    print(line)
print("  (구형 53bp) ICR53/FLK53:")
for smp in OUT:
    rows, q = OUT[smp]
    if not rows:
        continue
    r, n, p = pair_ratio(rows, "ICR53", "FLK53")
    print("     %-12s %.3f (n=%d, p=%.3g)" % (smp, r, n, p))

print("\n" + "=" * 100)
print("④ 유전자 안끼리 — W_ICR / W_gene5 (구역 교락 없음)")
print("=" * 100)
for smp in OUT:
    rows, q = OUT[smp]
    if not rows:
        continue
    r, n, p = pair_ratio(rows, "W_ICR", "W_gene5")
    print("  %-12s %.3f (n=%d, p=%.3g)" % (smp, r, n, p))

_TFOUT = os.path.join(_need("RDNA_OUT_DIR"),
                      "fsplace_tf_untreated_thr%d.tsv" % THR)
o = io.open(_TFOUT, "w",
            encoding="utf8", newline="\n")
o.write("sample\t" + "\t".join("%s_nA\t%s_prot" % (nm, nm) for nm, _, _ in WIN) + "\n")
for smp in OUT:
    for r in OUT[smp][0]:
        o.write(smp + "\t" + "\t".join("%d\t%d" % r[nm] for nm, _, _ in WIN) + "\n")
o.close()
print("\n기록: %s" % _TFOUT)
