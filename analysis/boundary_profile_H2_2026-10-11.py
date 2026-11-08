# -*- coding: utf-8 -*-
# 경계 프로파일 다시 — H2 관문 + 컷 0.508  (2026-10-11)
#
#   왜 고치나 (실측 2026-10-11)
#     ⒜ 디스크의 `RU.tsv.gz`·`CALLS.tsv.gz` 에 `data_source` 열이 없다
#        (`load_readunit_v2.py` 는 넣지만 파일은 v1 판) ⇒ 자료원을 거르지 못한 채 썼다.
#        되짚어 세니 샘플 303 중 HPRC 전용 289 · 1kGP+HPRC 이름충돌 9 ·
#        GIAB_SOMATIC(암 HG008-T) 1 · GIAB_ONT 3 · 1kGP 전용 1.
#        깨끗한 창 기준 : 대장 H2 안이 90.1% / 91.0% 뿐.
#     ⒝ `boundary_profile_wide.py` 의 컷이 0.5 인데 프로젝트 상수는 0.508 이다.
#
#   한 패스에서 다 낸다 (무거운 IO 규칙) — 788 MB 를 한 번만 읽는다
#     ① 새 자 (H2 · 0.508)  로 BOUNDARY_WIDE
#     ② 옛 자 (관문없음 · 0.5) 로도 같이 내서 재현검정을 공짜로 만든다
#     ③ 깨끗한 3+3 창 전수표 (새 자) — SuppFig S7 분자 다시 고르기용
#     ④ 대비 통과분만의 평균 — 회색이 완만한 까닭이 ㉠겹침인지 ㉡약한 창인지 가른다
#     ⑤ 창마다 전환 자리 — "AluY 4/4" 를 분포로 바꾼다
#
#   원본 `boundary_profile_wide.py` 와 그 산출물은 건드리지 않는다. 새 이름으로만 낸다.
import collections
import csv
import gzip
import io
import os
import sys

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


W = _need("RDNA_WORK_ROOT")
V = os.path.join(_need("RDNA_VERSION_ROOT"), "08_versions", "v11")
T = _need("RDNA_V11_ROOT")
RU = os.path.join(W, "out", "RU.tsv.gz")
CA = os.path.join(W, "out", "CALLS.tsv.gz")
SUM = os.path.join(V, "03_maps", "g8canonv10", "CONSENSUS_SUM_G8_full.tsv")
REG = os.path.join(T, "02_tables", "TIERS_H0H3_SLOTS_2026-10-05.tsv")
OUT = os.path.join(W, "out_H2_2026-10-11")
os.makedirs(OUT, exist_ok=True)

UL, MLTHR, MINCPG = 2240, 128, 20
NBEF, NAFT = 3, 3
NU = NBEF + NAFT
NL = chr(10)
#  자 두 개를 나란히 : (이름, 유닛 하드콜 컷, H2 관문을 거나)
RULERS = [("H2_508", 0.508, True), ("old_500", 0.500, False)]
#  대비 : 높은 쪽 3 유닛 최소 >= HI · 낮은 쪽 3 유닛 최대 <= LO
LO, HI = 0.25, 0.70
#  전환 자리 자 : 앞뒤 RW 개 CpG 중 RF 이상이 그 상태
RW, RF = 20, 0.90
EL = [("AluY", 1, 318), ("spacer-1", 319, 1178), ("GA-rich", 1179, 1551),
      ("spacer-2", 1552, 2031), ("5'-flank", 2032, 2121), ("5S gene", 2122, 2240)]


def elem(p):
    for nm, a, b in EL:
        if a <= p <= b:
            return nm
    return "?"


H2 = {(r["sample"], r["slot"]) for r in
      csv.DictReader(io.open(REG, encoding="utf8"), delimiter=chr(9)) if r["tier"] == "H2"}
sys.stderr.write("H2 슬롯 %d%s" % (len(H2), NL))

BRK = set()
for i, L in enumerate(io.open(SUM, encoding="utf8")):
    if i == 0:
        continue
    f = L.rstrip(NL).split("\t")
    if len(f) >= 8 and int(f[7]) > 0:
        BRK.add((f[0], f[1]))

#  유닛 하드콜 · beta — 자마다 따로
CLS = {k: {} for k, _, _ in RULERS}
BETA = {}
for i, L in enumerate(gzip.open(RU, "rt", encoding="utf8")):
    if i == 0:
        continue
    f = L.rstrip(NL).split("\t")
    if len(f) < 10 or int(f[7]) < MINCPG or (f[0], f[1]) in BRK:
        continue
    k3 = (f[0], f[1], f[2], int(f[3]))
    b = float(f[9])
    BETA[k3] = b
    for nm, thr, gate in RULERS:
        if gate and (f[0], f[1]) not in H2:
            continue
        CLS[nm][k3] = 1 if b >= thr else 0
for nm, _, _ in RULERS:
    sys.stderr.write("%s 유닛 %d%s" % (nm, len(CLS[nm]), NL))

WIN = {nm: {"H2L": collections.defaultdict(lambda: [0, 0]),
            "L2H": collections.defaultdict(lambda: [0, 0])} for nm, _, _ in RULERS}
#  ④ 대비 통과분만 (새 자에서만)
WINC = {"H2L": collections.defaultdict(lambda: [0, 0]),
        "L2H": collections.defaultdict(lambda: [0, 0])}
CTRL = {nm: {0: collections.defaultdict(lambda: [0, 0]),
             1: collections.defaultdict(lambda: [0, 0])} for nm, _, _ in RULERS}
N = collections.Counter()
ROWS = []        # ③ 깨끗한 창 전수 (새 자)
STEPS = []       # ⑤ 전환 자리
NEW = RULERS[0][0]


def step_pos(seq_pos, seq_m, A):
    """원시 CpG 로 전환 자리. seq_pos·seq_m 은 창 안 차례(원위→근위 걷는 차례) 그대로."""
    n = len(seq_m)
    if n < 2 * RW:
        return None
    B = 1 - A
    back = [0] * n
    fwd = [0] * n
    c = 0
    for i in range(n):
        c += 1 if seq_m[i] == A else 0
        if i >= RW:
            c -= 1 if seq_m[i - RW] == A else 0
        back[i] = c / float(min(i + 1, RW))
    c = 0
    for i in range(n - 1, -1, -1):
        c += 1 if seq_m[i] == B else 0
        if i + RW < n:
            c -= 1 if seq_m[i + RW] == B else 0
        fwd[i] = c / float(min(n - i, RW))
    la = [i for i in range(n) if back[i] >= RF]
    if not la:
        return None
    ia = max(la)
    fb = [i for i in range(ia + 1, n) if fwd[i] >= RF]
    if not fb:
        return None
    ib = min(fb)
    return ia, ib


def process(read, key3):
    for nm, thr, gate in RULERS:
        cl = CLS[nm]
        idxs = sorted(read, reverse=True)
        if len(idxs) < 3:
            continue
        c = [cl.get(key3 + (i,)) for i in idxs]
        #  대조 : 양옆이 같은 class 인 블록 내부 유닛
        for a in range(1, len(idxs) - 1):
            if idxs[a - 1] - idxs[a] != 1 or idxs[a] - idxs[a + 1] != 1:
                continue
            if c[a] is None or c[a - 1] != c[a] or c[a + 1] != c[a]:
                continue
            d = CTRL[nm][c[a]]
            for p, m in read[idxs[a]]:
                v = d[p]
                v[0] += 1
                v[1] += m
        #  전이 창
        for a in range(len(idxs) - 1):
            if c[a] is None or c[a + 1] is None or c[a] == c[a + 1]:
                continue
            lo, hi = a - NBEF + 1, a + NAFT
            if lo < 0 or hi >= len(idxs):
                continue
            if any(idxs[t] - idxs[t + 1] != 1 for t in range(lo, hi)):
                continue
            if any(c[t] != c[a] for t in range(lo, a + 1)):
                continue
            if any(c[t] != c[a + 1] for t in range(a + 1, hi + 1)):
                continue
            kind = "H2L" if c[a] == 1 else "L2H"
            N["%s|%s" % (nm, kind)] += 1
            units = [idxs[t] for t in range(lo, hi + 1)]
            bts = [BETA.get(key3 + (u,)) for u in units]
            for j, t in enumerate(range(lo, hi + 1)):
                base = j * UL
                for p, m in read[idxs[t]]:
                    v = WIN[nm][kind][base + p]
                    v[0] += 1
                    v[1] += m
            if nm != NEW:
                continue
            #  ---- 새 자에서만 : 창 표 · 대비 · 전환 자리
            ncs = [len(read[u]) for u in units]
            if None in bts:
                continue
            hi_s = min(bts[:3]) if kind == "H2L" else min(bts[3:])
            lo_s = max(bts[3:]) if kind == "H2L" else max(bts[:3])
            okc = (hi_s >= HI) and (lo_s <= LO)
            ROWS.append((kind,) + key3 + (units[0], units[-1],
                                          ",".join(str(x) for x in ncs),
                                          ",".join("%.4f" % x for x in bts),
                                          "%.4f" % (hi_s - lo_s), 1 if okc else 0))
            if okc:
                for j, t in enumerate(range(lo, hi + 1)):
                    base = j * UL
                    for p, m in read[idxs[t]]:
                        v = WINC[kind][base + p]
                        v[0] += 1
                        v[1] += m
            #  ⑤ 전환 자리 : 창을 CpG 수열로 펴고 원시 자로 가른다
            sp, sm, su = [], [], []
            for j, t in enumerate(range(lo, hi + 1)):
                for p, m in sorted(read[idxs[t]]):
                    sp.append(p)
                    sm.append(m)
                    su.append(j - (NBEF - 1))
            r = step_pos(sp, sm, c[a])
            if r:
                ia, ib = r
                STEPS.append((kind, key3[0], key3[1], su[ia], sp[ia], elem(sp[ia]),
                              su[ib], sp[ib], elem(sp[ib]), ib - ia - 1,
                              1 if okc else 0))


cur, buf, seen, skip, n = None, collections.defaultdict(list), set(), False, 0
for L in gzip.open(CA, "rt", encoding="utf8"):
    if L[0] == "s":
        continue
    f = L.rstrip(NL).split("\t")
    if len(f) < 6:
        continue
    rid = (f[0], f[1], f[2])
    if rid != cur:
        if cur is not None and not skip:
            process(buf, cur)
        buf = collections.defaultdict(list)
        cur = rid
        skip = rid in seen
        seen.add(rid)
    if not skip and ((f[0], f[1], f[2], int(f[3])) in BETA):
        buf[int(f[3])].append((int(f[4]), 1 if int(f[5]) >= MLTHR else 0))
    n += 1
    if n % 50000000 == 0:
        sys.stderr.write("  %d 백만%s" % (n // 1000000, NL))
        sys.stderr.flush()
if cur is not None and not skip:
    process(buf, cur)
print("창 수 : %s" % dict(N))


def dump(name, win, ctrl):
    o = io.open(os.path.join(OUT, name), "w", encoding="utf8", newline=NL)
    o.write("kind\tg\tunit_off\tunit_pos\tn\tbeta\tn_LL\tb_LL\tn_HH\tb_HH\tscore" + NL)
    k = 0
    for kind in ("H2L", "L2H"):
        for g in range(1, NU * UL + 1):
            v = win[kind].get(g)
            if not v or v[0] < 20:
                continue
            up = (g - 1) % UL + 1
            off = (g - 1) // UL - (NBEF - 1)
            cl, ch = ctrl[0].get(up), ctrl[1].get(up)
            if not cl or not ch or cl[0] < 50 or ch[0] < 50:
                continue
            bl, bh = cl[1] / cl[0], ch[1] / ch[0]
            b = v[1] / v[0]
            s = (b - bl) / (bh - bl) if bh - bl >= 0.15 else ""
            o.write("%s\t%d\t%d\t%d\t%d\t%.4f\t%d\t%.4f\t%d\t%.4f\t%s%s"
                    % (kind, g, off, up, v[0], b, cl[0], bl, ch[0], bh,
                       ("%.4f" % s) if s != "" else "", NL))
            k += 1
    o.close()
    print("  %-34s %s 줄" % (name, format(k, ",")))


for nm, _, _ in RULERS:
    dump("BOUNDARY_WIDE_%s.tsv" % nm, WIN[nm], CTRL[nm])
dump("BOUNDARY_WIDE_%s_대비통과만.tsv" % NEW, WINC, CTRL[NEW])

with io.open(os.path.join(OUT, "CLEAN_WINDOWS_%s.tsv" % NEW), "w", encoding="utf8",
             newline=NL) as o:
    o.write("kind\tsample\tslot\tread_id\tunit_first\tunit_last\tncpg_per_unit"
            "\tbeta_per_unit\tgap\tcontrast_ok" + NL)
    for r in ROWS:
        o.write("\t".join(str(x) for x in r) + NL)
print("  깨끗한 창 %s" % format(len(ROWS), ","))

with io.open(os.path.join(OUT, "STEP_POS_%s.tsv" % NEW), "w", encoding="utf8",
             newline=NL) as o:
    o.write("kind\tsample\tslot\tunit_off_before\tunit_pos_before\telem_before"
            "\tunit_off_after\tunit_pos_after\telem_after\tcpg_between\tcontrast_ok" + NL)
    for r in STEPS:
        o.write("\t".join(str(x) for x in r) + NL)
print("  전환 자리 %s" % format(len(STEPS), ","))
print("산출 : %s" % OUT)
