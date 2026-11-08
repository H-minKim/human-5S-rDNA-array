# -*- coding: utf-8 -*-
# The two tumour `HG008-T` slots are excluded on the array side of panels f and g.
#
#  f 는 `dispersed_vs_array_cut0508\DISP_PROFILES.tsv` 의 `array_v11 / ALL / low|high` 행을
#    그대로 그린다 ⇒ 그 두 행을 다시 계산해 갈아끼운다.
#  흩어진 쪽은 고칠 것이 없다 — `DISP_SITE_SAMPLE_ALL.tsv` 의 개체 307 에 `HG008-T` 가 없다
#    (있는 것은 `HG00867` 로 다른 정상 개체다 · 2026-10-16 확인).
#
#  한 패스로 둘을 낸다 : 전수(지금 표)와 `HG008-T` 만 따로 ⇒ 뺀 판 = 전수 − HG008T.
#    ⇒ 같은 패스에서 지금 표를 재현하는지도 함께 검산한다(재현되면 뺀 판도 믿을 수 있다).
#
#  관문(= f 의 계보 `disp_02_array_ref.py`) : ⒜ read id 중복 제거 ⒝ 사슬 끊긴 슬롯 제외
#    ⒞ 상태창 `spec` CpG >= 10 아니면 state=NA · 창 319-1178 + 1552-2240 · 컷 0.508
#    자 = 하드콜(ML >= 128) · 칸 = (unit_pos - 1) // 50
import collections
import gzip
import hashlib
import io
import os
import sys
import time

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


W = os.path.join(_need("RDNA_WORK_ROOT"))
V = os.path.join(_need("RDNA_VERSION_ROOT"), "08_versions", "v11")
SUM = os.path.join(V, "03_maps", "g8canonv10", "CONSENSUS_SUM_G8_full.tsv")
CA = os.path.join(W, "out", "CALLS.tsv.gz")
DV = os.path.join(W, "dispersed_vs_array_cut0508")
PROF = os.path.join(DV, "DISP_PROFILES.tsv")
OUTD = os.path.join(_need("RDNA_OUT_ROOT"), "02_output_tables", "Fig4_2026-10-09")
OUT = os.path.join(OUTD, "ARRAYPROF_noHG008T_2026-10-16.tsv")
LOG = os.path.join(OUTD, "ARRAYPROF_noHG008T_LOG_2026-10-16.txt")
NL = chr(10)
TAB = chr(9)
CUT = 0.508
TUM = "HG008-T"
EL = [("AluY", 1, 318), ("spacer-1", 319, 1178), ("GA-rich", 1179, 1551),
      ("spacer-2", 1552, 2031), ("5'-flank", 2032, 2121), ("5S gene", 2122, 2240)]
rep = []


def p(s=""):
    print(s)
    sys.stdout.flush()
    rep.append(s)


def spec(u):
    return 319 <= u <= 1178 or 1552 <= u <= 2240


def elem(mid):
    for n, a, z in EL:
        if a <= mid <= z:
            return n
    return "?"


BRK = set()
for i, L in enumerate(io.open(SUM, encoding="utf8")):
    if i == 0:
        continue
    f = L.rstrip(NL).split(TAB)
    if len(f) >= 8 and int(f[7]) > 0:
        BRK.add((f[0], f[1]))

#  (variant, state, bin) -> [n, meth] · variant : ALL / TUM
BINS = collections.defaultdict(lambda: [0, 0])
haps = {"ALL": set(), "TUM": set()}
inds = {"ALL": set(), "TUM": set()}
units = {"ALL": 0, "TUM": 0}
na = {"ALL": 0, "TUM": 0}
t0 = time.time()
cur, seen, skip, tum = None, set(), False, False
buf = collections.defaultdict(list)


def flush():
    for _pi, calls in buf.items():
        cs = [m for u, m in calls if spec(u)]
        for vr in (("ALL", "TUM") if tum else ("ALL",)):
            units[vr] += 1
        if len(cs) < 10:
            st = "NA"
            for vr in (("ALL", "TUM") if tum else ("ALL",)):
                na[vr] += 1
        else:
            st = "low" if (sum(cs) / float(len(cs))) < CUT else "high"
        for u, m in calls:
            b = (u - 1) // 50
            for vr in (("ALL", "TUM") if tum else ("ALL",)):
                v = BINS[(vr, st, b)]
                v[0] += 1
                v[1] += m


n = 0
for L in gzip.open(CA, "rt", encoding="utf8"):
    if L[0] == "s":
        continue
    f = L.rstrip(NL).split(TAB)
    if len(f) < 6:
        continue
    rid = (f[0], f[1], f[2])
    if rid != cur:
        if cur is not None and not skip:
            flush()
        buf = collections.defaultdict(list)
        cur = rid
        sl = (f[0], f[1])
        skip = (rid in seen) or (sl in BRK)
        tum = f[0] == TUM
        seen.add(rid)
        if not skip:
            for vr in (("ALL", "TUM") if tum else ("ALL",)):
                haps[vr].add(sl)
                inds[vr].add(f[0])
    if not skip:
        u = int(f[4])
        if 1 <= u <= 2240:
            buf[int(f[3])].append((u, 1 if int(f[5]) >= 128 else 0))
    n += 1
if cur is not None and not skip:
    flush()


def beta(vr, st, b):
    nn, mm = BINS[(vr, st, b)]
    return (nn, mm, (mm / float(nn)) if nn else None)


#  ── 산출 : 전수 · 종양만 · 뺀 판
rows = []
for st in ("low", "high", "NA"):
    for b in range(45):
        an, am, _ab = beta("ALL", st, b)
        tn, tm, _tb = beta("TUM", st, b)
        en, em = an - tn, am - tm
        if an == 0:
            continue
        rows.append([st, b, b * 50 + 1, b * 50 + 25, elem(b * 50 + 25),
                     an, am, round(am / an, 6),
                     tn, tm, (round(tm / tn, 6) if tn else ""),
                     en, em, (round(em / en, 6) if en else "")])
HDR = ["state", "bin", "bin_start", "bin_mid", "element",
       "n_all", "meth_all", "b_all", "n_tum", "meth_tum", "b_tum",
       "n_excl", "meth_excl", "b_excl"]
io.open(OUT, "w", encoding="utf8", newline=NL).write(
    NL.join([TAB.join(HDR)] + [TAB.join(str(x) for x in r) for r in rows]) + NL)

p("# 패널 f·g 의 array 프로필 — `HG008-T` 를 뺀 판 (2026-10-16)")
p()
p("입력 `out/CALLS.tsv.gz` %s 줄 · %.1f 분" % (format(n, ","), (time.time() - t0) / 60.0))
p()
p("## 1. 코호트")
p()
p("| | 하플로 | 개체 | read×유닛 | `state=NA` |")
p("|---|--:|--:|--:|--:|")
p("| 전수 (지금 그림) | **%d** | **%d** | %s | %s (%.2f%%) |"
  % (len(haps["ALL"]), len(inds["ALL"]), format(units["ALL"], ","),
     format(na["ALL"], ","), 100.0 * na["ALL"] / max(1, units["ALL"])))
p("| `HG008-T` 만 | **%d** | **%d** | %s | %s |"
  % (len(haps["TUM"]), len(inds["TUM"]), format(units["TUM"], ","), format(na["TUM"], ",")))
p("| **뺀 판 (새 그림)** | **%d** | **%d** | %s | %s |"
  % (len(haps["ALL"]) - len(haps["TUM"]), len(inds["ALL"]) - len(inds["TUM"]),
     format(units["ALL"] - units["TUM"], ","), format(na["ALL"] - na["TUM"], ",")))
p()

#  ── 지금 표를 재현하는지
PR = [l.rstrip(NL).split(TAB) for l in io.open(PROF, encoding="utf8") if l.strip()]
hh = PR[0]
bcol = [c for c in hh if c.startswith("b") and c[1:].isdigit()]
bcol.sort(key=lambda c: int(c[1:]))
cur_tab = {}
for r in PR[1:]:
    if r[hh.index("kind")] == "array_v11" and r[hh.index("name")] == "ALL":
        cur_tab[r[hh.index("state")]] = [r[hh.index(c)] for c in bcol]
p("## 2. 지금 표(`DISP_PROFILES.tsv` · `array_v11/ALL`)를 재현하나")
p()
p("| state | 칸 | 최대 |차| | 어느 칸 |")
p("|---|--:|--:|--:|")
ok = True
for st in ("low", "high"):
    if st not in cur_tab:
        continue
    worst, wb, cnt = 0.0, -1, 0
    for b in range(min(45, len(cur_tab[st]))):
        try:
            old = float(cur_tab[st][b])
        except ValueError:
            continue
        _n, _m, new = beta("ALL", st, b)
        if new is None:
            continue
        cnt += 1
        d = abs(new - old)
        if d > worst:
            worst, wb = d, b
    ok = ok and worst < 1e-4
    p("| %s | %d | **%.6f** | %d |" % (st, cnt, worst, wb))
p()
p("⇒ **%s**" % ("재현한다 (최대 |차| < 0.0001) — 뺀 판도 같은 코드로 낸 것이니 믿을 수 있다"
                if ok else "재현하지 않는다 — 뺀 판을 그림에 넣기 전에 원인을 찾아야 한다"))
p()
p("## 3. 종양을 빼면 얼마나 움직이나 (구획별)")
p()
p("| 구획 | state | 전수 | 뺀 판 | 차 |")
p("|---|---|--:|--:|--:|")
mx = 0.0
for nm, a, z in EL:
    for st in ("low", "high"):
        an = am = en = em = 0
        for b in range(45):
            if not (a <= b * 50 + 25 <= z):
                continue
            x = BINS[("ALL", st, b)]
            y = BINS[("TUM", st, b)]
            an += x[0]
            am += x[1]
            en += x[0] - y[0]
            em += x[1] - y[1]
        if an and en:
            d = em / en - am / an
            mx = max(mx, abs(d))
            p("| `%s` | %s | %.4f | **%.4f** | %+.4f |" % (nm, st, am / an, em / en, d))
p()
p("**가장 큰 움직임 %.4f**" % mx)
p()
p("표 `%s` md5 %s · %d 행"
  % (os.path.basename(OUT), hashlib.md5(io.open(OUT, "rb").read()).hexdigest()[:12], len(rows)))
io.open(LOG, "w", encoding="utf8", newline=NL).write(NL.join(rep) + NL)
p("기록 %s" % LOG)
