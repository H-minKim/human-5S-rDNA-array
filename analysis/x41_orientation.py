# -*- coding: utf-8 -*-
# array 방향(근위/원위)을 조립본에서 실측 — 두 갈래로 재고 서로 맞는지 본다
#   (가) 중심절 좌표 : tables/CENTROMERE.tsv — 중심절이 있는 쪽이 근위(proximal) · 반대가 원위(distal)
#
#   고침 — 그 전에는 "가까운 염색체 끝" 으로 정했고 6 줄이 틀렸다
#     옛 자 : `원위쪽 = 가까운 끝` (가정 : *"중심절 좌표 자료는 없다"*).
#     그 뒤 `tables/CENTROMERE.tsv` 가 생겼는데(2026-10-06) 이 글을 안 고쳤다.
#     `MEMORY.md` §2-17 *"가정을 세웠으면 그 가정을 잰다"* 를 가정을 적어 두고도 어긴 것이다.
#     실측(53 전수) : 맞음 46 · 어긋남 6 · 겹침 1
#       · 제브라피시 chr4 ×3(36.22·36.24·37.10 Mb) · chr18 ×1(31.47 Mb)
#         ⇒ array 가 긴 팔 에 있어 제 팔 텔로미어가 반대 팔 텔로미어보다 가까웠다
#       · 금화조 chr9 mat+pat ×2 ⇒ acrocentric — 중심절이 0.01~0.20 Mb(염색체 맨 앞)라
#         *"가까운 끝"* 이 곧 중심절 쪽 이었다
#       · 유인원 다섯은 전부 맞았다(중심절 99.8~117.4 Mb · array 16.5~21.0 Mb) ⇒ 그래서 안 들켰다
#     `끝비`(먼끝/가까운끝)는 싼 대리지표 다 — 금화조 chr9 은 `2.76` 로 강해 보였는데 틀렸다(§2-5).
#     여기에 "가까운 끝" 자를 되돌리지 말 것 — 근위/원위는 중심절 좌표로만 정한다.
#
#   (가-보조) 끝까지의 거리 : 00_catalog/5S_array_structure.tsv 의 contig_len · dist_contig_start/end
#        ⇒ `낮은끝까지`·`높은끝까지`·`가까운끝`·`가까운끝까지_Mb`·`먼끝까지_Mb` 는 기록용으로만 남긴다
#          (방향 판정에는 쓰지 않는다). 아래 글들이 그 열을 읽으므로 열 이름은 지우지 않는다.
#   (나) flank 서열 상동 : 종끼리 낮은 좌표 flank 끼리 가 낮은-높은 보다 더 닮았는지 (21-mer 자카드)
#        닮았으면 종들의 array 방향이 서로 같은 쪽으로 놓여 있다는 뜻(= 좌표 방향이 종끼리 뒤집히지 않았다)
#   산출 : tables/ORIENTATION.tsv (array 마다) · tables/ORIENTATION_FLANK_JACCARD.tsv (유인원 쌍) · 로그 logs/X41_ORIENT.txt
import io, os, sys, glob, re, itertools
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
LOG = io.open(NEW + "/logs/X41_ORIENT.txt", "w", encoding="utf8", newline="\n")
REASON = "array 방향(근위/원위) 조립본 실측(2026-10-01)"
K = 21


def P(s=""):
    print(s)
    LOG.write(s + "\n")
    LOG.flush()


src = io.open(OLD + "/code/x11_coding_readunit.py", encoding="utf8").read()
L11 = 'LOG = io.open(OUT + ("/X11_CODING_PILOT.txt" if PILOT else "/X11_CODING_BUILD.txt"), "a", encoding="utf8", newline="\\n")'
head = src.split("STAT = collections.Counter()")[0].replace(L11, "LOG = io.StringIO()")
G = {"__name__": "x11hdr"}
exec(compile(head, OLD + "/code/x11_coding_readunit.py", "exec"), G)
read_fa, COMP = G["read_fa"], G["COMP"]

ST = pd.read_csv(B + "/00_catalog/5S_array_structure.tsv", sep="\t")
MT = pd.read_csv(NEW + "/tables/PERARRAY_META.tsv", sep="\t")
SPK = dict(zip(MT.array_id, MT.종))

# ─────────────────────────── 중심절 좌표 (근위/원위를 정하는 유일한 자)
CEN = pd.read_csv(NEW + "/tables/CENTROMERE.tsv", sep="\t")
_cen_acc, _cen_sp = {}, {}
for _r in CEN.itertuples(index=False):
    if pd.isna(_r.시작) or pd.isna(_r.끝):
        continue
    _cen_acc[str(_r.접근번호)] = (float(_r.시작), float(_r.끝))
    _cen_sp[(str(_r.종), str(_r.염색체))] = (float(_r.시작), float(_r.끝))
# 접근번호 → 염색체이름 : `이름표`("chr4 · 36.2 Mb" · "chr9_mat") 와 CHROM_NAMES 둘로 짚는다
_CN = pd.read_csv(NEW + "/tables/CHROM_NAMES.tsv", sep="\t")
_chrno = dict(zip(_CN.접근번호.astype(str), _CN.염색체번호.astype(str)))
_lab = {}
for _t in MT.itertuples(index=False):
    _n = str(_t.이름표).split(" · ")[0].split("·")[0].strip()
    _lab.setdefault(str(_t.chrom), _n)


def cen_of(sp, acc):
    """Centromere interval (start, end, how it was located). If it is not found, None is returned and the script stops below rather than passing silently."""
    if acc in _cen_acc:
        return _cen_acc[acc] + ("접근번호",)
    nm = _lab.get(acc, "") or _chrno.get(acc, "")
    for cand, how in ((nm, "염색체이름"),
                      (nm.replace("_mat", "").replace("_pat", ""), "염색체이름(하플로 떼고)"),
                      (nm.replace("_pat", "_mat"), "부계→모계 상동")):
        if cand and (sp, cand) in _cen_sp:
            return _cen_sp[(sp, cand)] + (how,)
    return None


rows, missing = [], []
for t in ST.itertuples(index=False):
    if t.array_id not in SPK:
        continue
    sp = SPK[t.array_id]
    d0, d1 = float(t.dist_contig_start), float(t.dist_contig_end)
    tel = "낮은 좌표" if d0 < d1 else "높은 좌표"          # 기록용 — 방향 판정에 쓰지 않는다
    a0, a1 = float(t.array_start), float(t.array_end)
    got = cen_of(sp, str(t.chrom))
    if got is None:
        missing.append((sp, str(t.chrom), t.array_id))
        continue
    c0, c1, how = got
    # 근위 = 중심절이 있는 쪽
    prox = "높은 좌표" if (c0 + c1) / 2.0 > (a0 + a1) / 2.0 else "낮은 좌표"
    # array 가 중심절 구간과 겹치면 근위/원위가 뜻을 잃는다 ⇒ 표시해 두고 그림에서 글자를 빼게 한다
    ovl = not (a1 < c0 or a0 > c1)
    gap = 0.0 if ovl else (c0 - a1 if prox == "높은 좌표" else a0 - c1)
    rows.append(dict(종=sp, array_id=t.array_id, chrom=t.chrom, 염색체길이=int(t.contig_len),
                     array시작=int(a0), array끝=int(a1),
                     낮은끝까지=int(d0), 높은끝까지=int(d1),
                     가까운끝=tel, 원위쪽=("낮은 좌표" if prox == "높은 좌표" else "높은 좌표"), 근위쪽=prox,
                     가까운끝까지_Mb=round(min(d0, d1) / 1e6, 2), 먼끝까지_Mb=round(max(d0, d1) / 1e6, 2),
                     염색체내_상대위치=round((a0 + a1) / 2.0 / t.contig_len, 4),
                     중심절시작=int(c0), 중심절끝=int(c1),
                     중심절까지_Mb=round(gap / 1e6, 3), 중심절겹침=("예" if ovl else "아니오"),
                     근위근거="중심절 (%s)" % how,
                     끝자와같나=("같음" if prox != tel else "다름")))
assert not missing, "중심절을 못 찾은 array 가 있다 — 조용히 넘기지 않는다 : %r" % (missing,)
O = pd.DataFrame(rows)
assert len(O) == len(SPK), "array %d / 대장 %d — 수가 안 맞는다" % (len(O), len(SPK))
P("array %d · 근위(중심절) 쪽 분포 :" % len(O))
P(O.groupby(["종", "근위쪽"]).size().to_string())
P("")
P("옛 자('가까운 끝') 와 다른 줄 : %d / %d" % (int((O.끝자와같나 == "다름").sum()), len(O)))
for _r in O[O.끝자와같나 == "다름"].itertuples(index=False):
    P("  %-9s %-12s %-11s 중심절 %.2f-%.2f Mb · array %.2f-%.2f Mb ⇒ 근위 %s (옛 자 %s)"
      % (_r.종, _r.chrom, str(_r.array_id)[-11:], _r.중심절시작 / 1e6, _r.중심절끝 / 1e6,
         _r.array시작 / 1e6, _r.array끝 / 1e6, _r.근위쪽,
         "높은 좌표" if _r.가까운끝 == "낮은 좌표" else "낮은 좌표"))
P("")
P("중심절과 겹쳐 근위/원위를 적을 수 없는 array : %d" % int((O.중심절겹침 == "예").sum()))
for _r in O[O.중심절겹침 == "예"].itertuples(index=False):
    P("  %-9s %s" % (_r.종, _r.array_id))
P("")
P("종마다 가까운 끝까지 / 먼 끝까지 (Mb · 중앙값)")
P(O.groupby("종")[["가까운끝까지_Mb", "먼끝까지_Mb"]].median().to_string())

# ---------------------------------------------------------------- (나) flank 21-mer 상동 (유인원)
APE = ["gorilla_mGorGor1", "bonobo_mPanPan1", "chimp_mPanTro3", "orang_bornean_mPonPyg2", "orang_sumatran_mPonAbe1"]
FL = {}
for sp in APE:
    fs = glob.glob(B + "/01_canonical/%s/*.fa" % sp)
    fs = [f for f in fs if not f.replace("\\", "/").endswith(".unit.fa")]
    assert len(fs) == 1, (sp, len(fs))
    hd, s = read_fa(fs[0])
    m = re.search(r"flank(\d+)_array_flank(\d+) strand=([+-]) chrom=(\S+) array=(\d+)-(\d+)", hd)
    f5, f3, strand = int(m.group(1)), int(m.group(2)), m.group(3)
    if strand == "-":
        s = s.encode().translate(COMP)[::-1].decode()          # 게놈 정방향으로
        lo_len, hi_len = f3, f5
    else:
        lo_len, hi_len = f5, f3
    lo, hi = s[:lo_len].upper(), s[-hi_len:].upper()

    def kmers(x):
        out = set()
        for i in range(len(x) - K + 1):
            w = x[i:i + K]
            if "N" in w:
                continue
            rc = w.encode().translate(COMP)[::-1].decode()
            out.add(min(w, rc))                                  # 가닥 무관 정규형
        return out
    FL[sp] = dict(lo=kmers(lo), hi=kmers(hi), lo_len=lo_len, hi_len=hi_len)
    P("  %-24s 낮은쪽 flank %d bp (21-mer %s) · 높은쪽 %d bp (%s)"
      % (sp, lo_len, format(len(FL[sp]["lo"]), ","), hi_len, format(len(FL[sp]["hi"]), ",")))

jac = lambda a, b: (len(a & b) / len(a | b)) if (a or b) else np.nan
jr = []
for a, b in itertools.combinations(APE, 2):
    v = {("lo", "lo"): jac(FL[a]["lo"], FL[b]["lo"]), ("hi", "hi"): jac(FL[a]["hi"], FL[b]["hi"]),
         ("lo", "hi"): jac(FL[a]["lo"], FL[b]["hi"]), ("hi", "lo"): jac(FL[a]["hi"], FL[b]["lo"])}
    jr.append(dict(종A=a, 종B=b, 낮은_낮은=v[("lo", "lo")], 높은_높은=v[("hi", "hi")],
                   낮은_높은=v[("lo", "hi")], 높은_낮은=v[("hi", "lo")],
                   같은쪽_평균=np.mean([v[("lo", "lo")], v[("hi", "hi")]]),
                   엇갈린쪽_평균=np.mean([v[("lo", "hi")], v[("hi", "lo")]])))
J = pd.DataFrame(jr)
J["배수"] = J.같은쪽_평균 / J.엇갈린쪽_평균.replace(0, np.nan)
P("")
P("유인원 쌍마다 flank 21-mer 자카드 (같은 쪽 vs 엇갈린 쪽)")
P(J.round(4).to_string(index=False))
P("")
P("같은 쪽이 엇갈린 쪽보다 몇 배 : 중앙 %.1f 배 (최소 %.1f · 최대 %.1f) · 쌍 %d/%d 에서 같은 쪽이 큼"
  % (J.배수.median(), J.배수.min(), J.배수.max(), int((J.같은쪽_평균 > J.엇갈린쪽_평균).sum()), len(J)))

TMP = os.environ.get("TEMP", ".") + "/x41"
os.makedirs(TMP, exist_ok=True)
O.to_csv(TMP + "/ORIENTATION.tsv", sep="\t", index=False, encoding="utf8", lineterminator="\n")
J.round(6).to_csv(TMP + "/ORIENTATION_FLANK_JACCARD.tsv", sep="\t", index=False, encoding="utf8", lineterminator="\n")
arch.put(TMP + "/ORIENTATION.tsv", "tables/ORIENTATION.tsv", REASON)
arch.put(TMP + "/ORIENTATION_FLANK_JACCARD.tsv", "tables/ORIENTATION_FLANK_JACCARD.tsv", REASON)
P("")
P("표 : latest/tables/ORIENTATION.tsv · ORIENTATION_FLANK_JACCARD.tsv")
LOG.close()
