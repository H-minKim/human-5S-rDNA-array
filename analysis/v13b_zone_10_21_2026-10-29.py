# -*- coding: utf-8 -*-
# 유닛 안 자리별 메틸 — 구역(근위 / 내부 / 원위)별 로 쌓는다  2026-10-18
#
#   왜 새로 내나
#     전 판 (A)는 유닛 상태(저메틸/고메틸) 로 나눈 그림이었다. 구역으로 나누는 쪽 으로 갈았다.
#     있는 표로는 안 된다 :
#       · `UNIT_PROFILE_50bp_H2_FOR_FIG3_2026-10-16.tsv` (Fig 4 묶음) — 상태 3 × 칸 45 · 구역이 없다
#       · `02_tables/META_PROFILE.tsv` (v4) — 구역이 있으나 구역 정의가 다르다
#         (근위 ap<=15 · 원위 ad<=15 · 가운데 양쪽 >=30 · 나머지 other ⇒ 전수를 나누지 못한다)
#     ⇒ 이 그림의 다른 패널과 같은 공통 구역 으로 다시 쌓는다.
#
#   구역 (공통 5 구역을 쪽으로 묶은 것 — 새 경계가 아니다)
#     근위 = ap <= ad 이고 ap <= 17   (= 공통 '근위 1-3' + '근위 4-17')
#     원위 = ad <  ap 이고 ad <= 22   (= 공통 '원위 11-22' + '원위 1-10')
#     내부 = 나머지                   (= 공통 '내부')          셋이 전수를 가른다
#     (ad = array_pos · ap = N - array_pos + 1)
#
#   자·관문 — 다른 패널과 같게
#     하드콜 ML >= 128 · read×유닛마다 정본 CpG 호출 150 개 이상이면서 덮임 0.9 이상
#     게다가 이 그림의 다른 패널과 코호트를 맞추려고 L2_HAP_UNIT 의 유닛(n_read>=5) 45,544 개 로 한정한다.
#     칸 값 = Σ(메틸 호출) / Σ(호출)  — read × 유닛 을 모은 것이고 하플로타입 평균이 아니다.
#
#   되는 법 : python v13_unitprofile_zone_2026-10-18.py   (약 4 분 · raw 6 파일 순서대로)
import io, os, re, gzip, glob, time, collections, hashlib
import numpy as np
import pandas as pd
# 자료 뿌리는 환경변수 하나로 (2026-10-21 · 저장소 공개 준비)
#   안 주면 지금까지 쓰던 자리를 그대로 쓴다 — 돌아가는 것이 달라지지 않는다.
# 예치표의 열 이름과 값을 영문으로 낸다
#   쓰기 바로 전 한 자리 에서만 바꾼다 — 이 글 안의 내부 이름은 그대로 둔다.
_KO2EN_COL = {
    "FIG3_A_UNITPROFILE_ZONE.tsv": {"구역": "zone", "유닛안자리": "pos_in_unit",
                                    "유전자기준": "pos_gene_ref", "호출수": "n_call",
                                    "메틸호출수": "n_meth_call"},
    "FIG3_B_POSBETA.tsv": {"층": "stratum", "축": "axis", "자리": "pos"},
    "FIG3_C_LOWCOUNT.tsv": {"컷": "cut", "층": "stratum", "N중앙": "N_median",
                            "저메틸수": "n_low_meth", "저메틸수_sd": "n_low_meth_sd",
                            "비례했다면": "expected_if_proportional",
                            "저메틸몫": "low_meth_fraction"},
    "SUPP3_2_POSSD.tsv": {"분자수": "n_molecules", "평균_b_call": "mean_b_call",
                          "분자사이_sd": "sd_between_molecules"},
    "SUPP3_B_HAP3.tsv": {"하플로타입": "haplotype", "길이층": "length_class",
                         "유닛수": "n_units", "분자순번": "molecule_index",
                         "유닛자리": "unit_position"},
}
_KO2EN_VAL = {
    "FIG3_A_UNITPROFILE_ZONE.tsv": ("구역", {"근위": "proximal", "내부": "interior",
                                            "원위": "distal", "전체": "all"}),
    "FIG3_B_POSBETA.tsv": ("층", {"전체": "all"}),
    "FIG3_C_LOWCOUNT.tsv": ("컷", {"0.508(여기서 잴 골)": "0.508 (trough measured here)"}),
    "SUPP3_B_HAP3.tsv": ("길이층", {"짧은 array": "short array", "중간 길이": "medium",
                                   "긴 array": "long array"}),
}


def _to_english(n, df):
    b = os.path.basename(n)
    v = _KO2EN_VAL.get(b)
    if v and v[0] in df.columns:
        df = df.copy()
        df[v[0]] = df[v[0]].map(lambda x: v[1].get(x, x))
    c = _KO2EN_COL.get(b)
    if c:
        df = df.rename(columns=c)
    return df

DATA_ROOT = os.environ["DATA_ROOT"]
V_ROOT = os.environ["V_ROOT"]
CANON_FA = os.environ["CANON_FA"]

V = V_ROOT
OUT = DATA_ROOT
CANON = CANON_FA
MAP = os.path.join(V, r"03_maps\g8canonv10\CONSENSUS_MAP_G8_full.tsv")
RAW = os.path.join(V, r"08_methylation\raw")
DST = os.path.join(OUT, r"03_figures\paper_Fig\FIG3_A_UNITPROFILE_ZONE.tsv")
LOG = io.open(os.path.join(OUT, r"04_logs\V13_UNITPROFILE_ZONE.txt"), "w",
              encoding="utf8", newline="\n")
# 바로 위 두 줄 아래에 박혀 있던 임시폴더 경로 두 줄을 지웠다
#   지운 것 : `DST` 와 `LOG` 를 <임시폴더> 의 파일로 박아 눌 두 줄
#   까닭 : 그 두 줄이 위의 올바른 `DST`·`LOG` 를 덮어써 남이 돌리면 발행 표가 안 나왔고,
#   닫음 : 지운 판으로 다시 돌려 `FIG3_A_UNITPROFILE_ZONE.tsv` 를 바이트까지 재현 했다.
MLCUT = 128


def P(s=""):
    try:
        print(s)
    except UnicodeEncodeError:
        print(s.encode("ascii", "replace").decode())
    LOG.write(s + "\n")
    LOG.flush()


# ---------- 정본 CpG ----------
canon = "".join(x.strip() for x in io.open(CANON, encoding="utf8") if x[0] != ">")
CPG = sorted(i + 1 for i in range(len(canon) - 1) if canon[i:i + 2] == "CG")
CPGSET = set(CPG)

# ---------- 이 그림이 쓰는 유닛 집합 (다른 패널과 같은 코호트) ----------
L2 = pd.read_parquet(os.path.join(OUT, r"02_tables\L2_HAP_UNIT.parquet"),
                     columns=["sample", "slot", "ap", "ad", "n_read"])
L2 = L2[L2.n_read >= 5]
KEEP = set(zip(L2["sample"].astype(str), L2.slot.astype(str), L2.ad.astype(int)))
NHAP = L2.assign(h=L2["sample"].astype(str) + "|" + L2.slot.astype(str)).h.nunique()

# ---------- 유닛 좌표 (varlist 로 자리 옮김) ----------
IN = re.compile(r"^(\d+)I(\d+)$")
DE = re.compile(r"^(\d+)D(\d+)$")
UMAP, UN = {}, {}
with io.open(MAP, encoding="utf8") as h:
    hd = h.readline().rstrip("\n").split("\t")
    ix = {c: i for i, c in enumerate(hd)}
    for line in h:
        p = line.rstrip("\n").split("\t")
        if p[ix["seg_type"]] != "UNIT":
            continue
        smp = p[ix["sample"]]; sl = p[ix["slot"]]
        ap = int(p[ix["array_pos"]]); ul = int(p[ix["unit_len"]])
        UN[(smp, sl)] = int(p[ix["N"]])
        vl = p[ix["varlist"]]
        ins, dele = {}, {}
        if vl and vl != "-":
            for v in vl.split(","):
                m = IN.match(v)
                if m:
                    ins[int(m.group(1))] = int(m.group(2)); continue
                m = DE.match(v)
                if m:
                    dele[int(m.group(1))] = int(m.group(2))
        arr = np.zeros(ul + 2, dtype=np.int16)
        sp = 1; cp = 1; skip = 0
        while cp <= 2240 and sp <= ul:
            if skip > 0:
                skip -= 1; cp += 1; continue
            if cp in dele:
                skip = dele[cp] - 1; cp += 1; continue
            arr[sp] = cp; sp += 1
            if cp in ins:
                for _ in range(ins[cp]):
                    if sp <= ul:
                        arr[sp] = 0; sp += 1
            cp += 1
        UMAP[(smp, sl, ap)] = arr

P("=" * 104)
P("유닛 안 자리별 메틸 — 구역별 (Fig3 (A) 새 판)")
P("=" * 104)
P()
P("  정본 CpG %d · MLCUT %d" % (len(CPG), MLCUT))
P("  구역 : 근위(ap<=ad 이고 ap<=10) · 내부(나머지) · 원위(ad<ap 이고 ad<=21)")
P("  코호트 : L2_HAP_UNIT 의 n_read>=5 유닛 %s 개 · 하플로타입 %d" % (format(len(KEEP), ","), NHAP))
P()

ACC = collections.defaultdict(lambda: [0, 0])       # (구역, cpg) -> [n, n_meth]
NU = collections.Counter()                          # 구역별 read×유닛 수
USEEN = collections.defaultdict(set)                # 구역별 유닛 (중복 없이)
n_gate_out = 0
files = sorted(glob.glob(os.path.join(RAW, "METH_READUNIT_*.tsv.gz")))
for f in files:
    t1, nk = time.time(), 0
    with gzip.open(f, "rt", encoding="utf8") as h:
        h.readline()
        for line in h:
            p = line.rstrip("\n").split("\t")
            if p[1] != "HPRC_ONT_submissions":
                continue
            smp, sl = p[0], p[4]
            ad = int(p[5])                           # array_pos = 원위 번호
            if (smp, sl, ad) not in KEEP:            # 다른 패널과 같은 유닛만
                continue
            cs = p[11]
            if not cs or cs == ".":
                continue
            arr = UMAP.get((smp, sl, ad))
            if arr is None:
                continue
            N = UN.get((smp, sl), 0)
            if N <= 0:
                continue
            ap = N - ad + 1                          # 근위 번호
            L = len(arr) - 1
            hit, nraw = [], 0
            for it in cs.split(","):
                a, v = it.split(":")
                a = int(a); v = int(v)
                nraw += 1
                if a < 1 or a > L:
                    continue
                cp = int(arr[a])
                if cp == 0 or cp not in CPGSET:
                    continue
                hit.append((cp, 1 if v >= MLCUT else 0))
            if len(hit) < 150 or len(hit) / max(nraw, 1) < 0.9:
                n_gate_out += 1
                continue
            if ap <= ad and ap <= 10:
                zone = "근위"
            elif ad < ap and ad <= 21:
                zone = "원위"
            else:
                zone = "내부"
            NU[zone] += 1
            USEEN[zone].add((smp, sl, ad))
            for cp, m in hit:
                a_ = ACC[(zone, cp)]
                a_[0] += 1; a_[1] += m
                a2 = ACC[("전체", cp)]
                a2[0] += 1; a2[1] += m
            nk += 1
    P("  %-28s read×유닛 %9s · %.0f 초"
      % (os.path.basename(f)[14:-7], format(nk, ","), time.time() - t1))

rows = [(z, cp, cp - 2121, n, m, m / n) for (z, cp), (n, m) in ACC.items()]
D = _to_english("FIG3_A_UNITPROFILE_ZONE.tsv",
                pd.DataFrame(rows, columns=["구역", "유닛안자리", "유전자기준",
                                            "호출수", "메틸호출수", "beta"]))
D = D.sort_values(["zone", "pos_in_unit"])
D.to_csv(DST, sep="\t", index=False, float_format="%.5f")

P()
P("  행 %d (구역 4 × 자리 %d) · 저장 %s" % (len(D), D.유닛안자리.nunique(), os.path.basename(DST)))
P("  md5 %s" % hashlib.md5(io.open(DST, "rb").read()).hexdigest()[:12])
P()
P("  구역별 read×유닛 수 · 유닛 수")
tot_u = 0
for z in ("근위", "내부", "원위"):
    P("    %-4s read×유닛 %9s · 유닛 %7s" % (z, format(NU[z], ","), format(len(USEEN[z]), ",")))
    tot_u += len(USEEN[z])
P("    %-4s %s (= %s + %s + %s · 45,544 와 맞아야 한다)"
  % ("합", format(tot_u, ","), *[format(len(USEEN[z]), ",") for z in ("근위", "내부", "원위")]))
P("  관문에서 떨어진 read×유닛 %s" % format(n_gate_out, ","))
P()
P("  구역별 유닛 전체 평균 β (자리를 모두 합친 것)")
for z in ("전체", "근위", "내부", "원위"):
    s = D[D.구역 == z]
    P("    %-4s β %.4f (호출 %s)" % (z, s.메틸호출수.sum() / s.호출수.sum(), format(int(s.호출수.sum()), ",")))
