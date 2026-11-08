# -*- coding: utf-8 -*-
# v11 자료로 read x 유닛 메틸 표를 다시 만든다 (옛 2026-09-09_readvar 의 ru/ru_cmp 대체)
#
# 원천 : $RDNA_DATA_ROOT\for_paper\08_versions\v11\08_methylation\raw\METH_READUNIT_*.tsv.gz  (6 파일 · 1,019,408 행 · as-ran value in ENV_RDNA.tsv)
#   컬럼 : sample data_source cfg kind slot array_pos read_id strand n_cpg n_meth beta calls
#   calls = "유닛내좌표:ML값(0-255)" 목록
#
# 좌표틀 함정 (이번에 실측으로 잡은 것)
#   calls 의 유닛내좌표는 그 슬롯 자신의 유닛 프레임이다 (정본 2,240 이 아님).
#   유닛 길이가 슬롯마다 다르다(2,226~2,248 …). 그냥 정본 CpG 와 맞추면 적중률 0.27 뿐이고
#   한 자리 밀어도 0.44 다. 자리이동 하나로는 절대 안 맞는다 — 유닛 안에서 indel 이 쌓이기 때문.
#   ⇒ 지도의 varlist 로 유닛마다 슬롯좌표 -> 정본좌표 대응을 만들어야 한다.
#   검산: 2240 + ΣI − ΣD == unit_len 이 57,228 UNIT 중 57,151 (99.87%) 에서 성립 (실측).
#
# MLCUT = 128 (v11 노드 코드와 같게). beta_call 은 v11 의 beta 와 일치해야 한다 — 대조한다.
# 구획 (정본 2,240) : AluY 1-318(25) · NTS1 319-1178(78) · GA 1179-1551(29)
#                      NTS2 1552-2121(58, 끝 90bp 2032-2121 = 프로모터 12) · 유전자 2122-2240(9)
# 극성 : array_pos 1 = 게놈 원위. 이름(근위/원위)은 옛 판과 뒤집혔을 수 있으므로
#   여기서는 ap_dj(=array_pos) 와 ap_pj(=N-array_pos+1) 두 축을 다 저장하고 이름을 붙이지 않는다.
import io, os, re, gzip, sys, glob, time, collections
import numpy as np
import pandas as pd

V = __import__("os").environ["RDNA_DATA_ROOT"] + r"\for_paper\08_versions\v11"
OUT = __import__("os").environ["RDNA_DATA_ROOT"] + r"\v11_readvar"
CANON = __import__("os").environ["CANON_FA"]
MAP = os.path.join(V, r"03_maps\g8canonv10\CONSENSUS_MAP_G8_full.tsv")
RAW = os.path.join(V, r"08_methylation\raw")
SH = os.path.join(OUT, "shard")
os.makedirs(SH, exist_ok=True)
LOG = io.open(os.path.join(OUT, r"04_logs\V1_BUILD.txt"), "w", encoding="utf8", newline="\n")
MLCUT = 128


def P(s=""):
    try:
        print(s)
    except UnicodeEncodeError:
        print(s.encode("ascii", "replace").decode())
    LOG.write(s + "\n")
    LOG.flush()


canon = "".join(x.strip() for x in io.open(CANON, encoding="utf8") if x[0] != ">")
CPG = sorted(i + 1 for i in range(len(canon) - 1) if canon[i:i + 2] == "CG")
CPGSET = set(CPG)
CPGORD = {p: i for i, p in enumerate(CPG)}          # 정본 CpG 순번 (홀짝 반분할용)
SEG = [("alu", 1, 318), ("nts1", 319, 1178), ("ga", 1179, 1551),
       ("nts2", 1552, 2121), ("gene", 2122, 2240)]
UPS = (2032, 2121)                                   # 프로모터 (NTS2 끝 90bp)
SEGOF = {}
for nm, a, b in SEG:
    for p in CPG:
        if a <= p <= b:
            SEGOF[p] = nm
UPSET = set(p for p in CPG if UPS[0] <= p <= UPS[1])

P("=" * 104)
P("v11 read x 유닛 메틸 표 만들기")
P("=" * 104)
P()
P("  정본 CpG %d · 구획 %s · 프로모터 %d"
  % (len(CPG), " ".join("%s%d" % (n, sum(1 for p in CPG if a <= p <= b)) for n, a, b in SEG),
     len(UPSET)))
P("  calls 좌표는 슬롯 자신의 유닛 프레임 -> varlist 로 정본좌표에 옮긴다")
P()

# ---------------------------------------------------------------- 지도 읽기
IN = re.compile(r"^(\d+)I(\d+)$")
DE = re.compile(r"^(\d+)D(\d+)$")
UNITN = {}          # (sample, slot) -> N
UMAP = {}           # (sample, slot, ap) -> (슬롯좌표 -> 정본좌표) 배열
t0 = time.time()
nu = nbad = 0
with io.open(MAP, encoding="utf8") as h:
    hd = h.readline().rstrip("\n").split("\t")
    ix = {c: i for i, c in enumerate(hd)}
    for line in h:
        p = line.rstrip("\n").split("\t")
        if p[ix["seg_type"]] != "UNIT":
            continue
        smp = p[ix["sample"]]; sl = p[ix["slot"]]
        ap = int(p[ix["array_pos"]]); ul = int(p[ix["unit_len"]])
        UNITN[(smp, sl)] = int(p[ix["N"]])
        vl = p[ix["varlist"]]
        ins = {}; dele = {}
        if vl and vl != "-":
            for v in vl.split(","):
                m = IN.match(v)
                if m:
                    ins[int(m.group(1))] = int(m.group(2)); continue
                m = DE.match(v)
                if m:
                    dele[int(m.group(1))] = int(m.group(2))
        # 정본 1..2240 을 훑으며 슬롯좌표를 부여 -> 슬롯좌표별 정본좌표 배열
        arr = np.zeros(ul + 2, dtype=np.int16)
        sp = 1; cp = 1; skip = 0
        while cp <= 2240 and sp <= ul:
            if skip > 0:
                skip -= 1; cp += 1; continue
            if cp in dele:
                skip = dele[cp] - 1; cp += 1; continue
            arr[sp] = cp
            sp += 1
            if cp in ins:
                for _ in range(ins[cp]):
                    if sp <= ul:
                        arr[sp] = 0        # 삽입 염기는 정본에 대응 없음
                        sp += 1
            cp += 1
        if sp - 1 != ul:
            nbad += 1
        UMAP[(smp, sl, ap)] = arr
        nu += 1
P("  지도 UNIT %s · 길이 불일치 %d (%.4f) · %.0f 초"
  % (format(nu, ","), nbad, nbad / max(nu, 1), time.time() - t0))
P()

# ---------------------------------------------------------------- read-유닛
files = sorted(glob.glob(os.path.join(RAW, "METH_READUNIT_*.tsv.gz")))
P("  원천 파일 %d" % len(files))
P()
COLS = ["sample", "data_source", "cfg", "kind", "slot", "ap", "read_id", "strand",
        "N", "n_raw", "n_canon", "frac_canon", "b", "b_prob", "nA", "nB", "bA", "bB",
        "n_alu", "b_alu", "n_nts1", "b_nts1", "n_ga", "b_ga", "n_nts2", "b_nts2",
        "n_gene", "b_gene", "n_ups", "b_ups", "n_sp", "b_sp"]

for f in files:
    tag = os.path.basename(f)[len("METH_READUNIT_"):-len(".tsv.gz")]
    sp_out = os.path.join(SH, tag + ".parquet")
    if os.path.exists(sp_out):
        P("  건너뜀 (이미 있음) %s" % tag)
        continue
    t1 = time.time()
    rows = []
    nline = nmap0 = 0
    with gzip.open(f, "rt", encoding="utf8") as h:
        h.readline()
        for line in h:
            p = line.rstrip("\n").split("\t")
            smp, src, cfg, kind, sl = p[0], p[1], p[2], p[3], p[4]
            ap = int(p[5]); rid = p[6]; st = p[7]
            cs = p[11]
            if not cs or cs == ".":
                continue
            nline += 1
            arr = UMAP.get((smp, sl, ap))
            if arr is None:
                nmap0 += 1
                continue
            L = len(arr) - 1
            sumv = 0; nraw = 0
            tot = 0; met = 0
            segn = collections.Counter(); segm = collections.Counter()
            aN = aM = bN = bM = 0
            for it in cs.split(","):
                a, v = it.split(":")
                a = int(a); v = int(v)
                nraw += 1
                if a < 1 or a > L:
                    continue
                cp = int(arr[a])
                if cp == 0 or cp not in CPGSET:
                    continue
                m = 1 if v >= MLCUT else 0
                tot += 1; met += m; sumv += v
                sg = SEGOF[cp]
                segn[sg] += 1; segm[sg] += m
                if cp in UPSET:
                    segn["ups"] += 1; segm["ups"] += m
                if sg != "gene":
                    segn["sp"] += 1; segm["sp"] += m
                if CPGORD[cp] % 2 == 0:
                    aN += 1; aM += m
                else:
                    bN += 1; bM += m
            if tot == 0:
                continue
            r = [smp, src, cfg, kind, sl, ap, rid, st, UNITN.get((smp, sl), 0),
                 nraw, tot, tot / max(nraw, 1), met / tot, sumv / (255.0 * tot),
                 aN, bN, (aM / aN) if aN else np.nan, (bM / bN) if bN else np.nan]
            for k in ("alu", "nts1", "ga", "nts2", "gene", "ups", "sp"):
                n_ = segn[k]
                r += [n_, (segm[k] / n_) if n_ else np.nan]
            rows.append(r)
    D = pd.DataFrame(rows, columns=COLS)
    for c in ("ap", "N", "n_raw", "n_canon", "nA", "nB",
              "n_alu", "n_nts1", "n_ga", "n_nts2", "n_gene", "n_ups", "n_sp"):
        D[c] = D[c].astype("int32")
    for c in ("frac_canon", "b", "b_prob", "bA", "bB", "b_alu", "b_nts1", "b_ga",
              "b_nts2", "b_gene", "b_ups", "b_sp"):
        D[c] = D[c].astype("float32")
    for c in ("sample", "data_source", "cfg", "kind", "slot", "read_id", "strand"):
        D[c] = D[c].astype("category") if c != "read_id" else D[c].astype("string")
    D.to_parquet(sp_out, index=False)
    P("  %-28s 행 %9s · 지도없음 %6d · %.0f 초"
      % (tag, format(len(D), ","), nmap0, time.time() - t1))

P()
P("  샤드 합치기")
A = pd.concat([pd.read_parquet(x) for x in sorted(glob.glob(os.path.join(SH, "*.parquet")))],
              ignore_index=True)
P("  합계 %s 행" % format(len(A), ","))
P()

# ---------------------------------------------------------------- DUAL 정리
P("=" * 104)
P("[QC] 관문 점검")
P("=" * 104)
P()
P("  정본 CpG 적중 비율 (read-유닛당)")
P("     중앙 %.4f · 평균 %.4f · >=0.9 인 비율 %.4f"
  % (A.frac_canon.median(), A.frac_canon.mean(), (A.frac_canon >= 0.9).mean()))
P("     (옛 판 관문: >=0.9 가 98%% 이상이어야 정본 자료)")
P()
P("  유닛당 정본 CpG 개수: 중앙 %.0f · 사분위 %.0f-%.0f (최대 가능 199)"
  % (A.n_canon.median(), A.n_canon.quantile(.25), A.n_canon.quantile(.75)))
P()

dup = A.duplicated(subset=["sample", "slot", "ap", "read_id"], keep=False)
P("  같은 (개체,슬롯,자리,read) 가 두 번 (DUAL): %s 행" % format(int(dup.sum()), ","))
A["_pri"] = np.where(A.kind.astype(str) == "junction", 0, 1)
A = A.sort_values(["sample", "slot", "ap", "read_id", "_pri"])
A = A.drop_duplicates(subset=["sample", "slot", "ap", "read_id"], keep="first")
A = A.drop(columns=["_pri"])
P(" %s rows after junction-first consolidation" % format(len(A), ","))
P()

A["ap_dj"] = A.ap.astype("int32")
A["ap_pj"] = (A.N - A.ap + 1).astype("int32")
A.loc[A.N <= 0, "ap_pj"] = -1
A.to_parquet(os.path.join(OUT, r"02_tables\ru_v11.parquet"), index=False)
P("  저장: 02_표\\ru_v11.parquet")
P()
P("=" * 104)
P("[요약] 옛 판과 견주기")
P("=" * 104)
P()
H = A[A.data_source.astype(str) == "HPRC_ONT_submissions"]
P("   구분                 옛 (2026-09-09 · stageb)     새 (v11)")
P("   원천                 stageb\\tsv (정본 v7 이전)     METH_READUNIT (정본 v10 · 지도 g8canonv10)")
P("   read-유닛            555,895                      %s (HPRC %s)"
  % (format(len(A), ","), format(len(H), ",")))
P("   개체                 264                          %d (HPRC %d)"
  % (A["sample"].nunique(), H["sample"].nunique()))
P("   하플로타입/슬롯       528                          %d (HPRC %d)"
  % (A.groupby(["sample", "slot"], observed=True).ngroups,
     H.groupby(["sample", "slot"], observed=True).ngroups))
P("   read                 22,109                       %s (HPRC %s)"
  % (format(A.read_id.nunique(), ","), format(H.read_id.nunique(), ",")))
P("   array 안쪽 read      없음 (flank 만)                있음 (kind=array_only)")
P()
P("   출처별")
for k, g in A.groupby("data_source", observed=True):
    P("     %-24s 행 %9s · 개체 %3d · 슬롯 %3d · read %6s"
      % (k, format(len(g), ","), g["sample"].nunique(),
         g.groupby(["sample", "slot"], observed=True).ngroups,
         format(g.read_id.nunique(), ",")))
P()
P("   kind 별")
for k, g in A.groupby("kind", observed=True):
    P("     %-12s 행 %9s · read %6s · 자리 중앙 %.0f"
      % (k, format(len(g), ","), format(g.read_id.nunique(), ","), g.ap.median()))
P()
P("   N (array 길이) : 중앙 %.0f · 범위 %.0f-%.0f"
  % (A[A.N > 0].N.median(), A[A.N > 0].N.min(), A[A.N > 0].N.max()))
P()
LOG.close()
print("done")
