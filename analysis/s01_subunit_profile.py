# -*- coding: utf-8 -*-
"""S01 — 유닛 안 위치별(서브유닛) A 자리당 m6A 프로파일을 저메틸/고메틸 따로 만든다.

왜
  Figure 1 의 구획별 막대는 '평균'이라 구획 안 모양을 못 보여 준다. 원자료 calls(unit_pos:ML)를 그대로 세서
  유닛 좌표(실제 bp)에서의 윤곽을 L·H 따로 낸다.

입력 (읽기만)
  <RDNA_V12_CANON_DIR>/08_methylation_s6/<노드>/s6/{ont_meth,pb_meth}/M6A_READUNIT_*.tsv.gz
      calls 열 = "unit_pos:ML,..." (unit_pos = 유닛 안 1기반 · 두 가닥 A 자리 모두) · in_final=Y 만
  유닛 좌표·TSS : a_lib.load_units()  (MAP_V11 + 정본 서열 · TSS = 유닛 안 마지막 GTCTACGGC)
  유닛 상태 : s6 판 a13_unit_state_gmm.tsv.gz (L P>=0.9 · H P<=0.1 · U 제외)

방법
  tx = unit_pos - tss_u  (tss_u = 유닛 안 TSS 1기반 위치) → 유전자 tx 0..118 · AluY 는 tx 약 -2120..-1800
  BIN bp 칸마다 개체별 sum(m6a)/sum(A) → 개체 평균 ± 95% CI (개체 안에서 A>=30 인 칸만)
출력
  s01_subunit_profile.tsv     (src·state·tx_bin·n_sample·mean·ci95·n_a)
  s01_element_bounds.tsv      (구획 경계 중앙값 · 그림 x축 눈금용)
사용 : python s01_subunit_profile.py [BIN]      기본 BIN=20
"""
import io, os, sys, gzip, collections
import numpy as np

# The data and output roots are environment variables with no default (CONVENTIONS, "the deposited copy wins").
#  저장소 안 파일로 가는 길은 `__file__` 에서 상대로 푼다 — `a_lib`·`bamlite`·`p02_lib` 가 이 폴더에 같이 예치돼 있다.
#  변수 표는 `README.md` §2 에 있다.
_HERE = os.path.dirname(os.path.abspath(__file__))


def _need(k):
    try:
        return os.environ[k]
    except KeyError:
        raise SystemExit(
            "environment variable %s is not set; it is required by %s. "
            "The deposited copy has no default; see README.md section 2."
            % (k, os.path.basename(__file__)))

sys.path.insert(0, _HERE)
import a_lib as L

BIN = int(sys.argv[1]) if len(sys.argv) > 1 else 20
#  쓰는 폴더다 — 발행 입력 폴더를 가리키면 그 표를 덮는다.
OUT = _need("RDNA_RESULTMAP_DIR")
S6 = os.path.join(L.R, "08_methylation_s6")
MDIR = {"UCSC_ONT_FIBERSEQ": "ont_meth", "PACBIO_FIBERSEQ": "pb_meth"}
# PacBio 코호트 = jasmine 48 : ft0.4.2 2 명(HG02630·hg002) + m6A 없는 jasmine2.2.0 2 명(HG02624·HG03672) 제외
EXCL = {"PACBIO_FIBERSEQ": {"HG02630", "hg002", "HG02624", "HG03672"}, "UCSC_ONT_FIBERSEQ": set()}
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace", line_buffering=True)

print("정본 유닛·TSS 적재")
UNIT, _SEQ = L.load_units()
TSSU = {}                      # (sample, slot, apos) -> 유닛 안 TSS 1기반 · 유닛 길이
for k, (us, ue, tss) in UNIT.items():
    TSSU[k] = (tss - us + 1, ue - us + 1)
print("  유닛 %d 개 · 유닛 안 TSS 중앙 %d · 유닛 길이 중앙 %d"
      % (len(TSSU), int(np.median([v[0] for v in TSSU.values()])), int(np.median([v[1] for v in TSSU.values()]))))

print("s6 상태 적재")
ST = {}
with gzip.open(os.path.join(_need("RDNA_S6_DIR"), "a13_unit_state_gmm.tsv.gz"),
               "rt", encoding="utf8") as f:
    f.readline()
    for l in f:
        x = l.rstrip("\n").split("\t")
        ST[(x[1], x[0], x[2], x[3], x[4])] = x[10]     # (sample, src, read, slot, apos) -> L/H/U
print("  유닛 상태 %d" % len(ST))

acc = collections.defaultdict(lambda: [0, 0])          # (src, state, sample, bin) -> [A, m6A]
chk = {"row": 0, "no_unit": 0, "no_state": 0, "U": 0, "pos_gt_len": 0, "used": 0}
for src, sub in MDIR.items():
    for node in ("w3", "iB", "iA"):
        p = os.path.join(S6, node, "s6", sub, "M6A_READUNIT_%s.tsv.gz" % L.HOST[node])
        if not os.path.exists(p):
            print("  없음", p); continue
        with gzip.open(p, "rt", encoding="utf8") as f:
            ix = {k: i for i, k in enumerate(f.readline().rstrip("\n").split("\t"))}
            for l in f:
                x = l.rstrip("\n").split("\t")
                chk["row"] += 1
                if x[ix["data_source"]] != src or x[ix["in_final"]] != "Y":
                    continue
                smp = x[ix["sample"]]
                if smp in EXCL[src]:
                    continue
                key = (smp, x[ix["slot"]], x[ix["array_pos"]])
                if key not in TSSU:
                    chk["no_unit"] += 1; continue
                st = ST.get((smp, src, x[ix["read_id"]], x[ix["slot"]], x[ix["array_pos"]]))
                if st is None:
                    chk["no_state"] += 1; continue
                if st not in ("L", "H"):
                    chk["U"] += 1; continue
                tss_u, ulen = TSSU[key]
                calls = x[ix["calls"]]
                if not calls:
                    continue
                for c in calls.split(","):
                    q = c.split(":")
                    if len(q) != 2:
                        continue
                    pos = int(q[0])
                    if pos > ulen:
                        chk["pos_gt_len"] += 1
                    b = ((pos - tss_u) // BIN) * BIN
                    v = acc[(src, st, smp, b)]
                    v[0] += 1
                    if int(q[1]) >= 128:
                        v[1] += 1
                chk["used"] += 1
        print("  %s %s 끝 (누적 read×유닛 %d)" % (src, node, chk["used"]))
print("점검 :", chk)

# 개체 평균 ± 95% CI
per = collections.defaultdict(dict)
for (src, st, smp, b), (na, nm) in acc.items():
    if na >= 30:
        per[(src, st, b)][smp] = nm / float(na)
tot = collections.defaultdict(lambda: [0, 0])
for (src, st, smp, b), (na, nm) in acc.items():
    t = tot[(src, st, b)]
    t[0] += na; t[1] += nm

with io.open(os.path.join(OUT, "s01_subunit_profile.tsv"), "w", encoding="utf8", newline="\n") as g:
    g.write("src\tstate\ttx_bin\tbin_bp\tn_sample\tmean\tci95\tpooled_frac\tn_a\n")
    for (src, st, b) in sorted(per, key=lambda k: (k[0], k[1], k[2])):
        v = np.array(list(per[(src, st, b)].values()), dtype=float)
        if len(v) < 5:
            continue
        ci = 1.96 * v.std(ddof=1) / np.sqrt(len(v))
        na, nm = tot[(src, st, b)]
        g.write("%s\t%s\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%d\n"
                % (src, st, b, BIN, len(v), v.mean(), ci, nm / float(na), na))

# 구획 경계 (유닛 안 좌표 · tx 로 환산 · 중앙값)
tss_med = int(np.median([v[0] for v in TSSU.values()]))
bounds = [("AluY", 1, 318), ("NTS1", 319, 1178), ("GA", 1179, tss_med - 571),
          ("NTS2", tss_med - 570, tss_med - 91), ("상류90", tss_med - 90, tss_med - 1),
          ("유전자", tss_med, tss_med + 118)]
with io.open(os.path.join(OUT, "s01_element_bounds.tsv"), "w", encoding="utf8", newline="\n") as g:
    g.write("element\tunit_start\tunit_end\ttx_start\ttx_end\ttss_u_median\n")
    for e, a, b in bounds:
        g.write("%s\t%d\t%d\t%d\t%d\t%d\n" % (e, a, b, a - tss_med, b - tss_med, tss_med))
print("저장 : s01_subunit_profile.tsv · s01_element_bounds.tsv (BIN=%d · 유닛 안 TSS 중앙 %d)" % (BIN, tss_med))
