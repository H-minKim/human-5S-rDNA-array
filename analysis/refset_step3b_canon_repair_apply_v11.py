# -*- coding: utf-8 -*-
# 위 창을 정본에 끼워 cmp4d_v10 을 만든다
#   단계 : ③b 정본 수리 (v11 · 7 자리)
# Reason: no machine-specific path is left in the code. A default would leave such a path on this line, so none is set.
#       to run : export NODE_HOME=/home/<user>

# -*- coding: utf-8 -*-
# 정본 v10 — 수술식 교체 (2026-09-13)
#   v8(원조립본) 을 바탕으로, 끊긴 7 자리만 [조립본 앞앵커][read 합의 가운데][조립본 뒤앵커] 로 잇는다.
#   되돌림 없음 : v9 의 '앵커까지 합의로 교체' 가 문제였지 수리가 문제가 아니었다(SPLICE_MID2.tsv).
#   저장 방향 : cmp4d 의 .fa 는 뒤집기 전 서열. 지도 방향에서 만든 창은 되돌려(rc) 적는다.
import io
import os
import glob
import collections

H = os.environ["NODE_HOME"]   # the compute node's home directory
NL = chr(10)
TAB = chr(9)
CMP = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
SITES = [("HG01081", "2", 18), ("HG01106", "1", 37), ("HG02559", "2", 12),
         ("HG02624", "1", 53), ("NA20282", "2", 10), ("NA20850", "1", 161),
         ("HG00350", "2", 108)]


def rc(s):
    return "".join(CMP.get(c, "N") for c in reversed(s))


def rdfa(p):
    return "".join(x.strip() for x in io.open(p, encoding="utf8").read().split(NL)[1:]).upper()


SEQF = {}
for L in io.open(os.path.join(H, "hap_seq_node_v8.tsv"), encoding="utf8"):
    f = L.rstrip(NL).split(TAB)
    if len(f) >= 3:
        SEQF[(f[0], f[1])] = f[2]
FLIP = {}
f = io.open(os.path.join(H, "g8canonv8", "CONSENSUS_SUM_G8_full.tsv"), encoding="utf8")
f.readline()
for L in f:
    r = L.rstrip(NL).split(TAB)
    FLIP[(r[0], r[1])] = (r[4] == "Y")
U = collections.defaultdict(list)
f = io.open(os.path.join(H, "g8canonv8", "CONSENSUS_MAP_G8_full.tsv"), encoding="utf8")
ix = dict((n, i) for i, n in enumerate(f.readline().rstrip(NL).split(TAB)))
for L in f:
    r = L.rstrip(NL).split(TAB)
    if r[ix["seg_type"]] == "UNIT":
        U[(r[0], r[1])].append((int(r[ix["array_pos"]]), int(r[ix["unit_start"]]), int(r[ix["unit_end"]])))
for k in U:
    U[k].sort()

V10 = os.path.join(H, "cmp4d_v10")
os.makedirs(V10, exist_ok=True)
n = 0
for p in glob.glob(os.path.join(H, "cmp4d_v8", "*.fa")):
    b = os.path.basename(p)
    t = os.path.join(V10, b)
    if os.path.islink(t) or os.path.exists(t):
        os.remove(t)
    os.symlink(os.path.realpath(p), t)
    n += 1
print("v8 에서 이어받은 파일 %d" % n)

LOG = io.open(os.path.join(H, "CANON_V10_PATCH.tsv"), "w", encoding="utf8", newline=NL)
LOG.write(TAB.join(["sample", "slot", "after_unit", "file", "anchor_start", "anchor_end",
                    "asm_win_bp", "splice_win_bp", "hap_bp_old", "hap_bp_new", "flip"]) + NL)
for (smp, sl, ap) in SITES:
    fn = SEQF[(smp, sl)]
    raw = rdfa(os.path.join(H, "cmp4d_v8", fn))
    A = rc(raw) if FLIP.get((smp, sl)) else raw
    v = U[(smp, sl)]
    i = [j for j, x in enumerate(v) if x[0] == ap][0]
    a0, a1 = v[i][1], v[i + 1][2]
    spl = rdfa(os.path.join(H, "splice2", "spl_%s_%s_%d.fa" % (smp, sl, ap)))
    # 앵커 서열이 조립본과 글자 그대로 같은지 확인 (수술식의 전제)
    assert spl.startswith(A[v[i][1] - 1:v[i][2]]), "앞 앵커 불일치 %s %s" % (smp, sl)
    assert spl.endswith(A[v[i + 1][1] - 1:v[i + 1][2]]), "뒤 앵커 불일치 %s %s" % (smp, sl)
    B = A[:a0 - 1] + spl + A[a1:]
    out = rc(B) if FLIP.get((smp, sl)) else B
    t = os.path.join(V10, fn)
    if os.path.islink(t):
        os.remove(t)
    io.open(t, "w", newline=NL).write(">%s_splice_v10%s%s%s" % (fn.replace(".fa", ""), NL, out, NL))
    print("  %s h%s 유닛%d/%d : 창 %d -> %d bp · 전체 %d -> %d bp"
          % (smp, sl, v[i][0], v[i + 1][0], a1 - a0 + 1, len(spl), len(A), len(B)))
    LOG.write(TAB.join([smp, sl, str(ap), fn, str(a0), str(a1), str(a1 - a0 + 1),
                        str(len(spl)), str(len(A)), str(len(B)),
                        "Y" if FLIP.get((smp, sl)) else "N"]) + NL)
LOG.close()
real = sum(1 for x in os.listdir(V10) if not os.path.islink(os.path.join(V10, x)))
print("cmp4d_v10 : 전체 %d · 실제 파일 %d" % (len(os.listdir(V10)), real))
