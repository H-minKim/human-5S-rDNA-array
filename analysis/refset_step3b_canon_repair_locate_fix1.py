# -*- coding: utf-8 -*-
# read 조각 위치 잡기 (유전자끝 → 유닛 사슬 → 변이지문 → 유닛 DP) · 6 유닛 블록 minimap2 정렬은 실패해서 쓰지 않았다
#   단계 : ③b 정본 수리 (수정1 · 2 자리)

# -*- coding: utf-8 -*-
# v11-수정1 1단계 (2판) : 합의용 read 에서 [앞 앵커 유닛 시작 ~ 뒤 앵커 유닛 끝] 조각 — 유닛 단위 DP 로 위치를 잡는다
#
#   1판(seg_anchor.py) 실패 : 6 유닛 블록을 minimap2 로 read 에 붙였더니 일차·이차 점수차 0.7% ·
#                            앞/뒤 블록 일차 정렬이 read 같은 자리에 붙음 → ONT 오류가 유닛 간 변이보다 커서 염기 정렬로는 탠덤 위치 불가.
# second version : the functions of uindel2.py, the same call that chose the consensus reads, are used as they are.
#         read 유전자 끝(blastn) → 유닛 사슬 → 유닛 변이지문(asm20 --cs) → 지도와 유닛 DP → read 유닛 ↔ 지도 array_pos
#         조각 = 지도 유닛 aL 에 짝지어진 read 유닛의 시작 ~ aR 에 짝지어진 read 유닛의 끝 (read 방향 = 유전자 다수 가닥)
#   바꾸지 않은 것 : uindel2 의 gene_hits · gene_ends · chain_from · call_vars · score · dp · load_map (ast 로 함수만 추출)
#   인자 : uindel2.py경로 개체 슬롯 앞앵커 뒤앵커 read표 fastq.gz 출력폴더
import ast, io, os, sys, gzip, random, socket, statistics, collections
src_path, smp, sl, aL, aR, rtab, fq, OUT = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]), sys.argv[6], sys.argv[7], sys.argv[8]
tree = ast.parse(io.open(src_path, encoding="utf8").read())
keep_nodes = []
for node in tree.body:
    if isinstance(node, (ast.Import, ast.ImportFrom, ast.FunctionDef)):
        keep_nodes.append(node)
    elif isinstance(node, ast.Assign):
        names = [t.id for t in node.targets if isinstance(t, ast.Name)] + [e.id for t in node.targets if isinstance(t, ast.Tuple) for e in t.elts if isinstance(e, ast.Name)]
        if any(n in ("TGT", "FQ", "NP", "ONLY", "OUTD") for n in names):
            continue
        keep_nodes.append(node)
mod = ast.Module(body=keep_nodes, type_ignores=[])
G = {"__name__": "uindel2_funcs"}
exec(compile(mod, src_path, "exec"), G)
NL, TAB, H = G["NL"], G["TAB"], G["H"]
# uindel2 가 모듈 수준에서 채우는 표 두 개를 같은 규칙으로 채운다 (FN · SUMF)
for L in io.open(os.path.join(H, "hap_seq_node_v8.tsv"), encoding="utf8"):
    r = L.rstrip(NL).split(TAB)
    if len(r) >= 3 and r[0] != "sample": G["FN"][(r[0], r[1])] = r[2]
f = io.open(os.path.join(H, "g8canonv10", "CONSENSUS_SUM_G8_full.tsv"), encoding="utf8"); f.readline()
for L in f:
    r = L.rstrip(NL).split(TAB); G["SUMF"][(r[0], r[1])] = dict(N=int(r[3]), flip=(r[4] == "Y"), op=int(r[5]), od=int(r[6]))
os.makedirs(OUT, exist_ok=True)
W = os.path.join(OUT, "tmp_%s_%s" % (smp, sl)); os.makedirs(W, exist_ok=True)
rc, rdfa = G["rc"], G["rdfa"]
# 합의용 read
use = {}
f = io.open(rtab, encoding="utf8"); h = f.readline().rstrip(NL).split(TAB)
ci = [i for i, n in enumerate(h) if n.startswith("use_for_consensus")][0]
for L in f:
    r = L.rstrip(NL).split(TAB)
    if len(r) > ci and r[ci] == "Y": use[r[h.index("read16")]] = r[h.index("status")]
SEQ = {}
with gzip.open(fq, "rt") as fh:
    while True:
        hd = fh.readline()
        if not hd: break
        s = fh.readline().strip(); fh.readline(); fh.readline()
        rid = hd[1:].split()[0]
        if rid[:16] in use and rid[:16] not in SEQ: SEQ[rid[:16]] = (rid, s.upper())
seqs = dict((k, v[1]) for k, v in SEQ.items())
HT = G["gene_hits"](seqs, W, "a")
flip = {k for k, hs in HT.items() if sum(1 for x in hs if x[0] == "-") > sum(1 for x in hs if x[0] == "+")}
if flip:
    HT2 = G["gene_hits"](dict((k, rc(seqs[k])) for k in flip), W, "b")
    for k in flip: seqs[k] = rc(seqs[k]); HT[k] = HT2.get(k, [])
MAPS, FIX = G["load_map"](smp)
mv, aps = MAPS[sl]
CH = {k: G["chain_from"](G["gene_ends"](HT.get(k, []), -10 ** 9, 10 ** 9)) for k in seqs}
useqs, uidx = [], []
for k, ch in CH.items():
    for i, it in enumerate(ch):
        if it[0] == "U": useqs.append(seqs[k][max(0, it[1] - 1 - G["PAD"]):it[2] + G["PAD"]]); uidx.append((k, i))
VV = G["call_vars"](useqs, W, "u") if useqs else {}
FP = collections.defaultdict(dict)
for n_, (k, i) in enumerate(uidx): FP[k][i] = VV.get(n_, set()) - FIX
random.seed(20260916)
pool_r = [v for d in FP.values() for v in d.values() if v]; pool_m = [v for v in mv if isinstance(v, set) and v]
base = statistics.median([G["jac"](random.choice(pool_r), random.choice(pool_m)) for _ in range(2000)]) if pool_r and pool_m else 0.1
# 사이 유닛 : 지도에서 aL..aR 사이 정본 창 길이
U = {}
f = io.open(os.path.join(H, "g8canonv10", "CONSENSUS_MAP_G8_full.tsv"), encoding="utf8"); ix = dict((n, i) for i, n in enumerate(f.readline().rstrip(NL).split(TAB)))
for L in f:
    r = L.rstrip(NL).split(TAB)
    if (r[0], r[1]) == (smp, sl) and r[ix["seg_type"]] == "UNIT": U[int(r[ix["array_pos"]])] = (int(r[ix["unit_start"]]), int(r[ix["unit_end"]]))
asm_win = U[aR][1] - U[aL][0] + 1
host = socket.gethostname()
T = collections.Counter(); LOG = []
outf = os.path.join(OUT, "seg_%s_%s_%d.%s.fa" % (smp, sl, aL, host))
with io.open(outf, "w", newline=NL) as out:
    for k in sorted(seqs):
        T["read"] += 1
        ch = CH[k]
        # 끊김(B) 으로 나뉜 조각별 DP — uindel2 와 같은 분할
        segs, cur = [], []
        for i, it in enumerate(ch):
            if it[0] == "B":
                if cur: segs.append(cur)
                cur = []
            else:
                cur.append(i)
        if cur: segs.append(cur)
        hitL = hitR = None; nseg = 0; evs = []
        for sg in segs:
            R = [FP[k].get(i, set()) if ch[i][0] == "U" else set() for i in sg]
            if sum(1 for i in sg if ch[i][0] == "U") < G["MINSEG"]: continue
            nseg += 1
            sc, ops = G["dp"](R, mv, base)
            evs += [e for e in G["events"](ops, aps, R)]
            for o in ops:
                if o[0] != "M": continue
                ci_ = sg[o[1]]; ap = aps[o[2]]
                if ch[ci_][0] != "U": continue
                if ap == aL: hitL = ch[ci_]
                if ap == aR: hitR = ch[ci_]
        if hitL is None or hitR is None:
            T["앵커 유닛이 DP 짝에 없음"] += 1
            LOG.append((k, "앵커 짝 없음 L=%s R=%s · 사건 %s" % (hitL is not None, hitR is not None, evs))); continue
        s0, e1 = hitL[1], hitR[2]                     # 1 기반 : 앞 앵커 유닛 시작 · 뒤 앵커 유닛 끝
        if e1 - s0 < 500:
            T["조각 짧음/순서"] += 1; LOG.append((k, "짧음")); continue
        seg = seqs[k][s0 - 1:e1]
        out.write(">%s|%s|0|0%s%s%s" % (SEQ[k][0], use.get(k, "?"), NL, seg, NL))
        T["조각"] += 1
        LOG.append((k, "조각 %d bp · 조립본 창 %d 대비 %+.2f 유닛 · 조각 안 read 유닛 %d · 사건 %s" % (
            len(seg), asm_win, (len(seg) - asm_win) / 2240.0,
            sum(1 for it in ch if it[0] == "U" and it[1] >= s0 and it[2] <= e1),
            ",".join("%s%d@%s" % (e[0], e[1], e[2]) for e in evs) or "-")))
io.open(os.path.join(OUT, "SEG_LOG_%s_%s.tsv" % (smp, sl)), "w", encoding="utf8", newline=NL).write(
    "read16\tresult\n" + "".join("%s\t%s\n" % r for r in LOG))
for fn in os.listdir(W):
    try: os.remove(os.path.join(W, fn))
    except OSError: pass
os.rmdir(W)
print("%s h%s 앵커 %d/%d · 조립본 창 %d bp · 합의용 %d · fastq 찾음 %d · 지문 기준값 %.3f" % (smp, sl, aL, aR, asm_win, len(use), len(seqs), base))
for kk in sorted(T): print("   %-26s %d" % (kk, T[kk]))
for r in LOG: print("   %s  %s" % r)
