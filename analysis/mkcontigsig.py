#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 표지 contig 가 유의한가(진짜 자리 유래) 아니면 우연히 만들어졌나 — contig 를 이루는 read 자체의 증거로만 판정
#   다른 자료(PacBio↔ONT)의 '없음' 은 반증으로 쓰지 않는다 — 두 자료는 다른 연구 그룹 · 다른 세포 시료다
#   contig = mkmarkerasm v5 의 무리(read 2 개 이상 · 표지 유닛이 같은 칸) · read 분류 = mkreadcls v4
#   증거 (contig 마다)
#     ① 독립 read 수        : 서로 다른 read(분자) 수
#     ② 변이 염기 품질       : 각 read 에서 표지 대립형으로 읽힌 염기의 Phred Q (원 FASTQ) · Q>=20 인 read 수
#                              + 같은 read 의 같은 정규 자리(다른 유닛) Q 중앙과의 차 — 그 read 의 평소 품질 대비 국소 저품질인가
#     ③ read 끼리 맞음       : contig 구성 창끼리 표지 칸 고정 국소 정렬 점수(닮음 2J-1 · 틈 -0.8) · 쌍 중앙 · 음수 쌍 비율
#   판정 : 유의 = ①>=2 & ②(Q>=20 read)>=2 & ③ 쌍 중앙>=3(SMIN) · 아니면 유의 아님 + 이유
#   용법 : python3 mkcontigsig.py <CELL> <TAG> <fa>:<PB|ONT>:<자료> ... -- <자료>=<fastq.gz> ...
import sys, os, re, gzip, json, collections, statistics as st, random
argv = sys.argv; cut = argv.index("--"); FQ = dict(x.split("=", 1) for x in argv[cut + 1:]); sys.argv = argv[:cut]
exec(open(os.path.expanduser("~/p11_encode/mkcheck/mkmarkerasm.py")).read().split("RES = []; DIRECT = []")[0])
random.seed(20261014)
SPEC = {x.split(":")[2]: x.split(":")[0] for x in sys.argv[3:]}
# ---- paf 색인 : (read, qs, qe) → (가닥, ts, cigar)
PIDX = {}
for src, fa in SPEC.items():
    for L in open(fa + ".paf"):
        f = L.rstrip("\n").split("\t")
        cg = [x for x in f[12:] if x.startswith("cg:Z:")]
        if cg: PIDX[(f[0], int(f[2]), int(f[3]))] = (f[4], int(f[7]), cg[0][5:])
# ---- FASTQ 품질 (사슬 있는 read 만)
need = collections.defaultdict(set)
for q, r in READS.items(): need[r["src"]].add(q.split("|", 1)[1])
QUAL = {}
for src, path in FQ.items():
    with gzip.open(os.path.expanduser(path), "rt") as f:
        while True:
            h = f.readline()
            if not h: break
            f.readline(); f.readline(); qq = f.readline().rstrip("\n")
            nm = h[1:].split()[0]
            if nm in need[src]: QUAL[(src, nm)] = qq
log("품질 읽음", {s: sum(1 for k in QUAL if k[0] == s) for s in FQ})
def base_q(q, r, u, tpos):
    """read q 의 유닛 u 에서 정규 좌표 tpos 염기의 Q (없으면 None)"""
    rec = PIDX.get((q, u["qs"], u["qe"]))
    qq = QUAL.get((r["src"], q.split("|", 1)[1]))
    if not rec or qq is None: return None
    stn, ts, cg = rec; rp, tp = 0, ts
    for n, op in re.findall(r"(\d+)([MIDNSH=X])", cg):
        n = int(n)
        if op in "M=X":
            if tp < tpos <= tp + n:
                off = rp + (tpos - tp - 1); ro = u["qs"] + off if stn == "+" else u["qe"] - 1 - off
                return ord(qq[ro]) - 33 if 0 <= ro < len(qq) else None
            rp += n; tp += n
        elif op == "I": rp += n
        elif op in "DN":
            if tp < tpos <= tp + n: return None
            tp += n
    return None
# ---- contig · read 분류
G5 = json.load(open("%s/%s_markerasm_v5.json" % (W, TAG)))["results"]
RC = {(o["variant"], o["read"], o["k"]): o for o in json.load(open("%s/%s_readcls_v4.json" % (W, TAG)))}
VK = {m["variant"]: (p, a) for (p, a), m in MK.items()}
OUT = []
for gi, g in enumerate(G5):
    reads = {q for q, k in g["idx"]}
    if len(reads) < 2: continue
    v = g["variant"]; p, a = VK[v]; mz = g["mzone"]
    items = []
    for q, k in g["idx"]:
        r = READS[q]; u = r["ch"][k]
        qa = base_q(q, r, u, p)
        qref = [base_q(q, r, w, p) for i, w in enumerate(r["ch"]) if i != k and (p, a) not in w["raw"] and w["lo"] <= p <= w["hi"]]
        qref = [x for x in qref if x is not None]
        o = RC.get((v, q, k), {})
        lo, hi = max(0, k - HALF), min(len(r["ch"]), k + HALF + 1)
        items.append(dict(read=q, plat=r["plat"], k=k, q_alt=qa, q_ref_med=st.median(qref) if qref else None, cls=o.get("cls", ""), zone=o.get("zone", ""),
                          win=r["ch"][lo:hi], c=k - lo))
    # ③ 쌍 정렬 (최대 60 쌍)
    prs = [(i, j) for i in range(len(items)) for j in range(i + 1, len(items)) if items[i]["read"] != items[j]["read"]]
    random.shuffle(prs); sc = []
    for i, j in prs[:60]:
        A, Bw = items[i], items[j]
        sc.append(center_profile(A["win"], A["c"], Bw["win"], sim_f)[Bw["c"]])
    # 분류 (contig) — 구성 read 중 원래 자리 있으면 A
    cls = {it["cls"][:2].strip() for it in items}
    if "A" in cls: gcls = "A 원래 자리"
    else:
        zs = set()
        for it in items:
            if it["zone"] and it["zone"] != "미정":
                for z in it["zone"].split("|"): zs.add({"notD": "M", "notP": "M"}.get(z, z) if mz == "M" else z)
        zs -= {"D?", "P?"}
        gcls = "같은 구역 · 다른 자리/조립본 밖/모호" if (not zs or zs <= {mz} or mz in zs) else "다른 구역"
    # 독립 read (read 이름 · PacBio 는 ZMW 까지)
    mol = {re.sub(r"/ccs.*$", "", it["read"]) for it in items}
    qa = [it["q_alt"] for it in items if it["q_alt"] is not None]
    # 같은 read 가 여러 유닛을 가지면 read 당 최고 Q 하나
    best = {}
    for it in items:
        if it["q_alt"] is not None: best[it["read"]] = max(best.get(it["read"], -1), it["q_alt"])
    nq20 = sum(1 for x in best.values() if x >= 20)
    dq = [it["q_alt"] - it["q_ref_med"] for it in items if it["q_alt"] is not None and it["q_ref_med"] is not None]
    pmed = st.median(sc) if sc else None
    why = []
    if len(mol) < 2: why.append("독립 read <2")
    if nq20 < 2: why.append("변이 염기 Q>=20 read %d 개" % nq20)
    if pmed is None or pmed < SMIN: why.append("read 끼리 맞음 약함(쌍 중앙 %s)" % (None if pmed is None else round(pmed, 2)))
    OUT.append(dict(contig=gi, variant=v, mzone=mz, cls=gcls, n_read=len(reads), n_mol=len(mol), plat=dict(collections.Counter(READS[q]["plat"] for q in reads)),
                    q_alt=sorted(qa), n_q20=nq20, dq_med=st.median(dq) if dq else None, pair_med=None if pmed is None else round(pmed, 2),
                    pair_neg=round(sum(1 for x in sc if x < 0) / len(sc), 2) if sc else None, n_pair=len(sc),
                    sig="유의" if not why else "유의 아님", why="; ".join(why),
                    per_read=[dict(read=it["read"], plat=it["plat"], k=it["k"], q_alt=it["q_alt"], q_ref_med=it["q_ref_med"], cls=it["cls"], zone=it["zone"]) for it in items]))
json.dump(OUT, open("%s/%s_contigsig.json" % (W, TAG), "w"), ensure_ascii=False, indent=0)
c = collections.Counter((o["cls"], o["sig"]) for o in OUT)
log("contig", len(OUT), {"%s|%s" % k: n for k, n in sorted(c.items())})
