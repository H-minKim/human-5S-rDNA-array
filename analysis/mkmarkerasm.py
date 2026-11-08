#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Builds a new contig from the reads of each marker and decides whether the marker position falls in the original interval (read-only)
#   설계 :
#     ① 표지 대립형을 가진 read 를 전부 모은다 (조립본에 붙든 안 붙든)
#     ② 그 read 들끼리만 정렬해 contig 를 만든다 (조립본 안 씀) — 표지 유닛을 가운데로 한 창(±10 유닛)끼리, 표지 유닛끼리 같은 칸에 맞을 때만 잇는다
#     ③ contig 안에서 read 들의 표지 유닛이 같은 칸 = 가운데 칸 (잇는 조건 자체가 이것)
#     ④ contig 합의를 조립본과 비교
#        (가) 이웃 유닛으로 : 합의 창을 지도에 가운데 정렬 → 가운데가 놓인 지도 유닛·구역
#        (나) 배열 끝에서 잰 거리로 : contig 를 이루는 read 중 배열 끝(flank)까지 닿는 read 가 표지 유닛까지 센 유닛 수
#        (다) 둘 다 안 되면 '위치 미정'
#     ⑤ read 들이 표지 자리에 동의하는데 그 자리가 원래 구간이 아닌 contig 를 보고
#   추가 : contig 와 무관하게, 배열 끝까지 닿는 read 하나하나가 직접 잰 '끝에서 표지까지 거리'도 표지별로 낸다(조립본 없는 직접 측정)
#   방향 : 유닛 사슬은 정규 방향 차례 · 놓인 read 의 방향이 전부 +1 ⇒ 사슬 시작 = 원위 끝 · 사슬 끝 = 근위 끝 (실행 중 끝 닿는 놓인 read 로 다시 확인)
#   용법 : python3 mkmarkerasm.py <CELL> <TAG> <fa>:<PB|ONT>:<자료이름> [...]
import sys, os, re, json, math, collections, time, random
W = os.path.expanduser("~/p11_encode/mkcheck")
CELL, TAG = sys.argv[1], sys.argv[2]; SPECS = sys.argv[3:]
CAN = "".join(l.strip() for l in open(os.path.expanduser("~/s5s_tuc/ref/unit_canon.fa")) if not l.startswith(">")).upper()
EXCL = lambda p: p <= 318 or 1179 <= p <= 1551
S1R = re.compile(r"^(\d+)([ACGT])$")
SMIN, MMIN, GAP, HALF = 3.0, 1.0, -0.8, 10
FLANK_MIN = 3000                     # 마지막 5S 조각 너머 5S 아닌 서열이 이만큼 있어야 '끝에 닿음'
t0 = time.time()
def log(*a): print("[%5.0fs]" % (time.time() - t0), *a, flush=True)
random.seed(20261014)
# ---------------- 지도·표지·경계
FR = collections.defaultdict(list)
for i, line in enumerate(open("%s/mk_frames_%s.tsv" % (W, CELL))):
    if i == 0: continue
    s, N, ap, vl = line.rstrip("\n").split("\t")
    vs = {(int(m.group(1)), m.group(2)) for v in vl.split(",") for m in [S1R.match(v)] if m and not EXCL(int(m.group(1)))}
    FR[s].append((int(ap), vs))
for s in FR: FR[s].sort()
CNT = collections.Counter(v for s in FR for _, vs in FR[s] for v in vs); NU = sum(len(FR[s]) for s in FR); ALLF = set(CNT)
WMAX = math.log(NU + 1.0)
def wt(v): return math.log((NU + 1.0) / (CNT.get(v, 0) + 1.0)) if v in CNT else WMAX / 2
MK = {}
for i, line in enumerate(open("%s/mk_markers_%s.tsv" % (W, CELL))):
    if i == 0: continue
    v, p, a, z, units = line.rstrip("\n").split("\t"); car = collections.defaultdict(set)
    for t in units.split(";"):
        if t: sl, ap = t.split(":"); car[sl].add(int(ap))
    MK[(int(p), a)] = dict(variant=v, zone=z, car=car)
FRSET = collections.defaultdict(set)                     # 지도에 그 변이가 있는 유닛(HOT 포함)
for s in FR:
    for ap, vs in FR[s]:
        for v in vs: FRSET[v].add((s, ap))
B = json.load(open("%s/mk_bounds_%s.json" % (W, CELL)))
def zone(slot, ap):
    kp, kd, N = B["new"][slot]
    return "D" if ap <= kd else ("P" if N - ap + 1 <= kp else "M")
KD = {s: B["new"][s][1] for s in B["new"]}; KP = {s: B["new"][s][0] for s in B["new"]}
# ---------------- read → 유닛 사슬 (+ 끝 닿음)
COMP = str.maketrans("ACGTN", "TGCAN")
def rc(s): return s.translate(COMP)[::-1]
def load(fa, plat):
    seqs = {}; k = None
    for L in open(fa):
        if L.startswith(">"): k = L[1:].strip(); seqs[k] = []
        else: seqs[k].append(L.strip())
    seqs = {k: "".join(v) for k, v in seqs.items()}
    idmin = 0.90 if plat == "PB" else 0.80
    seg = collections.defaultdict(list); anyseg = collections.defaultdict(list)
    for L in open(fa + ".paf"):
        f = L.rstrip("\n").split("\t")
        q, qs, qe, st, ts, te, nm, bl = f[0], int(f[2]), int(f[3]), f[4], int(f[7]), int(f[8]), int(f[9]), int(f[10])
        if bl >= 300: anyseg[q].append((qs, qe))
        if bl < 1800 or nm / bl < idmin: continue
        cg = [x for x in f[12:] if x.startswith("cg:Z:")]
        if cg: seg[q].append((nm, qs, qe, st, ts, te, cg[0][5:]))
    out = {}
    for q, v in seg.items():
        keep = []
        for x in sorted(v, reverse=True):
            if all(min(x[2], kk[2]) - max(x[1], kk[1]) < 0.5 * (x[2] - x[1]) for kk in keep): keep.append(x)
        stm = collections.Counter(kk[3] for kk in keep).most_common(1)[0][0]
        keep = [kk for kk in keep if kk[3] == stm]; keep.sort(key=lambda kk: kk[1], reverse=(stm == "-"))
        if len(keep) < 3: continue
        ch = []
        for nm, qs, qe, st, ts, te, cg in keep:
            sub = seqs[q][qs:qe]; sub = rc(sub) if st == "-" else sub
            rp, tp, d = 0, ts, set()
            for n, op in re.findall(r"(\d+)([MIDNSH=X])", cg):
                n = int(n)
                if op in "M=X":
                    for i in range(n):
                        b = sub[rp + i] if rp + i < len(sub) else "N"
                        if b != CAN[tp + i] and b in "ACGT": d.add((tp + i + 1, b))
                    rp += n; tp += n
                elif op == "I": rp += n
                elif op in "DN": tp += n
                else: rp += n
            raw = {x for x in d if not EXCL(x[0])}
            ch.append(dict(lo=ts + 1, hi=te, raw=raw, f=({x for x in raw if x in ALLF} if plat == "ONT" else raw), qs=qs, qe=qe, plat=plat))
        # 5S 조각 슬롯(느슨 · >=1000bp) 으로 끝까지 센다 · 사슬 방향(정규) 기준 시작=원위 끝
        L = len(seqs[q]); iv = sorted(x for x in anyseg[q] if x[1] - x[0] >= 1000)
        slots = []
        for a, b in iv:
            if slots and a < slots[-1][1] - 0.5 * (b - a): slots[-1] = (slots[-1][0], max(slots[-1][1], b))
            else: slots.append((a, b))
        if stm == "-": slots = [(L - b, L - a) for a, b in slots][::-1]     # 정규 방향 read 좌표로
        lo5 = min(a for a, b in anyseg[q]); hi5 = max(b for a, b in anyseg[q])
        if stm == "-": lo5, hi5 = L - hi5, L - lo5
        sq = seqs[q] if stm == "+" else rc(seqs[q])          # 정규 방향 read 서열
        tails = {}
        if lo5 >= FLANK_MIN: tails["시작"] = sq[max(0, lo5 - 20000):lo5]
        if L - hi5 >= FLANK_MIN: tails["끝"] = sq[hi5:hi5 + 20000]
        start_flank = end_flank = False                       # flank 정렬로 나중에 정한다
        for u in ch:                                   # 유닛의 슬롯 번호(정규 방향 · 0 부터)
            c = ((L - u["qe"] + L - u["qs"]) / 2) if stm == "-" else (u["qs"] + u["qe"]) / 2
            u["slot_i"] = min(range(len(slots)), key=lambda i: abs((slots[i][0] + slots[i][1]) / 2 - c)) if slots else None
        out["%s" % q] = dict(ch=ch, nslot=len(slots), sf=start_flank, ef=end_flank, plat=plat, tails=tails)
    return out
# ---------------- 닮음 · 가운데 정렬 (mkcheck 와 같은 원리)
SFC = {}
def sim_f(a, b):                                       # read 유닛끼리 (둘 중 ONT 가 있으면 지도 어휘만) · 캐시
    key = (id(a), id(b))
    if key not in SFC: SFC[key] = _sim_f(a, b)
    return SFC[key]
def _sim_f(a, b):
    lo, hi = max(a["lo"], b["lo"]), min(a["hi"], b["hi"])
    if hi - lo < 1500: return 0.0
    voc = a["plat"] == "ONT" or b["plat"] == "ONT"
    A = {v for v in a["f"] if lo <= v[0] <= hi and (not voc or v in ALLF)}; Bs = {v for v in b["f"] if lo <= v[0] <= hi and (not voc or v in ALLF)}
    U = A | Bs
    if not U: return 0.0
    return 2.0 * sum(wt(v) for v in A & Bs) / sum(wt(v) for v in U) - 1.0
def sim_map(u, fs):                                    # 합의 유닛 ↔ 지도 유닛 (합의 특징 중 지도 어휘만 비교 · 새 변이는 따로 보고)
    Fs = {v for v in fs if u["lo"] <= v[0] <= u["hi"]}; R = {v for v in u["f"] if v in ALLF}; U = R | Fs
    if not U: return 0.0
    return 2.0 * sum(wt(v) for v in R & Fs) / sum(wt(v) for v in U) - 1.0
def center_profile(win, c, tgt, sfun):
    n, m = len(win), len(tgt)
    S = [[sfun(win[i], tgt[j]) for j in range(m)] for i in range(n)]
    Hf = [[0.0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            Hf[i][j] = max(0.0, Hf[i - 1][j - 1] + S[i - 1][j - 1], Hf[i - 1][j] + GAP, Hf[i][j - 1] + GAP)
    Hb = [[0.0] * (m + 2) for _ in range(n + 2)]
    for i in range(n, 0, -1):
        for j in range(m, 0, -1):
            Hb[i][j] = max(0.0, Hb[i + 1][j + 1] + S[i - 1][j - 1], Hb[i + 1][j] + GAP, Hb[i][j + 1] + GAP)
    ci = c + 1
    return [max(0.0, Hf[ci - 1][j - 1]) + S[c][j - 1] + max(0.0, Hb[ci + 1][j + 1]) for j in range(1, m + 1)]
def place_on_map(win, c):
    cand = []
    for s in sorted(FR):
        units = FR[s]
        for dirn in (1, -1):
            ii = list(range(len(units))) if dirn == 1 else list(range(len(units)))[::-1]
            prof = center_profile(win, c, [units[j][1] for j in ii], sim_map)
            for k, sc in enumerate(prof): cand.append((sc, s, units[ii[k]][0], dirn))
    cand.sort(reverse=True)
    best = cand[0]; second = next((x for x in cand[1:] if (x[1], x[2]) != (best[1], best[2])), (0.0, None, None, 0))
    LAST_AMB[:] = sorted({(x[1], x[2], zone(x[1], x[2])) for x in cand if x[0] >= best[0] - MMIN and x[0] > 0})
    return best, second
LAST_AMB = []
def sw_path(X, F):
    """X : 합의 칸들 · F : 지도 유닛 [(ap, fs)] → (점수, [(칸 i, 지도 j)], S)"""
    n, m = len(X), len(F)
    S = [[sim_map(X[i], F[j][1]) for j in range(m)] for i in range(n)]
    H = [[0.0] * (m + 1) for _ in range(n + 1)]; P = [[0] * (m + 1) for _ in range(n + 1)]; best = (0.0, 0, 0)
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            a_, b_, c_ = H[i - 1][j - 1] + S[i - 1][j - 1], H[i - 1][j] + GAP, H[i][j - 1] + GAP
            h = max(0.0, a_, b_, c_); H[i][j] = h; P[i][j] = 0 if h == 0 else (1 if h == a_ else (2 if h == b_ else 3))
            if h > best[0]: best = (h, i, j)
    sc, i, j = best; pairs = []
    while i > 0 and j > 0 and P[i][j]:
        if P[i][j] == 1: pairs.append((i - 1, j - 1)); i -= 1; j -= 1
        elif P[i][j] == 2: i -= 1
        else: j -= 1
    return sc, pairs[::-1], S
def informative(u, fs):
    return bool({vv for vv in u["f"] if vv in ALLF}) or bool({vv for vv in fs if u["lo"] <= vv[0] <= u["hi"]})
# ---------------- 실행
READS = {}
for spec in SPECS:
    fa, plat, src = spec.split(":")
    for q, r in load(fa, plat).items(): r["src"] = src; READS[q] = r
    log(src, "사슬", sum(1 for r in READS.values() if r["src"] == src))
# 끝 닿음 = 꼬리가 flank 서열에 정렬 (>=1000bp · 일치 >=0.8) · 사슬 쪽(시작/끝) ↔ flank 파일 대응을 센다
FLK = os.path.expanduser("~/p11_encode/ref/_flank_both.fa")
with open(FLK, "w") as o:
    for nm in ("FLANK_distal_50000", "FLANK_proximal_50000"):
        o.write(open(os.path.expanduser("~/p11_encode/ref/%s.fa" % nm)).read().rstrip("\n") + "\n")
chk = collections.Counter(); HIT = {}
for plat in ("PB", "ONT"):
    tf = "%s/%s_tails_%s.fa" % (W, TAG, plat); n = 0
    with open(tf, "w") as o:
        for q, r in READS.items():
            if r["plat"] != plat: continue
            for side, s in r["tails"].items(): o.write(">%s@@%s\n%s\n" % (q, side, s)); n += 1
    if not n: continue
    paf = tf + ".paf"
    os.system("minimap2 -x %s -t 12 %s %s > %s 2>/dev/null" % ("map-hifi" if plat == "PB" else "map-ont", FLK, tf, paf))
    for L in open(paf):
        f = L.split("\t"); bl, nm_ = int(f[10]), int(f[9])
        if bl >= 1000 and nm_ / bl >= 0.8:
            q, side = f[0].split("@@"); k = (q, side)
            if k not in HIT or bl > HIT[k][1]: HIT[k] = (f[5], bl)
for (q, side), (fn, bl) in HIT.items():
    r = READS[q]
    if side == "시작": r["sf"] = True
    else: r["ef"] = True
    chk[(side, fn)] += 1
log("끝 닿음(flank 정렬) 사슬쪽 × flank 파일", dict(chk), "· 꼬리 후보", sum(len(r["tails"]) for r in READS.values()))
RES = []; DIRECT = []
for (p, a), mk in sorted(MK.items(), key=lambda kv: kv[1]["variant"]):
    v = mk["variant"]; items = []
    for q, r in READS.items():
        ch = r["ch"]
        for k, u in enumerate(ch):
            if (p, a) in u["raw"]:
                lo, hi = max(0, k - HALF), min(len(ch), k + HALF + 1)
                items.append(dict(read=q, src=r["src"], plat=r["plat"], win=ch[lo:hi], c=k - lo, k=k, r=r))
                # 직접 측정 : 끝에서 표지 유닛까지(슬롯 번호로 · 1 = 끝 유닛)
                if u["slot_i"] is not None:
                    if r["sf"]: DIRECT.append(dict(variant=v, mzone=mk["zone"], src=r["src"], read=q, end="원위", dist=u["slot_i"] + 1))
                    if r["ef"]: DIRECT.append(dict(variant=v, mzone=mk["zone"], src=r["src"], read=q, end="근위", dist=r["nslot"] - u["slot_i"]))
    if not items: continue
    # ② 표지 유닛끼리 같은 칸에 맞는 짝 → 탐욕 합치기(합칠 때 새로 겹치는 짝 검사)
    n = len(items); E = []
    for i in range(n):
        for j in range(i + 1, n):
            A, Bw = items[i], items[j]
            if A["read"] == Bw["read"]: continue
            prof = center_profile(A["win"], A["c"], Bw["win"], sim_f)
            cc = prof[Bw["c"]]
            if cc >= SMIN and cc >= max(prof) - 0.5: E.append((cc, i, j))
    comp = {}; mem = {}; NEXT = [0]
    def ok_merge(X, Y):
        prs = [(x, y) for x in X for y in Y if items[x]["read"] != items[y]["read"]]; random.shuffle(prs); bad = good = 0
        for x, y in prs[:40]:
            A, Bw = items[x], items[y]
            sc = center_profile(A["win"], A["c"], Bw["win"], sim_f)[Bw["c"]]
            if sc > 0: good += 1
            else: bad += 1
        return not (bad >= 3 or (bad + good and bad / (bad + good) > 0.2))
    for sc, i, j in sorted(E, reverse=True):
        if i not in comp and j not in comp:
            cid = NEXT[0]; NEXT[0] += 1; mem[cid] = [i, j]; comp[i] = comp[j] = cid
        elif i in comp and j not in comp:
            if ok_merge([j], mem[comp[i]]): comp[j] = comp[i]; mem[comp[i]].append(j)
        elif j in comp and i not in comp:
            if ok_merge([i], mem[comp[j]]): comp[i] = comp[j]; mem[comp[j]].append(i)
        elif comp[i] != comp[j]:
            ci, cj = comp[i], comp[j]
            if ok_merge(mem[cj], mem[ci]):
                for x in mem[cj]: comp[x] = ci
                mem[ci] += mem[cj]; del mem[cj]
    groups = list(mem.values()) + [[i] for i in range(n) if i not in comp]
    for g in groups:
        reads = sorted({items[x]["read"] for x in g})
        # 합의 창 : 가운데 기준 칸 −10..+10 · 덮음>=2(홀로면 그 read 그대로) · 과반
        cols = collections.defaultdict(list)
        for x in g:
            it = items[x]
            for k, u in enumerate(it["win"]): cols[k - it["c"]].append(u)
        need = 1 if len(reads) == 1 else 2
        rel = sorted(c for c in cols if len(cols[c]) >= need)
        cons = []
        for c in rel:
            us = cols[c]; cnt = collections.Counter(vv for u in us for vv in u["f"])
            cons.append(dict(lo=int(sorted(u["lo"] for u in us)[len(us) // 2]), hi=int(sorted(u["hi"] for u in us)[len(us) // 2]),
                             f={vv for vv, k in cnt.items() if k > len(us) / 2}, rel=c, n=len(us)))
        ci = rel.index(0)
        # (가) 이웃으로
        (s1, sl, ap, dirn), (s2, sl2, ap2, _) = place_on_map(cons, ci)
        amb = list(LAST_AMB); amb_z = sorted({z for _, _, z in amb})
        placed = s1 >= SMIN and s1 - s2 >= MMIN
        how, where, verdict = "", "", ""
        if placed:
            zz = zone(sl, ap); how = "가_이웃"
            where = "slot%s ap%d (%s)" % (sl, ap, zz)
            if ap in mk["car"].get(sl, set()): verdict = "기대 유닛"
            elif (sl, ap) in FRSET[(p, a)]: verdict = "지도엔 있음(HOT 등)·같은 구역" if zz == mk["zone"] else "지도엔 있음·다른 구역"
            elif zz == mk["zone"]: verdict = "같은 구역·다른 유닛"
            else: verdict = "다른 구역"
        # (나) 끝에서 잰 거리
        dd = [d for d in DIRECT if d["variant"] == v and d["read"] in reads]
        dist_d = sorted(d["dist"] for d in dd if d["end"] == "원위"); dist_p = sorted(d["dist"] for d in dd if d["end"] == "근위")
        def zone_from_dist():
            zs = set()
            for d in dist_d: zs.add("D" if d <= min(KD.values()) else ("D?" if d <= max(KD.values()) else "not-D"))
            for d in dist_p: zs.add("P" if d <= min(KP.values()) else ("P?" if d <= max(KP.values()) else "not-P"))
            return zs
        zd = zone_from_dist()
        if not placed:
            if zd:
                how = "나_끝거리"; where = "원위끝에서 %s · 근위끝에서 %s" % (dist_d, dist_p)
                z = mk["zone"]
                contra = (z == "D" and ("not-D" in zd or "P" in zd)) or (z == "P" and ("not-P" in zd or "D" in zd)) or (z == "M" and ("D" in zd or "P" in zd))
                verdict = "끝 거리가 다른 구역" if contra else "끝 거리가 구역과 맞음(또는 판정 범위 밖)"
            else:
                how = "다_미정"; verdict = "위치 미정"
        novel_center = sorted("%d%s" % vv for vv in cons[ci]["f"] if vv not in CNT)
        # v4 : 칸 단위 정렬로 3 분류
        bestp = None
        for s_ in sorted(FR):
            for dirn_ in (1, -1):
                F_ = FR[s_] if dirn_ == 1 else FR[s_][::-1]
                sc_, pr_, S_ = sw_path(cons, F_)
                if bestp is None or sc_ > bestp[0]: bestp = (sc_, s_, F_, pr_, S_)
        psc, ps, PF, ppairs, PS = bestp; pm = dict(ppairs)
        status = []
        for k_, u_ in enumerate(cons):
            nov_ = [vv for vv in u_["f"] if vv not in CNT] if u_["n"] >= 2 else []
            j_ = pm.get(k_)
            if j_ is None: st_ = "불일치(짝없음)" if (u_["f"] & ALLF or nov_) else "정보없음"
            elif PS[k_][j_] >= 0.5: st_ = "맞음"
            elif informative(u_, PF[j_][1]): st_ = "불일치"
            else: st_ = "정보없음"
            if nov_ and st_ == "맞음": st_ = "맞음+새변이"
            status.append((st_, PF[j_][0] if j_ is not None else None))
        mst, map_ = status[ci]
        nmis = sum(1 for s_, _ in status if s_.startswith("불일치") or s_ == "맞음+새변이")
        if mst == "맞음":
            cls = "A 원래 자리" if (map_ in mk["car"].get(ps, set()) or (ps, map_) in FRSET[(p, a)]) else "B 조립본 다른 자리"
            zf = zone(ps, map_); zhow = "그 지도 유닛"
            anchor = "slot%s ap%d" % (ps, map_)
        else:
            cls = "C 조립본 밖"
            L_ = next(((k_, status[k_][1]) for k_ in range(ci - 1, -1, -1) if status[k_][0] == "맞음"), None)
            R_ = next(((k_, status[k_][1]) for k_ in range(ci + 1, len(cons)) if status[k_][0] == "맞음"), None)
            zs_ = [zone(ps, x[1]) for x in (L_, R_) if x]
            anchor = "왼쪽 %s · 오른쪽 %s" % ("slot%s ap%d(%d칸)" % (ps, L_[1], ci - L_[0]) if L_ else "없음", "slot%s ap%d(%d칸)" % (ps, R_[1], R_[0] - ci) if R_ else "없음")
            if L_ and R_ and zs_[0] == zs_[1]: zf, zhow = zs_[0], "양쪽 닻"
            elif L_ and R_: zf, zhow = "%s|%s" % tuple(zs_), "닻이 두 구역에 걸침"
            elif zs_: zf, zhow = zs_[0], "한쪽 닻"
            elif zd: zf, zhow = ",".join(sorted(zd)), "끝 거리"
            else: zf, zhow = "미정", "닻·끝 없음"
        zres = "같은 구역" if zf == mk["zone"] else ("미정" if zf == "미정" else ("걸침" if "|" in zf or "," in zf or "?" in zf else "다른 구역"))
        if zhow == "닻이 두 구역에 걸침" and mk["zone"] in zf.split("|"): zres = "걸침(표지 구역 포함)"
        RES.append(dict(variant=v, mzone=mk["zone"], n_read=len(reads), n_item=len(g), src=dict(collections.Counter(items[x]["src"] for x in g)),
                        cons_len=len(cons), how=how, where=where, verdict=verdict, s1=round(s1, 2), s2=round(s2, 2),
                        second="slot%s ap%s (%s)" % (sl2, ap2, zone(sl2, ap2) if sl2 else ""), dist_distal=dist_d, dist_prox=dist_p,
                        novel_center=novel_center, reads=reads, best="slot%s ap%s (%s)" % (sl, ap, zone(sl, ap)), amb_n=len(amb), amb_zones=amb_z,
                        amb=["s%s:%d" % (a_, b_) for a_, b_, _ in amb][:30], idx=[[items[x]["read"], items[x]["k"]] for x in g],
                        cls=cls, zone_final=zf, zone_how=zhow, zone_res=zres, anchor=anchor, marker_col=mst, n_mis=nmis,
                        path_score=round(psc, 2), cols=["%s%s" % (s_[:2], "" if a_ is None else a_) for s_, a_ in status],
                        n_pb=sum(1 for q_ in reads if READS[q_]["plat"] == "PB")))
    log("%-6s %s : read 창 %d · 무리 %d (2 read 이상 %d)" % (v, mk["zone"], n, len(groups), sum(1 for g in groups if len({items[x]['read'] for x in g}) >= 2)))
json.dump(dict(direction_check={"|".join(k): c for k, c in chk.items()}, results=RES, direct=DIRECT),
          open("%s/%s_markerasm_v5.json" % (W, TAG), "w"), ensure_ascii=False, indent=1, default=str)
tot = collections.Counter((r["verdict"], r["n_read"] >= 2) for r in RES)
log("완료", {("%s|%s" % (k[0], "2read+" if k[1] else "1read")): c for k, c in sorted(tot.items())})
