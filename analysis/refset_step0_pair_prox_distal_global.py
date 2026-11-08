# -*- coding: utf-8 -*-
# 근위↔원위 짝 전역 배정
#   단계 : ① 하플로타입 확정

# -*- coding: utf-8 -*-
# flank read 하플로타입 전역 배정 (2026-09-05)
#
#  규칙은 consensus_v2.py 의 정본을 따른다:
#    "근위<->원위 짝은 관통 read 가 있으면 그 분자가 정한다 (번호 가정 금지)"
#
#  근위/원위 hap 번호는 서로 독립이다.
#    실측: 양쪽 다 가진 read 575 개 중 번호 일치 45.4% = 무작위.
#    번호로 join 하면 절반이 틀린다. 반드시 분자(관통) 또는 부모로 이어야 한다.
#
#  우선순위
#    ① 부모 유래 (~/parent_hap_full.tsv) — 전역 일관 (양쪽 read 13/13 일치 실측)
#    ② 관통(SPAN) read 가 정한 prox<->dist 짝
#       소거법: 하플로타입은 2 개뿐이므로 한쪽 연결만 알아도 나머지가 결정된다
#         (앞 판에서 이걸 빠뜨려 139 개체가 "일부만 확정" 으로 남았다)
#    ③ 둘 다 없으면 미확정 -> 하플로타입 서열로 잇는다 (다음 단계)
#
#  전역 라벨은 근위 번호를 기준으로 삼는다 (prox hap 1/2 를 그대로 전역 1/2 로).
import io
import os
import collections

S = os.environ["RDNA_S"]
NL = "\n"


def load(p):
    R = collections.defaultdict(dict)
    for i, L in enumerate(io.open(os.path.join(S, p), encoding="utf8")):
        if i == 0:
            continue
        f = L.rstrip(NL).split("\t")
        if len(f) < 4:
            continue
        q = f[2].split("|")
        rid = q[1] if len(q) > 1 else f[2]
        R[(f[0], rid)][f[1]] = f[3]
    return R


V = load("hap_reads_v5.tsv")
P = load("parent_hap_full.tsv")

SPAN = set()
for i, L in enumerate(io.open(os.path.join(S, "units_all.tsv"), encoding="utf8")):
    if i == 0:
        continue
    f = L.rstrip(NL).split("\t")
    if len(f) >= 3 and f[2] == "SPAN":
        SPAN.add((f[0], f[1]))

# ---- 개체별 prox<->dist 연결 ----
LINK = {}          # smp -> {prox_hap: dist_hap}
CONF = set()
raw = collections.defaultdict(collections.Counter)
for (s, r), v in V.items():
    if (s, r) not in SPAN or len(v) < 2:
        continue
    raw[s][(v["prox"], v["dist"])] += 1
for s, c in raw.items():
    d = {}
    bad = False
    for (ph, dh), n in c.most_common():
        if ph in d and d[ph] != dh:
            bad = True
        elif dh in d.values() and d.get(ph) != dh:
            bad = True
        else:
            d[ph] = dh
    if bad:
        CONF.add(s)
        continue
    # 소거법 — 하플로타입은 2 개뿐이다
    if len(d) == 1:
        ph, dh = list(d.items())[0]
        op = "2" if ph == "1" else "1"
        od = "2" if dh == "1" else "1"
        d[op] = od
    LINK[s] = d

# ---- 부모 라벨을 근위 기준으로 맞춘다 ----
# 부모 라벨은 전역이지만 번호가 v5 와 다를 수 있다. 개체·side 별 대응을 다수결로 만든다.
PMAP = {}
for side in ("prox", "dist"):
    m = collections.defaultdict(collections.Counter)
    for (s, r), v in P.items():
        w = V.get((s, r))
        if not w or side not in v or side not in w:
            continue
        m[s][(w[side], v[side])] += 1
    for s, c in m.items():
        d = collections.defaultdict(collections.Counter)
        for (a, b), n in c.items():
            d[a][b] += n
        mp = {}
        pur = []
        for a in d:
            b, n = d[a].most_common(1)[0]
            mp[a] = b
            pur.append(n / float(sum(d[a].values())))
        PMAP[(s, side)] = (mp, min(pur) if pur else 0.0, sum(c.values()))

# 부모로 연결 : prox 라벨 -> 부모라벨 -> dist 라벨
PLINK = {}
for s in set(k[0] for k in PMAP):
    a = PMAP.get((s, "prox"))
    b = PMAP.get((s, "dist"))
    if not a or not b or a[1] < 0.8 or b[1] < 0.8 or a[2] < 4 or b[2] < 4:
        continue
    inv = {}
    ok = True
    for dl, pl in b[0].items():
        if pl in inv:
            ok = False
        inv[pl] = dl
    if not ok:
        continue
    d = {}
    for pl, par in a[0].items():
        if par in inv:
            d[pl] = inv[par]
    if len(d) == 1:
        ph, dh = list(d.items())[0]
        d["2" if ph == "1" else "1"] = "2" if dh == "1" else "1"
    if d:
        PLINK[s] = d

OUT = io.open(os.path.join(S, "HAP_READ_GLOBAL.tsv"), "w", encoding="utf8", newline=NL)
OUT.write("sample\tread\tside\tside_hap\tglobal_hap\tevidence\n")
n = collections.Counter()
per = collections.defaultdict(collections.Counter)
for (s, r), v in sorted(V.items()):
    lk = LINK.get(s)
    pk = PLINK.get(s)
    for side, h in sorted(v.items()):
        g = "-"
        ev = "none"
        use = lk if lk else pk
        src = "span" if lk else ("parent" if pk else None)
        if use:
            if side == "prox":
                g, ev = h, src
            else:
                inv = {dv: pv for pv, dv in use.items()}
                if h in inv:
                    g, ev = inv[h], src
        n[ev] += 1
        per[s][ev] += 1
        OUT.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (s, r, side, h, g, ev))
OUT.close()

smps = set(k[0] for k in V)
print("=" * 78)
print("  flank read 전역 하플로타입 배정 — 개체 %d · read x side %d" %
      (len(smps), sum(n.values())))
print("=" * 78)
print()
print("  개체별 연결 근거")
print("     관통(SPAN)으로 연결 : %d  (모순 %d 제외)" % (len(LINK), len(CONF)))
print("     부모로 연결        : %d  (관통 없는 개체 중 %d)" %
      (len(PLINK), len(set(PLINK) - set(LINK))))
print("     둘 중 하나       : %d / %d" % (len(set(LINK) | set(PLINK)), len(smps)))
print("     미확정           : %d" % len(smps - set(LINK) - set(PLINK)))
print()
print("  read x side 배정")
for k, v in n.most_common():
    print("     %-8s %6d (%5.1f%%)" % (k, v, 100.0 * v / sum(n.values())))
print()
st = collections.Counter()
for s in smps:
    c = per[s]
    st["전부 확정" if c["none"] == 0 else
       ("전부 미확정" if c["none"] == sum(c.values()) else "일부만")] += 1
print("  개체 단위 :", dict(st))
un = sorted(s for s in smps if per[s]["none"] == sum(per[s].values()))
io.open(os.path.join(S, "hap_unlinked.txt"), "w", newline=NL).write(NL.join(un) + NL)
print("  저장 : HAP_READ_GLOBAL.tsv · hap_unlinked.txt (%d 개체)" % len(un))
