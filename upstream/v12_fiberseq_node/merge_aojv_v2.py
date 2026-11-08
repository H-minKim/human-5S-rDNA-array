# -*- coding: utf-8 -*-
"""array only 배치 + junction 배치를 규칙대로 합친다 (2026-09-16)

  sort -u 로 합치면 안 된다 — 같은 read 에 두 출처가 다른 슬롯을 주면 두 줄이 다 살아남고
    ao_g8cut 이 먼저 읽은 줄(사전순이라 늘 슬롯 "1")을 잡는다. 규칙이 아니라 우연이다.
    (measured : of 803 overlapping reads 6 conflict, and all 6 are cut to slot 1 in the output)

  충돌 해소 규칙 (2026-09-16 검정 근거)
    ① 판별 자리 < 10 인 투표는 배제 — 라벨 뒤섞기 귀무에서 자리<5 는 위양성 14.3%,
       자리 10+ 는 0.0% (자리 4 개짜리 표차 1.000 은 1/7 로 우연히 난다)
    ② 남은 충돌은 array only 사슬을 따른다 — 충돌 read 는 전부 유닛 >=14 로
       사슬 게이트(K=12)가 검정된 범위 안이다. 투표는 사슬이 안 되는 구간을 메우려 만든 것이지
       사슬을 뒤집으라고 만든 것이 아니다.
    ③ 결과와 근거를 반드시 기록한다 (조용히 정해지는 일이 없도록)

  쓰기 : merge_aojv.py <array only tsv…> --jx <junction tsv> --vote <slotvote tsv> --out <합본> --log <기록>
"""
import io
import os
import sys
import collections

NL = "\n"
TAB = "\t"
# the threshold is the environment variable MINSITE, default 10; the value used is 5 :
#   귀무 위양성 자리<5 14.3% · 5-9 2.4% · 10+ 0.0% -> 5-9 는 살린다 (되살린 것 중 약 3% 오배정 가능)
MINSITE = int(os.environ.get("MINSITE", "10"))
# read 열쇠 길이 : ONT 16 · PacBio 0(전체) — 앞 16 자가 movie 라 전건 충돌
KEYLEN = int(os.environ.get("KEYLEN", "16"))


def rkey(x):
    return x[:KEYLEN] if KEYLEN else x

args = sys.argv[1:]
ao_files, jx, vote, out, log = [], None, None, None, None
i = 0
while i < len(args):
    a = args[i]
    if a == "--jx":
        i += 1
        jx = args[i]
    elif a == "--vote":
        i += 1
        vote = args[i]
    elif a == "--out":
        i += 1
        out = args[i]
    elif a == "--log":
        i += 1
        log = args[i]
    else:
        ao_files.append(a)
    i += 1

# 투표의 판별 자리 수 : (개체, read16) -> 자리수
NSITE = {}
if vote and os.path.exists(vote):
    for k, L in enumerate(io.open(vote, encoding="utf8")):
        if k == 0:
            continue
        f = L.rstrip(NL).split(TAB)
        if len(f) < 8:
            continue
        NSITE[(f[0], rkey(f[1]))] = int(f[4])


def load(paths, tag):
    """(개체, read16) -> (슬롯, [행…])"""
    D = collections.defaultdict(lambda: [None, []])
    for p in paths:
        if not os.path.exists(p):
            continue
        for L in io.open(p, encoding="utf8"):
            f = L.rstrip(NL).split(TAB)
            if len(f) < 5:
                continue
            k = (f[0], rkey(f[1]))
            if D[k][0] is None:
                D[k][0] = f[2]
            D[k][1].append(L.rstrip(NL))
    return D


AO = load(ao_files, "AO")
JX = load([jx], "JX") if jx else {}

lg = io.open(log, "w", encoding="utf8", newline=NL) if log else sys.stderr
lg.write("# 합치기 기록 (충돌 해소 규칙 : 자리<%d 투표 배제 · 남은 충돌은 array only)\n" % MINSITE)
lg.write("sample\tread16\tao_slot\tjx_slot\tn_site\t채택\t사유\n")

T = collections.Counter()
rows = []
for k in sorted(set(list(AO.keys()) + list(JX.keys()))):
    a = AO.get(k)
    b = JX.get(k)
    if a and not b:
        rows.extend(a[1])
        T["array only 만"] += 1
        continue
    if b and not a:
        ns = NSITE.get(k, 0)
        if ns < MINSITE:
            T["junction 만인데 자리부족 → 버림"] += 1
            lg.write("%s\t%s\t-\t%s\t%d\t버림\t자리<%d\n" % (k[0], k[1], b[0], ns, MINSITE))
            continue
        rows.extend(b[1])
        T["junction 만"] += 1
        continue
    # 둘 다 있다
    if a[0] == b[0]:
        # 같은 슬롯 — array only 를 본으로 쓰고 junction 이 더 찾은 유닛을 보탠다
        seen = set(x.split(TAB)[3] for x in a[1])
        extra = [x for x in b[1] if x.split(TAB)[3] not in seen]
        rows.extend(a[1])
        rows.extend(extra)
        T["겹침·같은슬롯"] += 1
        T["  junction 이 더한 유닛"] += len(extra)
        continue
    ns = NSITE.get(k, 0)
    rows.extend(a[1])                      # 규칙 ② array only 를 따른다
    T["슬롯충돌 → array only 채택"] += 1
    lg.write("%s\t%s\t%s\t%s\t%d\tAO\t%s\n"
             % (k[0], k[1], a[0], b[0], ns,
                ("자리<%d 라 투표 배제" % MINSITE) if ns < MINSITE else "사슬 범위(유닛>=12) 우선"))

with io.open(out, "w", encoding="utf8", newline=NL) as o:
    for r in rows:
        o.write(r + NL)

lg.write("\n# 요약\n")
for kk, vv in sorted(T.items()):
    lg.write("#   %-34s %d\n" % (kk, vv))
lg.write("#   %-34s %d\n" % ("최종 행", len(rows)))
if log:
    lg.close()
print("  " + " · ".join("%s %d" % kv for kv in sorted(T.items())) + " · 행 %d" % len(rows))
