# -*- coding: utf-8 -*-
# tier1 우선 → tier2 · 중복 제거 · 개체당 최대 2
#   단계 : 참조 확정
#   Method document : analysis/METHOD_CONSENSUS_MAP_G8_REBUILD.md

# -*- coding: utf-8 -*-
# 하플로타입별 참조 확정 — 메틸 맵을 붙일 자를 정한다
#
#  tier 규칙 (사용자 확정 2026-09-02):
#    tier1 = 관통 read 의 array 구간      (있으면 이것)
#    tier2 = hifiasm 재조립 array          (관통 read 가 없을 때)
#    tier3 = 둘 다 없음                    (참조를 만들 수 없다)
#
#  hifiasm 이 여러 contig 를 낼 때 고르는 법
#    HG002 실측: hifiasm 이 3 개를 냈고 hap1=68(정답) · hap2=77(가짜) · pri=103(정답).
#    ⇒ 라벨(hap1/hap2)로 고르면 틀린다.
#    가짜 77 은 진짜 103 과 99.84% 로 겹치는 부분 중복본이었다.
#    ⇒ 규칙: 길이 내림차순으로 보되, 이미 고른 것과 서열이 크게 겹치면(>=90% 덮임) 버린다.
#    개체당 최대 2 개(하플로타입 2 개)만 남긴다.
import io
import os
import sys
import subprocess
import collections

W = os.path.expanduser("~/refbuild")
os.makedirs(W, exist_ok=True)
T1 = os.path.expanduser("~/tier1ref")
HA = os.path.expanduser("~/hiasm/arr")
TH = os.environ.get("TH", "4")
TAG = os.environ.get("TAG", "x")
NL = "\n"


def read_fa1(p):
    L = io.open(p, encoding="utf8").read().split(NL)
    return L[0][1:], "".join(L[1:])


def dedup(cands):
    """길이 내림차순 · 이미 고른 것과 크게 겹치면 버린다"""
    keep = []
    for nm, sq, src in sorted(cands, key=lambda z: -len(z[1])):
        dup = False
        for knm, ksq, ksrc in keep:
            a = os.path.join(W, "a.fa")
            b = os.path.join(W, "b.fa")
            io.open(a, "w", newline=NL).write(">a" + NL + sq + NL)
            io.open(b, "w", newline=NL).write(">b" + NL + ksq + NL)
            paf = subprocess.run(["minimap2", "-x", "asm10", "-t", TH,
                                  "--secondary=no", b, a],
                                 capture_output=True, text=True)
            cov = 0
            iv = []
            for L in paf.stdout.splitlines():
                f = L.split("\t")
                if len(f) >= 11:
                    iv.append((int(f[2]), int(f[3])))
            iv.sort()
            mg = []
            for s0, e0 in iv:
                if mg and s0 <= mg[-1][1]:
                    mg[-1][1] = max(mg[-1][1], e0)
                else:
                    mg.append([s0, e0])
            cov = sum(e - s for s, e in mg)
            if cov >= 0.90 * len(sq):
                dup = True
                break
        if not dup:
            keep.append((nm, sq, src))
        if len(keep) >= 2:
            break
    return keep


OUT = io.open(os.path.join(W, "refs_%s.tsv" % TAG), "a", encoding="utf8", newline=NL)
for smp in sys.argv[1:]:
    # tier1 후보
    t1 = []
    for f in sorted(os.listdir(T1)) if os.path.isdir(T1) else []:
        if f.startswith("t1_%s_" % smp) and f.endswith(".fa"):
            nm, sq = read_fa1(os.path.join(T1, f))
            t1.append((nm, sq, "tier1"))
    # tier2 후보
    t2 = []
    if os.path.isdir(HA):
        for f in sorted(os.listdir(HA)):
            if f.startswith("hiarr_%s_" % smp) and f.endswith(".fa"):
                nm, sq = read_fa1(os.path.join(HA, f))
                t2.append((nm, sq, "tier2"))
    if not t1 and not t2:
        OUT.write("%s\t-\tNONE\t-\t-%s" % (smp, NL))
        OUT.flush()
        continue
    # tier1 우선. tier1 이 하나뿐이면 나머지 한 자리는 tier2 에서 채운다.
    picked = dedup(t1) if t1 else []
    if len(picked) < 2 and t2:
        more = dedup(picked + t2)
        for x in more:
            if len(picked) >= 2:
                break
            if not any(x[0] == y[0] for y in picked):
                picked.append(x)
    od = os.path.join(W, "ref")
    os.makedirs(od, exist_ok=True)
    for k, (nm, sq, src) in enumerate(picked):
        p = os.path.join(od, "%s_h%d.fa" % (smp, k + 1))
        io.open(p, "w", newline=NL).write(
            ">%s_h%d|%s|%s|%d%s%s%s" % (smp, k + 1, src, nm, len(sq), NL, sq, NL))
        OUT.write("%s\th%d\t%s\t%s\t%d%s" % (smp, k + 1, src, nm, len(sq), NL))
    OUT.flush()
    sys.stderr.write("  %s : tier1후보 %d · tier2후보 %d -> 채택 %s\n"
                     % (smp, len(t1), len(t2),
                        " ".join("%s(%d bp)" % (x[2], len(x[1])) for x in picked)))
    sys.stderr.flush()
OUT.close()
print("done")
