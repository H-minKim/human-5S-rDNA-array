# -*- coding: utf-8 -*-
"""read 의 하플로타입을 정본 변이 자리에서만 투표로 정한다 (2026-09-15 사용자 설계)

  원리 : read 의 나머지 불일치는 오류로 본다. 정본 두 슬롯이 서로 다른 자리만 진짜 변이로 쓰고,
          그 자리에서 read 가 어느 슬롯 염기를 갖는지 센다.
          (array 배치의 |S∩V|-|S-V| 와 같은 원리 — 정체율 비교는 ONT 5% 오류에 묻힌다)

  flank 만 쓴다 : array 는 탠덤이라 두 슬롯이 같은 자리에 다 붙어 판별력이 없다
                   (조각 합산 실험에서 정확도 100% -> 84.7% 로 떨어졌다)

  절차
    ① 슬롯1 flank <-> 슬롯2 flank 를 정렬해 다른 자리(판별 자리)를 뽑는다
    ② read 를 두 슬롯에 정렬해, 판별 자리에서 read 염기를 읽는다
    ③ 슬롯1 염기와 같으면 1표, 슬롯2 염기와 같으면 2표. 둘 다 아니면 오류로 버린다
    ④ 표차가 크면 배정, 아니면 미배정

  출력 : sample read units flank n_site v1 v2 margin verdict
"""
import io
import os
import sys
import subprocess
import collections
import re

CS = re.compile(r"([:*+\-~])([A-Za-z0-9\[\]]+)")
HSN = os.path.expanduser(os.environ.get("HSN", "~/hap_seq_node.tsv"))
CMPD = os.path.expanduser(os.environ.get("CMPD", "~/cmp4d_flank"))
BOUND = os.path.expanduser(os.environ.get("BOUND", "~/fs41/out"))
READS = os.path.expanduser(os.environ.get("READS", "~/fs41/out"))
TH = os.environ.get("TH", "6")
MINV = int(os.environ.get("MINV", "3"))        # 최소 판별 자리 수
OUT = sys.argv[2] if len(sys.argv) > 2 else "/tmp/slotvote.tsv"

SEQF = {}
for x in io.open(HSN, encoding="utf8"):
    f = x.rstrip("\n").split("\t")
    if len(f) >= 3:
        SEQF[(f[0], f[1])] = f[2]


def load(p):
    hdr = ""
    buf = []
    for L in io.open(p, encoding="utf8"):
        if L.startswith(">"):
            hdr = L.strip()
        else:
            buf.append(L.strip())
    s = "".join(buf).upper()
    d = dict(t.split("=", 1) for t in hdr.split("|")[1:] if "=" in t)
    return s, int(d.get("flank5", 0)), int(d.get("flank3", 0))


def paf(ref, q, preset, extra=()):
    r = subprocess.run(["minimap2", "-x", preset, "-c", "--cs", "-t", TH] + list(extra) + [ref, q],
                       capture_output=True, text=True)
    return r.stdout.splitlines()


o = io.open(OUT, "w", encoding="utf8", newline="\n")
o.write("sample\tread\tunits\tflank\tn_site\tv1\tv2\tmargin\tverdict\n")

for smp in [x.strip() for x in io.open(sys.argv[1], encoding="utf8") if x.strip()]:
    try:
        s1, a5, a3 = load(os.path.join(CMPD, SEQF[(smp, "1")]))
        s2, b5, b3 = load(os.path.join(CMPD, SEQF[(smp, "2")]))
    except Exception:
        continue

    # ── ① 판별 자리 : 두 슬롯 flank 를 서로 정렬해 다른 자리를 뽑는다
    fa = "/tmp/%s.f1.fa" % smp
    fb = "/tmp/%s.f2.fa" % smp
    io.open(fa, "w", newline="\n").write(">s1_5\n%s\n>s1_3\n%s\n" % (s1[:a5], s1[len(s1) - a3:]))
    io.open(fb, "w", newline="\n").write(">s2_5\n%s\n>s2_3\n%s\n" % (s2[:b5], s2[len(s2) - b3:]))
    # 판별 자리를 (슬롯1 전체좌표, 염기1, 슬롯2 전체좌표, 염기2) 로 모은다
    SITE = []
    for L in paf(fa, fb, "asm20", ("--secondary=no",)):
        f = L.split("\t")
        if len(f) < 12:
            continue
        qn, tn = f[0], f[5]
        if qn[3:] != tn[3:]:                     # 5'는 5'끼리, 3'는 3'끼리만
            continue
        cs = ""
        for x in f[12:]:
            if x.startswith("cs:Z:"):
                cs = x[5:]
        tp = int(f[7])
        qp = int(f[2])
        rev = (f[4] == "-")
        if rev:
            continue                             # 방향 뒤집힌 정렬은 쓰지 않는다
        for g_ in CS.finditer(cs):
            op, v = g_.group(1), g_.group(2)
            if op == ":":
                n_ = int(v)
                tp += n_
                qp += n_
            elif op == "*":
                # v = "ab" : a=참조(슬롯1) 염기, b=질의(슬롯2) 염기
                t_abs = tp if tn.endswith("_5") else (len(s1) - a3 + tp)
                q_abs = qp if qn.endswith("_5") else (len(s2) - b3 + qp)
                SITE.append((t_abs, v[0].upper(), q_abs, v[1].upper()))
                tp += 1
                qp += 1
            elif op == "-":
                tp += len(v)
            elif op == "+":
                qp += len(v)
    for p in (fa, fb):
        try:
            os.remove(p)
        except OSError:
            pass
    if len(SITE) < MINV:
        sys.stderr.write("  %-10s 판별자리 %d — 건너뜀\n" % (smp, len(SITE)))
        continue
    S1 = dict((t, (b1, b2)) for (t, b1, _q, b2) in SITE)
    S2 = dict((q, (b1, b2)) for (_t, b1, q, b2) in SITE)

    # ── ② read 를 두 슬롯에 정렬
    ref = "/tmp/%s.two.fa" % smp
    io.open(ref, "w", newline="\n").write(">%s_h1\n%s\n>%s_h2\n%s\n" % (smp, s1, smp, s2))
    bp = os.path.join(BOUND, "%s.bound.tsv" % smp)
    if not os.path.exists(bp):
        continue
    want = {}
    fh = io.open(bp, encoding="utf8")
    hh = fh.readline().rstrip("\n").split("\t")
    jx = dict((n, i) for i, n in enumerate(hh))
    for L in fh:
        r = L.rstrip("\n").split("\t")
        if len(r) < 10:
            continue
        want[r[jx["read"]]] = (int(float(r[jx["units"]])),
                               max(int(r[jx["flankL"]]), int(r[jx["flankR"]])))
    sub = "/tmp/%s.sub.fa" % smp
    n = 0
    with io.open(sub, "w", newline="\n") as g:
        it = io.open(os.path.join(READS, "%s.hit.fq" % smp), encoding="utf8")
        while True:
            h = it.readline()
            if not h:
                break
            sq = it.readline().strip()
            it.readline()
            it.readline()
            rid = h[1:].split()[0].split("\t")[0]
            if rid in want:
                g.write(">%s\n%s\n" % (rid, sq))
                n += 1
    if n == 0:
        continue

    V = collections.defaultdict(lambda: [0, 0, 0])     # read -> [자리수, 1표, 2표]
    for L in paf(ref, sub, "map-ont", ("-N", "5", "-p", "0.1")):
        f = L.split("\t")
        if len(f) < 12:
            continue
        rid = f[0]
        sl = f[5].rsplit("_h", 1)[1]
        TAB = S1 if sl == "1" else S2
        cs = ""
        for x in f[12:]:
            if x.startswith("cs:Z:"):
                cs = x[5:]
        tp = int(f[7])
        for g_ in CS.finditer(cs):
            op, v = g_.group(1), g_.group(2)
            if op == ":":
                n_ = int(v)
                for t in range(tp, tp + n_):
                    hit = TAB.get(t)
                    if hit:
                        base = hit[0] if sl == "1" else hit[1]   # 그 슬롯의 염기와 같다
                        V[rid][0] += 1
                        V[rid][1 if base == hit[0] else 2] += 1
                tp += n_
            elif op == "*":
                hit = TAB.get(tp)
                if hit:
                    rb = v[1].upper()                            # read 염기
                    V[rid][0] += 1
                    if rb == hit[0]:
                        V[rid][1] += 1
                    elif rb == hit[1]:
                        V[rid][2] += 1
                tp += 1
            elif op == "-":
                tp += len(v)
    for p in (ref, sub):
        try:
            os.remove(p)
        except OSError:
            pass

    nw = 0
    for rid, (ns, v1, v2) in V.items():
        if v1 + v2 < MINV:
            continue
        u, fk = want.get(rid, (0, 0))
        mg = abs(v1 - v2) / float(v1 + v2)
        o.write("%s\t%s\t%d\t%d\t%d\t%d\t%d\t%.4f\t%s\n"
                % (smp, rid, u, fk, v1 + v2, v1, v2, mg,
                   "SEP" if mg >= 0.60 else "AMB"))
        nw += 1
    sys.stderr.write("  %-10s 판별자리 %d · 투표된 read %d\n" % (smp, len(SITE), nw))
o.close()
