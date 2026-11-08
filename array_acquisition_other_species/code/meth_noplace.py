# -*- coding: utf-8 -*-
# 배치 없이 메틸을 잰다 — read 를 array 에 배치하지 않고 5S 메틸 상태를 구한다
#
#   왜 : array 는 탠덤이라 read 가 mapq 0 로 놓인다(금화조 배치 read 의 69~72%).
#        "어느 유닛에 놓였나"가 못 미더우면 결론도 못 미덥다.
#        그런데 메틸 호출은 read 자신에 실려 있다(MM/ML). 참조에 놓을 필요가 없다.
#        필요한 것은 "read 의 이 구간이 5S 인가, 어느 종류의 5S 인가" 뿐이고,
#        그것은 서열 대조로 바로 안다. array 좌표도 mapq 도 쓰지 않는다.
#
#   방법
#     1) read 를 WIN bp 창으로 쪼갠다 (정렬·좌표 일절 안 씀)
#     2) 창을 참조 유닛 종류들에만 정렬해 최고 점수 종류로 딱지를 붙인다(문맥판 사용)
#     3) 그 창 안의 C 변형에 대해 MM/ML 로 메틸 비율을 낸다 (mod_qual >= 0.5 를 메틸로)
#        문맥은 read 서열에서 직접 판정한다 — 다음 염기가 G 면 CG, 그 다음이 G 면 CHG, 아니면 CHH.
#        (식물은 세 문맥을 다 부르므로 섞으면 CG 0.84 가 0.5 로 희석된다. 실측으로 드러난 함정)
#        단 호출 모델이 CpG 전용이면(유인원·금화조·제브라피시) MM 에 실린 C 는 전부 CpG 다.
#          그때 read 서열로 문맥을 나누면 read 오류만 나누는 꼴이 된다(고릴라 실측 : CHH 로 18% 샘).
#          그래서 MM 자리의 '다음 염기가 G' 비율로 모델 종류를 자료에서 판정하고,
#          0.75 이상이면 전부 CG 로 둔다.
#     4) 종류마다 · read 마다 비율을 모은다
#   BAM 은 MM/ML 을 꺼내는 그릇으로만 쓴다 — 정렬 좌표·mapq·유닛 번호는 읽지 않는다.
#
#   산출 : <out>/<species>__<sample>.noplace_window.tsv  (창 한 줄)
#          <out>/<species>__<sample>.noplace_read.tsv    (read 한 줄 · 종류별)
#   인자 : 종키 샘플 bam 참조유닛종류.fa 출력폴더 [WIN=500] [MINID=0.70] [CTX=120]
# MIT License · Copyright (c) 2026 Hungmin Kim
import io
import os
import sys
import csv
import subprocess
import collections
import pysam

NL = chr(10)
TAB = chr(9)
SP, SAMPLE, BAM, REFT, OUT = sys.argv[1:6]
WIN = int(sys.argv[6]) if len(sys.argv) > 6 else 500
MINID = float(sys.argv[7]) if len(sys.argv) > 7 else 0.70
CTX = int(sys.argv[8]) if len(sys.argv) > 8 else 120
# 호출 문턱 : modkit 은 기본으로 하위 10% 신뢰도를 버리고 남은 것에 0.5 를 적용한다.
#   여기서 0.5 를 그대로 쓰면 애매한 호출이 전부 '메틸' 로 들어가 저메틸 문맥(CHH)이 부풀려진다.
#   같은 잣대로 맞추려면 문턱을 올린다(감도로 여러 값을 같이 본다).
MODTHR = float(sys.argv[9]) if len(sys.argv) > 9 else 0.5
os.makedirs(OUT, exist_ok=True)
TMP = os.path.join(OUT, "_t_%s_%s" % (SP, SAMPLE))
os.makedirs(TMP, exist_ok=True)

# 참조 유닛 종류 → 문맥판 (탠덤이라 창이 유닛 경계를 걸친다)
refs = collections.OrderedDict()
name = None
for L in io.open(REFT, encoding="utf8"):
    L = L.rstrip(NL)
    if L.startswith(">"):
        name = L[1:].split()[0]
        refs[name] = []
    elif name:
        refs[name].append(L.strip())
refs = collections.OrderedDict((k, "".join(v).upper()) for k, v in refs.items())
cp = os.path.join(TMP, "ctx.fa")
with io.open(cp, "w", encoding="utf8", newline=NL) as fh:
    for k, s in refs.items():
        r = (s * 4)[:max(len(s) * 2, WIN + 2 * CTX)]   # 탠덤으로 이어 붙인 문맥판
        fh.write(">%s%s%s%s" % (k, NL, r, NL))

bam = pysam.AlignmentFile(BAM, "rb", check_sq=False)
wf = io.open(os.path.join(TMP, "win.fa"), "w", encoding="utf8", newline=NL)
META = []
MOD = {}
nread = 0
SEQ = {}
for rd in bam.fetch(until_eof=True):
    if rd.is_secondary or rd.is_supplementary:
        continue
    q = rd.query_sequence
    if not q:
        continue
    try:
        mb = rd.modified_bases
    except (ValueError, KeyError):
        mb = None
    if not mb:
        continue
    nread += 1
    # CpG 5mC : (canonical, strand, code) 키에서 C+m / C+h 를 모은다
    call = {}
    for k, v in mb.items():
        code = k[2]
        if code not in ("m", "h", 27551):
            continue
        st = k[1]          # 0 = read 가닥 · 1 = 반대 가닥(그 자리의 read 염기는 G 다)
        for pos, ql in v:
            p = call.get(pos, [0.0, 0.0, st])
            p[2] = st
            if code == "m":
                p[0] = max(p[0], ql / 255.0)
            else:
                p[1] = max(p[1], ql / 255.0)
            call[pos] = p
    if not call:
        continue
    MOD[rd.query_name] = call
    SEQ[rd.query_name] = q
    for off in range(0, len(q) - WIN + 1, WIN):
        META.append((rd.query_name, off))
        wf.write(">w%d%s%s%s" % (len(META) - 1, NL, q[off:off + WIN], NL))
wf.close()
bam.close()
best = {}
rm = subprocess.run("minimap2 -cx map-ont -N 5 -p 0.6 -t 8 %s %s 2>/dev/null"
                    % (cp, os.path.join(TMP, "win.fa")), shell=True, capture_output=True, text=True)
for L in rm.stdout.splitlines():
    f = L.split(TAB)
    if len(f) < 12:
        continue
    nm, tgt = f[0], f[5]
    nmatch, blk = int(f[9]), int(f[10])
    sc = 0
    for x in f[12:]:
        if x.startswith("AS:i:"):
            sc = int(x[5:])
    if nm not in best or sc > best[nm][0]:
        best[nm] = (sc, tgt, nmatch / float(blk) if blk else 0, int(f[2]), int(f[3]))

NG = [0, 0]
for rid, call in MOD.items():
    sq = SEQ.get(rid, "")
    for p, v in call.items():
        NG[0] += 1
        if v[2] == 1:      # 반대 가닥 : read 에서는 G 이고 앞 염기가 C 면 CpG
            if p - 1 >= 0 and sq[p - 1] == "C":
                NG[1] += 1
        elif p + 1 < len(sq) and sq[p + 1] == "G":
            NG[1] += 1
CPGONLY = (NG[1] / float(NG[0])) >= 0.75 if NG[0] else False


COMP2 = {"A": "T", "C": "G", "G": "C", "T": "A"}


def ctx_of(seq, p, st):
    if CPGONLY:
        return "CG"
    # 배치 없이 문맥 판정 : read 서열만 본다. 반대 가닥이면 거꾸로 읽는다
    if st == 1:
        a = COMP2.get(seq[p - 1], "N") if p - 1 >= 0 else "N"
        b = COMP2.get(seq[p - 2], "N") if p - 2 >= 0 else "N"
    else:
        a = seq[p + 1] if p + 1 < len(seq) else "N"
        b = seq[p + 2] if p + 2 < len(seq) else "N"
    if a == "G":
        return "CG"
    if b == "G":
        return "CHG"
    return "CHH"


HW = ["species_key", "sample", "read_id", "win_offset", "ref_type", "context", "identity", "aln_frac",
      "n_c", "n_5mC", "frac_5mC", "n_5hmC"]
ow = io.open(os.path.join(OUT, "%s__%s.noplace_window.tsv" % (SP, SAMPLE)), "w",
             encoding="utf8", newline=NL)
ow.write(TAB.join(HW) + NL)
byread = collections.defaultdict(lambda: collections.defaultdict(lambda: [0, 0, 0]))  # read -> (종류,문맥)
tot = collections.Counter()
for i, (rid, off) in enumerate(META):
    b = best.get("w%d" % i)
    if not b or b[2] < MINID:
        tot["unassigned"] += 1
        continue
    sc, tgt, idt, qs, qe = b
    call = MOD.get(rid, {})
    sq = SEQ.get(rid, "")
    acc = collections.defaultdict(lambda: [0, 0, 0])
    for p, (pm, ph, st) in call.items():
        if off <= p < off + WIN:
            c = acc[ctx_of(sq, p, st)]
            if max(pm, ph) < MODTHR and max(pm, ph) > 1.0 - MODTHR:
                continue          # 애매한 호출은 버린다(modkit 의 저신뢰 제거에 해당)
            c[0] += 1
            if pm >= MODTHR:
                c[1] += 1
            if ph >= MODTHR:
                c[2] += 1
    if not acc:
        continue
    for c, (n, m, h) in sorted(acc.items()):
        if n == 0:          # 애매한 호출만 있던 문맥은 건너뛴다(문턱을 올리면 생긴다)
            continue
        ow.write(TAB.join([SP, SAMPLE, rid, str(off), tgt, c, "%.3f" % idt,
                           "%.3f" % ((qe - qs) / float(WIN)), str(n), str(m),
                           "%.3f" % (m / float(n)), str(h)]) + NL)
        a = byread[rid][(tgt, c)]
        a[0] += n
        a[1] += m
        a[2] += h
    tot[tgt] += 1
ow.close()
HR = ["species_key", "sample", "read_id", "ref_type", "context", "n_c", "n_5mC", "frac_5mC", "n_5hmC"]
orr = io.open(os.path.join(OUT, "%s__%s.noplace_read.tsv" % (SP, SAMPLE)), "w",
              encoding="utf8", newline=NL)
orr.write(TAB.join(HR) + NL)
cnt = collections.Counter()
agg = collections.defaultdict(lambda: [0, 0, 0, 0])
for rid, d in byread.items():
    for (tgt, c), (n, m, h) in d.items():
        orr.write(TAB.join([SP, SAMPLE, rid, tgt, c, str(n), str(m),
                            "%.3f" % (m / float(n)) if n else "-", str(h)]) + NL)
        a = agg[(tgt, c)]
        a[0] += 1
        a[1] += n
        a[2] += m
        a[3] += h
orr.close()
for f in os.listdir(TMP):
    os.remove(os.path.join(TMP, f))
os.rmdir(TMP)
print("%-22s %-14s 문턱 %.2f · 모델 %s(MM 자리 중 다음염기 G %.3f) · read %5d · 창 %7d · 미배정 %6d"
      % (SP, SAMPLE, MODTHR, "CpG전용" if CPGONLY else "전문맥",
         NG[1] / float(NG[0]) if NG[0] else 0, nread, len(META), tot["unassigned"]))
for t in sorted(agg):
    a = agg[t]
    print("   %-46s %-4s read %5d · C %9d · 5mC %.3f · 5hmC %.4f"
          % (t[0][:46], t[1], a[0], a[1], a[2] / float(a[1]) if a[1] else 0,
             a[3] / float(a[1]) if a[1] else 0))
