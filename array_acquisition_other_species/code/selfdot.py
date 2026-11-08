# -*- coding: utf-8 -*-
# 자기-도트플롯으로 주기를 직접 본다 (자·주석·내 유닛 모두 안 씀)
#   array 조각을 자기 자신에 blastn -> HSP 대각선 옵셋(sstart-qstart) 을 정렬길이로 가중해 모은다
#   진짜 탠덤이면 옵셋 히스토그램이 주기의 배수마다 서는 사다리가 된다
#   대조군으로 사람·금화조도 같이 돌려 눈금을 본다
#   산출 : SELFDOT.tsv (array 마다 상위 옵셋 8 개와 그 정렬 총량)
# MIT License · Copyright (c) 2026 Hungmin Kim
import io
import os
import sys
import glob
import subprocess
import collections

NL = chr(10)
TAB = chr(9)
H = os.path.expanduser("~/lit5s")
PICK = sys.argv[1] if len(sys.argv) > 1 else "rhesus,human_CHM13v2,orang_bornean,zfinch_v03_dip__CM109753,zebrafish_GRCz12tu__CM111185,zebrafish_GRCz12tu__CM111196,mouse"
TAGS = PICK.split(",")
BIN = 20          # 옵셋 바구니
MINALN = 200      # 이 길이 이상 HSP 만 (짧은 조각 반복 배제)
w = io.open(os.path.join(H, "SELFDOT.tsv"), "w", encoding="utf8", newline=NL)
w.write(TAB.join(["array_id", "slice_len", "n_hsp", "aln_total_bp", "top_offsets(bp:정렬총량)"]) + NL)
for sfa in sorted(glob.glob(os.path.join(H, "slices", "*.arr.fa"))):
    aid = os.path.basename(sfa).replace(".arr.fa", "")
    if not any(t in aid for t in TAGS):
        continue
    hd = io.open(sfa, encoding="utf8").readline()
    kv = dict(x.split("=", 1) for x in hd.split() if "=" in x)
    slen = int(kv["slice_len"])
    r = subprocess.run(["blastn", "-task", "blastn", "-word_size", "16", "-dust", "no", "-evalue", "1e-20",
                        "-query", sfa, "-subject", sfa, "-outfmt", "6 qstart qend sstart send pident length"],
                       capture_output=True, text=True)
    off = collections.Counter()
    n, tot = 0, 0
    for L in r.stdout.splitlines():
        f = L.split(TAB)
        qs, qe, ss, se, pid, aln = int(f[0]), int(f[1]), int(f[2]), int(f[3]), float(f[4]), int(f[5])
        if se < ss:
            continue                      # 역방향은 제외 (주기 아님)
        d = ss - qs
        if d <= 0 or aln < MINALN:
            continue
        n += 1
        tot += aln
        off[int(round(d / float(BIN))) * BIN] += aln
    top = off.most_common(8)
    print("%-50s 조각 %6d · HSP %6d · 정렬 %9d" % (aid, slen, n, tot))
    print("      " + " | ".join("%d:%d" % (k, v) for k, v in top))
    w.write(TAB.join([aid, str(slen), str(n), str(tot),
                      ",".join("%d:%d" % (k, v) for k, v in top)]) + NL)
    sys.stdout.flush()
w.close()
print("표 : ~/lit5s/SELFDOT.tsv")
