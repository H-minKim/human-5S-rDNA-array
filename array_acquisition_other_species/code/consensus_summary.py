# -*- coding: utf-8 -*-
# 07 합본 요약 + 체계적 불일치 진단
#
#   read 합의가 조립본과 다를 때 두 가지를 갈라야 한다.
#     · 흩어진 불일치 = 잡음 (자리마다 몇 유닛에서만 불린다)
#     · 체계적 불일치 = 조립본(또는 참조 유닛)이 그 자리에서 틀렸거나 read 가 다른 계통이다
#       → 같은 자리가 거의 모든 유닛에서 불린다
#   지표 : read 전용 변이 자리마다 "몇 유닛에서 불렸나 / 전체 유닛" 의 중앙값 = systematic_frac
#          0.5 이상이면 체계적으로 본다 (실측 : 시아망 0.9 · 고릴라 0.08)
#   인자 : 07 폴더 출력파일
# MIT License · Copyright (c) 2026 Hungmin Kim
import io
import os
import sys
import csv
import glob
import collections

TAB = chr(9)
NL = chr(10)
D, OUTF = sys.argv[1], sys.argv[2]
HDR = ["species_key", "sample", "n_units", "depth_median", "n_var_asm", "n_var_cons", "n_var_strict",
       "n_shared", "n_read_only", "n_asm_only", "concordance", "concordance_strict",
       "read_only_sites", "systematic_frac", "systematic_sites", "verdict"]
o = io.open(OUTF, "w", encoding="utf8", newline=NL)
o.write(TAB.join(HDR) + NL)
print("%-22s %-14s %6s %6s %8s %8s %7s %7s %6s %s"
      % ("종", "샘플", "유닛", "깊이", "조립본", "합의", "일치", "엄격", "체계", "판정"))
for f in sorted(glob.glob(os.path.join(D, "*.unit_consensus.tsv"))):
    b = os.path.basename(f)[:-len(".unit_consensus.tsv")]
    sp, sam = b.split("__", 1)
    R = list(csv.DictReader(io.open(f, encoding="utf8"), delimiter=TAB))
    if not R:
        continue
    n = len(R)
    dep = sorted(int(r["depth_median"]) for r in R)
    A = sum(int(r["n_var_asm"]) for r in R)
    C = sum(int(r["n_var_cons"]) for r in R)
    ST = sum(int(r["n_var_strict"]) for r in R)
    SH = sum(int(r["n_shared"]) for r in R)
    RO = sum(int(r["n_read_only"]) for r in R)
    AO = sum(int(r["n_asm_only"]) for r in R)
    SHs = sum(int(r["n_shared_strict"]) for r in R)
    ASs = sum(int(r["n_asm_strict"]) for r in R)
    conc = 2.0 * SH / (2.0 * SH + RO + AO) if (SH + RO + AO) else 0
    dcs = ST + ASs
    concs = 2.0 * SHs / dcs if dcs else 0
    sf = os.path.join(D, b + ".site_support.tsv")
    per = collections.defaultdict(set)
    units = set()
    if os.path.exists(sf):
        for r in csv.DictReader(io.open(sf, encoding="utf8"), delimiter=TAB):
            units.add((r["array_id"], r["array_pos"]))
            if r["called"] == "Y" and r["in_assembly"] == "N":
                per[(r["array_id"], r["pos_in_refunit"])].add(r["array_pos"])
    # array 마다 유닛 수
    uc = collections.Counter(a for a, _p in units)
    fr = sorted((len(v) / float(uc[k[0]])) for k, v in per.items() if uc.get(k[0]))
    med = fr[len(fr) // 2] if fr else 0.0
    nsys = sum(1 for x in fr if x >= 0.5)
    v = ("체계적 불일치(조립본/참조 재검토)" if med >= 0.5 else
         ("일부 체계적" if nsys >= 5 else "흩어짐(잡음)"))
    o.write(TAB.join([sp, sam, str(n), str(dep[len(dep) // 2]), str(A), str(C), str(ST), str(SH),
                      str(RO), str(AO), "%.3f" % conc, "%.3f" % concs, str(len(fr)),
                      "%.3f" % med, str(nsys), v]) + NL)
    print("%-22s %-14s %6d %6d %8d %8d %7.3f %7.3f %6.2f %s"
          % (sp[:22], sam[:14], n, dep[len(dep) // 2], A, C, conc, concs, med, v))
o.close()
