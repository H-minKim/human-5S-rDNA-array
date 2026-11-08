# -*- coding: utf-8 -*-
# 07 합본 — read 를 유닛마다 모아 자리별 합의 변이를 부른다 (타 생물종 전용)
#
#   왜 필요한가 : read 한 개의 변이 목록은 대부분 ONT 오류다(고릴라 실측 중앙 94 개/2,248 bp = 4.2%).
#     유닛마다 그 유닛을 덮은 read 들을 모아 자리마다 깊이로 나눠야 진짜 변이가 남는다.
#
#   방법 (사람 파이프라인의 깊이 문턱 합의와 같은 규칙)
#     1) 유닛 = (array_id, array_pos). 그 유닛을 덮은 read 조각을 모은다(NOALN 제외)
#     2) 자리마다 덮은 깊이 = 그 자리를 cov_start_in_ref~cov_end_in_ref 안에 넣은 조각 수
#     3) 변이 토큰마다 지지 수 / 그 자리 깊이 = 지지율. 지지율 >= THR 이고 깊이 >= MINDEP 이면 채택
#     4) 채택된 토큰을 자리순으로 묶어 유닛의 합의 변이 목록으로 삼는다
#   채움 금지 — 안 덮인 자리는 그냥 비운다(cons_cov_start/end 로 어디까지 덮었는지 적는다)
#   mapq 로 거르지 않는다(탠덤 array 안은 mapq0 이 정상)
#
#   엄격판(varlist_strict) : 치환 + 3 bp 초과 indel 만. ONT 는 호모폴리머 결실을 흔히 내서
#     1~3 bp 결실이 문턱을 넘어버린다(고릴라 실측 : 채택 결실 204 개 중 121 개가 조립본에 없다.
#     치환은 1,514/1,593 = 95% 가 조립본과 일치). 유닛 단위 결론은 엄격판으로 낼 것.
#   검산 : 같은 유닛의 조립본 변이(g8 지도)와 얼마나 맞는지 세어 함께 낸다(전체판·엄격판 따로)
#   감도 : THR 0.20 / 0.35 / 0.50 세 값의 변이 수를 같이 적는다
#
#   산출 : <out>/<species>__<sample>.unit_consensus.tsv   (유닛 한 줄)
#          <out>/<species>__<sample>.site_support.tsv      (지지율 >= REPORT 인 변이 한 줄)
#   인자 : read유닛표 참조유닛폴더 지도폴더 출력폴더 [THR=0.35] [MINDEP=3] [REPORT=0.10]
# MIT License · Copyright (c) 2026 Hungmin Kim
import io
import os
import re
import sys
import csv
import glob
import collections

NL = chr(10)
TAB = chr(9)
RUF, RDIR, MAPD, OUT = sys.argv[1:5]
THR = float(sys.argv[5]) if len(sys.argv) > 5 else 0.35
MINDEP = int(sys.argv[6]) if len(sys.argv) > 6 else 3
REPORT = float(sys.argv[7]) if len(sys.argv) > 7 else 0.10
os.makedirs(OUT, exist_ok=True)
NUM = re.compile(r"^([0-9]+)")
pos_of = lambda v: int(NUM.match(v).group(1)) if NUM.match(v) else -1
base = os.path.basename(RUF)[:-len(".read_unit_seq.tsv")]
SP, SAMPLE = base.split("__", 1)


def rdlen(p):
    return len("".join(x.strip() for x in io.open(p, encoding="utf8").read().split(NL)[1:]))


ULEN, ASM = {}, {}
for f in glob.glob(os.path.join(RDIR, "*.unit.fa")):
    ULEN[os.path.basename(f)[:-len(".unit.fa")]] = rdlen(f)
if os.path.isdir(MAPD):
    for f in glob.glob(os.path.join(MAPD, "*.g8map.tsv")):
        for r in csv.DictReader(io.open(f, encoding="utf8"), delimiter=TAB):
            ASM[(r["array_id"], r["array_pos"])] = (
                set() if r["varlist"] in ("-", "NOALN", "REV", "") else set(r["varlist"].split(",")),
                r["verify"])

HU = ["species_key", "sample", "array_id", "array_pos", "refunit_len", "n_reads", "n_reads_aln",
      "n_reads_span", "depth_median", "cons_cov_start", "cons_cov_end", "n_var_cons",
      "n_sub", "n_ins", "n_del", "varlist_cons",
      "n_var_strict", "varlist_strict", "n_var_thr20", "n_var_thr50", "asm_verify", "n_var_asm",
      "n_shared", "n_read_only", "n_asm_only", "concordance",
      "n_asm_strict", "n_shared_strict", "concordance_strict"]
HS = ["species_key", "sample", "array_id", "array_pos", "var", "kind", "pos_in_refunit",
      "support", "depth", "support_frac", "called", "in_assembly"]
ou = io.open(os.path.join(OUT, "%s.unit_consensus.tsv" % base), "w", encoding="utf8", newline=NL)
ou.write(TAB.join(HU) + NL)
os_ = io.open(os.path.join(OUT, "%s.site_support.tsv" % base), "w", encoding="utf8", newline=NL)
os_.write(TAB.join(HS) + NL)
tot = collections.Counter()


# 입력은 read 순서라 유닛별로 모여 있지 않다. 한 번 훑으며 유닛마다 필요한 것만 쌓는다
#   (행 전체를 들고 있으면 수 GB 가 된다)
COV = collections.defaultdict(list)      # (array_id, array_pos) -> [(덮은 시작, 끝), ...]
VAR = collections.defaultdict(collections.Counter)
NUMS = collections.defaultdict(lambda: [0, 0, 0])   # 전체 · 정렬됨 · 완전덮음


def flush(key):
    aid, apos = key
    U = ULEN.get(aid)
    if not U:
        return
    iv = COV[key]
    if not iv:
        return
    nall, naln, nspan = NUMS[key]
    dif = [0] * (U + 2)
    cnt = VAR[key]
    for a, b in iv:
        dif[a] += 1
        dif[b + 1] -= 1
    dep = [0] * (U + 2)
    run = 0
    for i in range(1, U + 1):
        run += dif[i]
        dep[i] = run
    covered = [i for i in range(1, U + 1) if dep[i] >= MINDEP]
    c0, c1 = (covered[0], covered[-1]) if covered else (0, 0)
    dm = sorted(dep[1:U + 1])
    asm, ver = ASM.get(key, (set(), "-"))
    call, n20, n50 = [], 0, 0
    for v, c in cnt.items():
        p = pos_of(v)
        if p < 1 or p > U:
            continue
        d = dep[p]
        if d < MINDEP:
            continue
        fr = c / float(d)
        if fr >= 0.20:
            n20 += 1
        if fr >= 0.50:
            n50 += 1
        ok = fr >= THR
        if ok:
            call.append((p, v))
        if fr >= REPORT:
            kind = "del" if "D" in v[len(str(p)):] else ("ins" if "I" in v[len(str(p)):] else "sub")
            os_.write(TAB.join([SP, SAMPLE, aid, apos, v, kind, str(p), str(c), str(d),
                                "%.3f" % fr, "Y" if ok else "N",
                                "Y" if v in asm else "N"]) + NL)
    call.sort()
    cv = [v for _p, v in call]
    cs = set(cv)

    def strict(vs):
        out = []
        for v in vs:
            t = v[len(str(pos_of(v))):]
            if "D" in t or "I" in t:
                if int(t[1:]) > 3:
                    out.append(v)
            else:
                out.append(v)
        return out

    cst = strict(cv)
    ast = set(strict(sorted(asm, key=pos_of)))
    shs = len(set(cst) & ast)
    dcs = len(cst) + len(ast)
    concs = "%.3f" % (2.0 * shs / dcs) if dcs else "-"
    sh = len(cs & asm)
    conc = "%.3f" % (2.0 * sh / (len(cs) + len(asm))) if (cs or asm) else "-"
    nsub = sum(1 for v in cv if "D" not in v[len(str(pos_of(v))):] and "I" not in v[len(str(pos_of(v))):])
    nins = sum(1 for v in cv if "I" in v[len(str(pos_of(v))):])
    ndel = sum(1 for v in cv if "D" in v[len(str(pos_of(v))):])
    ou.write(TAB.join([SP, SAMPLE, aid, apos, str(U), str(nall), str(naln),
                       str(nspan),
                       str(dm[len(dm) // 2]), str(c0), str(c1), str(len(cv)),
                       str(nsub), str(nins), str(ndel), ",".join(cv) if cv else "-",
                       str(len(cst)), ",".join(cst) if cst else "-",
                       str(n20), str(n50), ver, str(len(asm)), str(sh),
                       str(len(cs - asm)), str(len(asm - cs)), conc,
                       str(len(ast)), str(shs), concs]) + NL)
    tot["unit"] += 1
    tot["var"] += len(cv)
    tot["shared"] += sh
    tot["ronly"] += len(cs - asm)
    tot["aonly"] += len(asm - cs)
    tot["svar"] += len(cst)
    tot["sshared"] += shs
    tot["sronly"] += len(set(cst) - ast)
    tot["saonly"] += len(ast - set(cst))


for r in csv.DictReader(io.open(RUF, encoding="utf8"), delimiter=TAB):
    k = (r["array_id"], r["array_pos"])
    n = NUMS[k]
    n[0] += 1
    if r["spans_unit"] == "Y":
        n[2] += 1
    if r["aln_to_refunit"] == "NOALN" or not r["cov_start_in_ref"].isdigit():
        continue
    U = ULEN.get(r["array_id"])
    if not U:
        continue
    a, b = max(1, int(r["cov_start_in_ref"])), min(U, int(r["cov_end_in_ref"]))
    if b < a:
        continue
    n[1] += 1
    COV[k].append((a, b))
    if r["varlist"] not in ("-", "NOALN", ""):
        VAR[k].update(r["varlist"].split(","))
for k in sorted(COV, key=lambda x: (x[0], int(x[1]))):
    flush(k)
ou.close()
os_.close()
d = tot["shared"] + tot["ronly"] + tot["aonly"]
ds = 2.0 * tot["sshared"] + tot["sronly"] + tot["saonly"]
print("%-22s %-14s 유닛 %6d · 합의 %7d (일치 %.3f) · 엄격 %7d (일치 %.3f)"
      % (SP, SAMPLE, tot["unit"], tot["var"],
         (2.0 * tot["shared"] / (2.0 * tot["shared"] + tot["ronly"] + tot["aonly"])) if d else 0,
         tot["svar"], (2.0 * tot["sshared"] / ds) if ds else 0))
