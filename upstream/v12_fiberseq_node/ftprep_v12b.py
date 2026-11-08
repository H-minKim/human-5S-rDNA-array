# -*- coding: utf-8 -*-
# m6A called afresh from PacBio kinetics and stored separately
#   대상 : 원자료에 m6A 호출이 없는 개체 (HG03672 · HG02624 = Y5 제출분 jasmine 2.2.0, 5mC 만)
#   다른 개체의 m6A 는 기기 jasmine 2.7.x -> 호출기가 다르므로 본 분석과 합치지 않는다 (출처 PACBIO_FIBERSEQ_FT)
#   순서
#     1. 배치 SAM(~/aog8v12pb/sam/<sample>.AOJV_PB.sam) 의 주 정렬 기록을 원본 read 방향 비정렬 기록으로 되돌림
#        (역정렬이면 SEQ 역상보·QUAL 뒤집기; kinetics fi/fp/ri/rp/fn/rn 과 MM/ML 은 원본 방향 그대로라 손대지 않음)
#        @RG = 원자료 BAM 머리줄(movie 별) · 각 read 에 RG:Z 붙임
#     2. ft predict-m6a (fibertools-rs v0.13.0 · 기본값)
#     3. 새 MM/ML 을 원래 정렬 기록에 넣은 새 SAM (~/pbft_v12/sam) — 다른 열·태그는 그대로
#     검정 : read 마다 ft 전후 C+m 블록(자리·값) 동일 · SEQ 동일
#   인자 : 개체
import io
import os
import sys
import glob
import subprocess
import collections

NL = chr(10)
TAB = chr(9)
H = os.path.expanduser
S = sys.argv[1]
SRC = H(os.environ.get("SRCD", "~/aog8v12pb/sam") + "/%s.AOJV_PB.sam" % S)   # v12b : SAM 폴더·작업폴더 환경변수 (기본값 = 옛 동작)
W = H(os.environ.get("FTW", "~/pbft_v12"))
FT = H("~/bin/ft")
os.makedirs(os.path.join(W, "sam"), exist_ok=True)
os.makedirs(os.path.join(W, "ft"), exist_ok=True)
CMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
KEEP = ("MM:Z:", "ML:B:", "fi:B:", "fp:B:", "ri:B:", "rp:B:", "fn:i:", "rn:i:")


def rc(s):
    return "".join(CMP.get(c, "N") for c in reversed(s))


def blocks(mm, ml):
    """{블록머리: (deltas 문자열, 값 목록)}"""
    vals = ml.split(",")[1:] if ml else []
    out = {}
    base = 0
    for b in [p for p in mm.split(";") if p]:
        f = b.split(",")
        n = len(f) - 1
        out[f[0]] = (",".join(f[1:]), vals[base:base + n])
        base += n
    return out


def tag(fields, pre):
    for x in fields:
        if x.startswith(pre):
            return x[len(pre):]
    return None


# ── 1. 비정렬 BAM
RG = {}
hdr = ["@HD" + TAB + "VN:1.6" + TAB + "SO:unknown"]
for fp in sorted(glob.glob(os.path.join(W, "rg", "*.hdr"))):
    mv = os.path.basename(fp)[:-4]
    for L in io.open(fp, encoding="utf8"):
        if L.startswith("@RG"):
            hdr.append(L.rstrip(NL))
            RG[mv] = [x for x in L.rstrip(NL).split(TAB) if x.startswith("ID:")][0][3:]
recs = []
orig = {}
for L in io.open(SRC, encoding="utf8"):
    if L.startswith("@"):
        continue
    r = L.rstrip(NL).split(TAB)
    if int(r[1]) & 0x900:
        continue
    mv = r[0].split("/")[0]
    seq, qual = r[9], r[10]
    if int(r[1]) & 16:
        seq, qual = rc(seq), qual[::-1]
    tg = [x for x in r[11:] if x.startswith(KEEP)]
    recs.append(TAB.join([r[0], "4", "*", "0", "0", "*", "*", "0", "0", seq, qual] + tg + ["RG:Z:" + RG[mv]]))
    orig[r[0]] = (seq, tag(r[11:], "MM:Z:"), tag(r[11:], "ML:B:"))
ub = os.path.join(W, "ft", "%s.unaligned.bam" % S)
sp = subprocess.run(["samtools", "view", "-b", "-o", ub, "-"], input=NL.join(hdr + recs) + NL, text=True, capture_output=True)
if sp.returncode:
    sys.exit("samtools 실패 " + sp.stderr)

# ── 2. ft predict-m6a
ob = os.path.join(W, "ft", "%s.ft.bam" % S)
p = subprocess.run([FT, "predict-m6a", "-t", "8", ub, ob], capture_output=True, text=True)
io.open(os.path.join(W, "ft", "%s.ft.log" % S), "w", encoding="utf8").write(p.stdout + p.stderr)
if p.returncode:
    sys.exit("ft 실패 " + p.stderr[-500:])

# ── 3. 새 MM/ML + 검정
T = collections.Counter()
NEW = {}
heads = collections.Counter()
pr = subprocess.run(["samtools", "view", ob], capture_output=True, text=True).stdout.splitlines()
for L in pr:
    r = L.split(TAB)
    mm, ml = tag(r[11:], "MM:Z:"), tag(r[11:], "ML:B:")
    T["ft 기록"] += 1
    o = orig.get(r[0])
    if o is None:
        T["원본에 없는 read"] += 1
        continue
    if r[9] != o[0]:
        T["SEQ 다름"] += 1
    a, b = blocks(o[1] or "", o[2] or ""), blocks(mm or "", ml or "")
    ca = [(k.rstrip(".?"), v) for k, v in a.items() if k.startswith("C+m")]
    cb = [(k.rstrip(".?"), v) for k, v in b.items() if k.startswith("C+m")]
    T["C+m 같음" if ca == cb else "C+m 다름"] += 1
    for k in b:
        heads[k] += 1
    if any(k.startswith(("A+a", "T-a")) and v[0] for k, v in b.items()):
        T["m6A 호출 있음"] += 1
    NEW[r[0]] = (mm, ml)
T["원본 주정렬 read"] = len(orig)

out = io.open(os.path.join(W, "sam", "%s.AOJV_PB.sam" % S), "w", encoding="utf8", newline=NL)
for L in io.open(SRC, encoding="utf8"):
    if L.startswith("@"):
        out.write(L)
        continue
    r = L.rstrip(NL).split(TAB)
    if not (int(r[1]) & 0x900) and r[0] in NEW and NEW[r[0]][0] is not None:
        mm, ml = NEW[r[0]]
        r = [x for x in r if not x.startswith(("MM:Z:", "ML:B:"))] + ["MM:Z:" + mm, "ML:B:" + ml]
        T["새 MM/ML 넣은 기록"] += 1
    out.write(TAB.join(r) + NL)
out.close()
for k, v in sorted(T.items()):
    print("  %s %-20s %d" % (S, k, v))
print("  %s MM 블록머리 %s" % (S, dict(heads)))
