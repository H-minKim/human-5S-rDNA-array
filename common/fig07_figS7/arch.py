# -*- coding: utf-8 -*-
# Version-archiving helper
#   규칙 : 무엇을 고치든 ① 이전 파일은 절대 지우지 않고 archive\<date>_v<edition>\ 로 옮겨 적고(복사)
#          ② 새 파일은 latest\ 에만 둔다 ③ 대장에 md5·날짜·사유를 남긴다.
#
#   쓰는 법
#     python arch.py put <새파일> <최신안_상대경로> --사유 "..."     새 파일을 최신에 넣는다(있으면 옛 판을 보관으로)
#     python arch.py archive <최신안_상대경로> --사유 "..."          최신에 있는 파일을 보관으로만 넘긴다(대체 없음)
#     python arch.py index                                        대장(최신·보관) 다시 만든다
#   폴더 : latest\{코드,표,그림,로그,문서}  ·  archive\<YYYY-MM-DD>_v<N>_<reason>\  ·  90_ledger\
import io, os, sys, time, shutil, hashlib, glob, argparse
import pandas as pd

import os as _os_need, sys as _sys_need


def _need(k):
    """No defaults are set: if an environment variable is missing the script stops instead of silently reading the wrong location."""
    v = _os_need.environ.get(k)
    if not v:
        _sys_need.exit("environment variable %s is not set (the deposited copy has no default; "
                       "see README.md, section 'How to run it')" % k)
    return v


ROOT = _os_need.path.dirname(_need("RDNA_OUT_DIR"))
NEW = ROOT + "/latest"
ARC = ROOT + "/archive"
LED = ROOT + "/90_ledger"
for d in (NEW, ARC, LED):
    os.makedirs(d, exist_ok=True)
for d in ("코드", "표", "그림", "로그", "문서"):
    os.makedirs(NEW + "/" + d, exist_ok=True)


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def mt(p):
    return time.strftime("%Y-%m-%d %H:%M", time.localtime(os.path.getmtime(p)))


def log_line(row):
    p = LED + "/CHANGE_LOG.tsv"
    new = not os.path.exists(p)
    with io.open(p, "a", encoding="utf8", newline="\n") as f:
        if new:
            f.write("시각\t동작\t최신경로\t보관경로\t원본에서\tmd5\t사유\n")
        f.write("\t".join(str(x) for x in row) + "\n")


def next_ver(rel):
    """이 파일이 몇 번째로 보관되는가"""
    n = 0
    for d in glob.glob(ARC + "/*/"):
        if os.path.exists(os.path.join(d, rel)):
            n += 1
    return n + 1


def archive(rel, reason):
    cur = NEW + "/" + rel
    if not os.path.exists(cur):
        return None
    v = next_ver(rel)
    box = "%s/%s_v%d_%s" % (ARC, time.strftime("%Y-%m-%d"), v, (reason or "수정")[:24].replace("/", "_").replace(" ", "_"))
    dst = box + "/" + rel
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    shutil.copy2(cur, dst)          # 복사 : 최신 쪽은 곧 새 파일로 덮이지만 보관본이 남는다
    log_line([time.strftime("%Y-%m-%d %H:%M"), "보관", cur, dst, "", md5(dst), reason])
    print("보관 : %s → %s (md5 %s · 수정일 %s)" % (rel, dst, md5(dst)[:12], mt(dst)))
    return dst


def put(src, rel, reason):
    cur = NEW + "/" + rel
    os.makedirs(os.path.dirname(cur), exist_ok=True)
    old = None
    if os.path.exists(cur):
        if md5(cur) == md5(src):
            print("건너뜀(내용 같음) : %s" % rel)
            return
        old = archive(rel, reason)
    shutil.copy2(src, cur)
    log_line([time.strftime("%Y-%m-%d %H:%M"), "새판" if old else "새로", cur, old or "", src, md5(cur), reason])
    print("최신 : %s ← %s (md5 %s)" % (rel, src, md5(cur)[:12]))


def index():
    rows = []
    for p in sorted(glob.glob(NEW + "/**/*", recursive=True)):
        if os.path.isfile(p):
            rows.append(dict(구분="최신", 경로=p.replace("\\", "/"), 상대=os.path.relpath(p, NEW).replace("\\", "/"),
                             md5=md5(p), 크기=os.path.getsize(p), 수정일=mt(p)))
    pd.DataFrame(rows).to_csv(LED + "/LATEST_FILES.tsv", sep="\t", index=False)
    rows = []
    for p in sorted(glob.glob(ARC + "/**/*", recursive=True)):
        if os.path.isfile(p):
            box = os.path.relpath(p, ARC).split(os.sep)[0]
            rows.append(dict(구분="보관", 보관함=box, 경로=p.replace("\\", "/"),
                             상대=os.path.relpath(p, ARC + "/" + box).replace("\\", "/"),
                             md5=md5(p), 크기=os.path.getsize(p), 수정일=mt(p)))
    pd.DataFrame(rows).to_csv(LED + "/ARCHIVE_FILES.tsv", sep="\t", index=False)
    print("대장 : 90_ledger/LATEST_FILES.tsv · ARCHIVE_FILES.tsv")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["put", "archive", "index"])
    ap.add_argument("a", nargs="?")
    ap.add_argument("b", nargs="?")
    ap.add_argument("--사유", default="")
    x = ap.parse_args()
    if x.cmd == "put":
        put(x.a, x.b or os.path.basename(x.a), x.사유)
    elif x.cmd == "archive":
        archive(x.a, x.사유)
    else:
        index()
