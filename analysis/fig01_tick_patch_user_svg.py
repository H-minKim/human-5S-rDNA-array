# -*- coding: utf-8 -*-
# 사용자 정본 SVG 에 빠진 유닛 눈금을 넣는다 (2026-10-19)
#
# 무엇이 빠졌나 : 왼쪽 보라 read 2·3 번째 줄의 눈금이 x 86.4 에서 멈췄다.
#   막대는 104.9 · 119.5 까지 가므로 18.5 pt · 33.1 pt 가 민둥이었다.
# 까닭(제 코드) : 눈금을 그린 뒤에 라벨 덮어쓰기가 `graphics` 인자 없이 돌아
#   닿은 밑그림(= 그 눈금)을 같이 지웠다. 로그는 21·26 개를 그렸다고 적는데
#   결과물엔 17 개뿐이었다 — 그린 뒤에 지워진 것이다.
# 좌표는 고친 PDF 에서 그대로 가져왔다(짐작이 아니다). 간격 3.08 pt.
# 사용자 판의 꼴을 그대로 따른다 : `<path class="st53" d="M{x},{y}v2.9"/>` · 소수 1 자리.
import io, os, re, sys, hashlib
sys.stdout.reconfigure(encoding="utf-8")
# -- (repo) paths come from the environment --------------------------------
#   저장소판에서 더한 줄 (2026-10-21). 전에는 돌리는 폴더의 상대 이름이었다 — 뜻은 같다.
SRC = os.environ.get("USER_SVG", "USER_Fig1.svg")
OUT = os.environ.get("PATCHED_SVG", "Fig1_tickfix.svg")

NEW = {
    "229.5": [89.5, 92.6, 95.7, 98.7, 101.8],                       # 막대 끝 104.9
    "231.8": [89.5, 92.6, 95.7, 98.7, 101.8, 104.9, 108.0, 111.1, 114.2, 117.2],  # 막대 끝 119.5
}

s0 = io.open(SRC, encoding="utf8", errors="replace").read()
s = s0
added = 0
for y, xs in NEW.items():
    anchor = '<path class="st53" d="M86.4,%sv2.9"/>' % y
    assert s.count(anchor) == 1, (y, s.count(anchor))
    ins = anchor + "".join('<path class="st53" d="M%s,%sv2.9"/>' % (("%.1f" % x), y) for x in xs)
    s = s.replace(anchor, ins)
    added += len(xs)
    print("   y %s : %d 개 넣음 → x %s" % (y, len(xs), xs))

assert len(s) > len(s0)
io.open(OUT, "w", encoding="utf8", newline="\n").write(s)
print()
print("넣은 눈금 : %d 개" % added)
print("%s  %s -> %s" % (OUT, hashlib.md5(s0.encode("utf8")).hexdigest()[:12],
                         hashlib.md5(s.encode("utf8")).hexdigest()[:12]))

# 검산 — 줄마다 눈금이 막대 끝까지 가나
import collections
BAR = {"229.5": 104.9, "231.8": 119.5, "227.0": 87.3, "234.2": 133.2, "236.6": 123.4}
rows = collections.defaultdict(list)
for m in re.finditer(r'<path class="st53" d="M([0-9.]+),([0-9.]+)v[-0-9.]+"/>', s):
    rows[m.group(2)].append(float(m.group(1)))
print()
print("검산 : 줄마다 마지막 눈금 → 막대 끝")
for y in sorted(BAR, key=float):
    xs = sorted(rows.get(y, []))
    gap = BAR[y] - xs[-1] if xs else 999
    print("   y %-6s 눈금 %2d  마지막 %7.2f  막대끝 %7.2f  빈 구간 %5.2f pt  %s"
          % (y, len(xs), xs[-1], BAR[y], gap, "민둥" if gap > 5 else "ok"))
