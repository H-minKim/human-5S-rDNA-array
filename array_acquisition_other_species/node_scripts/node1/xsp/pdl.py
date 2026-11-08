# -*- coding: utf-8 -*-
# 병렬 · 순서 보존 HTTP 범위 내려받기 (EBI 가 연결당 ~0.7 MB/s 라서 여러 연결을 동시에 쓴다)
#   python3 pdl.py URL OUT [연결수=16] [조각MB=32]
#   OUT = 파일 경로 -> 이어받기 가능 (이미 있는 크기가 조각 배수이면 거기서부터)
#   OUT = -        -> 표준출력으로 순서대로 흘림 (tar 에 바로 물림 · 디스크에 전체를 쓰지 않는다)
#   조각마다 받은 길이 = 요청 길이 확인, 실패 시 재시도 (최대 20 회, 지수 대기)
#   진행 : 1 GB 마다 표준오류에 한 줄
# MIT License · Copyright (c) 2026 Hungmin Kim
import os
import sys
import time
import queue
import threading
import urllib.request

URL = sys.argv[1]
OUT = sys.argv[2]
N = int(sys.argv[3]) if len(sys.argv) > 3 else 16
CH = int(sys.argv[4] if len(sys.argv) > 4 else 32) * 1024 * 1024
AHEAD = 2 * N


def head_size():
    for k in range(10):
        try:
            r = urllib.request.urlopen(urllib.request.Request(URL, method="HEAD"), timeout=60)
            return int(r.headers["Content-Length"])
        except Exception:
            time.sleep(2 ** k)
    raise RuntimeError("HEAD 실패")


SIZE = head_size()
start = 0
if OUT != "-" and os.path.exists(OUT):
    have = os.path.getsize(OUT)
    start = (have // CH) * CH
    with open(OUT, "r+b") as g:
        g.truncate(start)
CHUNKS = [(i, a, min(a + CH, SIZE) - 1) for i, a in enumerate(range(start, SIZE, CH))]
RES = {}
CV = threading.Condition()
NXT = [0]
WQ = queue.Queue()
for c in CHUNKS:
    WQ.put(c)
ERR = []


def fetch(a, b):
    for k in range(20):
        try:
            rq = urllib.request.Request(URL, headers={"Range": "bytes=%d-%d" % (a, b)})
            d = urllib.request.urlopen(rq, timeout=180).read()
            if len(d) == b - a + 1:
                return d
        except Exception:
            pass
        time.sleep(min(60, 2 ** k))
    raise RuntimeError("조각 실패 %d-%d" % (a, b))


def worker():
    while not ERR:
        try:
            i, a, b = WQ.get_nowait()
        except queue.Empty:
            return
        with CV:
            while i - NXT[0] > AHEAD and not ERR:
                CV.wait(5)
        try:
            d = fetch(a, b)
        except Exception as e:
            ERR.append(str(e))
            with CV:
                CV.notify_all()
            return
        with CV:
            RES[i] = d
            CV.notify_all()


TH = [threading.Thread(target=worker, daemon=True) for _ in range(N)]
for t in TH:
    t.start()
g = sys.stdout.buffer if OUT == "-" else open(OUT, "ab")
t0 = time.time()
done = start
mark = done // 1000000000
for i in range(len(CHUNKS)):
    with CV:
        while i not in RES and not ERR:
            CV.wait(5)
        if ERR:
            sys.stderr.write("pdl 실패 : %s\n" % ERR[0])
            sys.exit(2)
        d = RES.pop(i)
        NXT[0] = i + 1
        CV.notify_all()
    g.write(d)
    done += len(d)
    if done // 1000000000 > mark:
        mark = done // 1000000000
        el = time.time() - t0
        sys.stderr.write("pdl %d/%d GB · %.1f MB/s\n" % (done // 1000000000, SIZE // 1000000000, (done - start) / 1e6 / el))
        sys.stderr.flush()
g.flush()
if OUT != "-":
    g.close()
    if os.path.getsize(OUT) != SIZE:
        sys.stderr.write("pdl 크기 불일치 %d != %d\n" % (os.path.getsize(OUT), SIZE))
        sys.exit(3)
sys.stderr.write("pdl 완료 %d bytes · %.0f s\n" % (SIZE, time.time() - t0))
