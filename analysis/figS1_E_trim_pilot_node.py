# -*- coding: utf-8 -*-
# The script that ran the trimming test for Fig. S1 E on the compute node
#
#   이 글은 AWS 노드에서 돌았고 로컬에 사본이 없었다.
#   본문은 한 글자도 바꾸지 않았다 — 이 머릿글만 더했다.
#
#   이 글이 낸 것 : `ctrl.log` · `t45.log` · `goctrl_*.log` · `go45_*.log`
#     → 그 로그를 `analysis/figS1_E_trim_gate_table.py` 가 읽어 `FIG_data_table.tsv` 를 만든다.
#
#   쓰는 도구 — `blastn` 은 쓰지 않는다
#     · `minimap2 -x map-ont -t <n> --secondary=no`  (flank 판정 · 노드 기록 `2.26-r1175`)
#     · `hifiasm --ont -o <smp> -t <n>`              (재조립 · `0.25.0-r726`)
#
#   난수 — 씨앗이 박혀 있다 : `random.Random(20260928)` (135 행)
#     array 만 있는 read 는 양 끝 중 한쪽을 무작위로 고르므로 씨앗이 없으면 재현이 안 된다.
#
#   다시 돌리려면 노드가 필요하다 (read 원본이 거기 있다) — 로컬에서는 안 돌아간다.

# -*- coding: utf-8 -*-
# 자르기 파일럿 — read 를 빼지 않고 길이를 자른다 (깊이 유지 · 장거리 연결만 제거)
#   규칙 (사용자 지정)
#     · 한도 L = THR * 2240 * c_max   (THR 기본 0.45)
#     · L 을 넘는 read 만 자른다
#     · flank 가 붙은 read 는 flank 가 아닌 끝을 자른다 (flank 쪽 L bp 를 남긴다)
#     · array only read 는 양 끝 중 무작위로 한쪽을 골라 자른다 (시드 고정)
#   flank 판정 : CHM13 근위·원위 flank 50 kb 탐침에 minimap2 정렬 → read 좌표 어느 쪽에 붙나
#   사용 : python3 trimpilot.py <목록(개체 c_max)> <THR>
import io, os, sys, gzip, random, subprocess, shutil, tempfile, collections

H = os.path.expanduser("~")
W = os.path.join(H, "trimpilot")
os.makedirs(os.path.join(W, "arr"), exist_ok=True)
os.makedirs(os.path.join(W, "log"), exist_ok=True)
SUM = os.path.join(W, "summary.tsv")
HIFIASM = os.path.join(H, "hifiasm")
HIARR = os.path.join(H, "hiarr.py")
PY = os.path.join(H, "miniconda3", "bin", "python")
FL = os.path.join(H, "s5s_tuc", "ref")
UNIT = 2240
TH = os.environ.get("TH", "8")
LIST = sys.argv[1]
THR = float(sys.argv[2]) if len(sys.argv) > 2 else 0.45

PROBE = os.path.join(W, "flank_probe.fa")
if not os.path.exists(PROBE):
    with io.open(PROBE, "w", newline="\n") as g:
        for nm in ("proximal", "distal"):
            p = os.path.join(FL, "FLANK_%s_50000.fa" % nm)
            seq = "".join(x.strip() for x in io.open(p) if not x.startswith(">"))
            g.write(">%s\n%s\n" % (nm, seq))
    print("탐침 만듦 : %s" % PROBE, flush=True)


def find_reads(s):
    for p in (os.path.join(H, "kgp_ont/out/%s.array.fastq.gz" % s),
              os.path.join(H, "s5s_ont/hits/%s.hits.fastq.gz" % s),
              os.path.join(H, "giab_hits/%s.hits.fastq.gz" % s),
              os.path.join(H, "cancer5s/hits/%s.hits.fastq.gz" % s)):
        if os.path.exists(p):
            return p
    return None


def flank_side(fq, tmp):
    """read 이름 -> 'L'(read 앞쪽에 flank) / 'R' / None"""
    paf = os.path.join(tmp, "fl.paf")
    with io.open(paf, "w") as g:
        subprocess.call(["minimap2", "-x", "map-ont", "-t", TH, "--secondary=no", PROBE, fq],
                        stdout=g, stderr=subprocess.DEVNULL)
    best = {}
    for L in io.open(paf):
        f = L.rstrip("\n").split("\t")
        if len(f) < 12:
            continue
        q, ql, qs, qe = f[0], int(f[1]), int(f[2]), int(f[3])
        bl, mq = int(f[10]), int(f[11])
        if bl < 2000:                      # 2 kb 이상 붙은 것만 flank 로 본다
            continue
        if q not in best or bl > best[q][0]:
            best[q] = (bl, (qs + qe) // 2, ql)
    side = {}
    for q, (bl, mid, ql) in best.items():
        side[q] = "L" if mid < ql / 2.0 else "R"
    return side


def trim_fastq(src, dst, lim, side, seed):
    rnd = random.Random(seed)
    n = ntr = 0
    cut_flank = cut_rand = 0
    with gzip.open(src, "rt") as f, gzip.open(dst, "wt", compresslevel=1) as g:
        while True:
            h = f.readline()
            if not h:
                break
            sq = f.readline().rstrip("\n")
            pl = f.readline()
            ql = f.readline().rstrip("\n")
            n += 1
            if len(sq) > lim:
                nm = h[1:].split()[0]
                s = side.get(nm)
                if s is None:
                    s = "L" if rnd.random() < 0.5 else "R"
                    cut_rand += 1
                else:
                    cut_flank += 1
                if s == "L":               # flank 가 앞 -> 앞 L bp 를 남기고 뒤를 자름
                    sq, ql = sq[:lim], ql[:lim]
                else:                      # flank 가 뒤 -> 뒤 L bp 를 남김
                    sq, ql = sq[-lim:], ql[-lim:]
                ntr += 1
            g.write(h); g.write(sq + "\n"); g.write(pl); g.write(ql + "\n")
    return n, ntr, cut_flank, cut_rand


def run_asm(smp, cond, fq):
    d = tempfile.mkdtemp(prefix="trp_", dir=H)
    try:
        with io.open(os.path.join(d, "hifiasm.log"), "w") as lg:
            rc = subprocess.call([HIFIASM, "--ont", "-o", smp, "-t", TH, fq],
                                 cwd=d, stdout=lg, stderr=subprocess.STDOUT)
        gfa = [x for x in os.listdir(d) if x.endswith(".p_ctg.gfa")]
        if rc != 0 or not gfa:
            shutil.copy(os.path.join(d, "hifiasm.log"),
                        os.path.join(W, "log", "%s.%s.hifiasm.log" % (smp, cond)))
            return "FAIL_rc%d" % rc
        out = subprocess.run([PY, HIARR, "%s_%s" % (smp, cond), d, os.path.join(W, "arr")],
                             capture_output=True, text=True)
        with io.open(os.path.join(W, "log", "%s.%s.hiarr.log" % (smp, cond)), "w") as lg:
            lg.write(out.stdout + "\n" + out.stderr)
        return out.stdout.strip() or "NO_ARR"
    finally:
        shutil.rmtree(d, ignore_errors=True)


for Ln in io.open(LIST):
    Ln = Ln.strip()
    if not Ln:
        continue
    smp, cmax = Ln.split()[0], int(Ln.split()[1])
    fq = find_reads(smp)
    if not fq:
        print("%s\tNO_READS" % smp, flush=True)
        continue
    lim = int(THR * UNIT * cmax)
    tmp = tempfile.mkdtemp(prefix="trs_", dir=H)
    try:
        side = flank_side(fq, tmp)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    sub = os.path.join(W, "%s.trim.fastq.gz" % smp)
    n, ntr, cf, cr = trim_fastq(fq, sub, lim, side, 20260928)
    print("%s\tread %d\tc_max %d\t한도 %d bp\t자른 read %d (%.2f%%)\tflank기준 %d\t무작위 %d\tflank붙은 read %d"
          % (smp, n, cmax, lim, ntr, 100.0 * ntr / max(1, n), cf, cr, len(side)), flush=True)
    res = run_asm(smp, "trim", sub)
    os.remove(sub)
    line = "%s\ttrim\t%d\t%d\t%d\t%d\t%s" % (smp, n, ntr, cf, cr, res)
    print(line, flush=True)
    with io.open(SUM, "a", encoding="utf8") as g:
        g.write(line + "\n")
print("DONE")
