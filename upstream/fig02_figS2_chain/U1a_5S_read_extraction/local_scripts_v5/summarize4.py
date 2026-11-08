#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ONT summary v4: tRNA removal + CHM13 alignment replaced by the TUC tandem map.

Reasons for the change (all measured)
  - CHM13 alignment took 44-54% of the time per individual, and its PAF was used only for locus assignment and
    for choosing 'which reads to keep'.
  - Replaced by the TUC tandem map (90kb): on the same 8,000 ONT reads, CHM13 312s vs TUC 60s = 5.2x.
    At the keep threshold tandem_span>=500bp, recall against CHM13 1.000, false positives 0.79%.
    (Raising the threshold to 4,480 gives 0 false positives but recall falls to 0.719, losing short array reads.
     Recall comes first at the keep step, so 500 is used.)
  - tRNA is excluded from this project. Removing BBDuk at the extraction step gives 1.60x
    (same 1,394MB input 20.58s -> 12.85s), and the tRNA identification blastn of the summary step also goes away.
  - tRNA hit read IDs are still recorded only (re-extraction possible later). Their sequences are not kept.
"""
import sys, os, gzip, argparse, subprocess, collections

ARR_MIN = 500          # keep threshold (recall 1.000 measured)
NL = chr(10)

ap = argparse.ArgumentParser()
ap.add_argument("--sample", required=True)
ap.add_argument("--tech", required=True)
ap.add_argument("--hits", required=True)          # combined <S>.5s.tsv
ap.add_argument("--paf", required=True)           # tandem map PAF (not CHM13)
ap.add_argument("--fastq", required=True)
ap.add_argument("--nts", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--array", required=True)
ap.add_argument("--trnaids", default="")          # list of tRNA hit IDs (recorded only, if given)
ap.add_argument("--threads", default="8")
a = ap.parse_args()


def die(m):
    sys.stderr.write("FATAL: %s%s" % (m, NL)); raise SystemExit(9)


for f in (a.hits, a.paf, a.fastq, a.nts):
    if not os.path.exists(f): die("입력 없음 %s" % f)

# ---- 5S hits ----
per = {}
for L in open(a.hits):
    if L.startswith("read_id"): continue
    c = L.rstrip(NL).split("\t")
    if len(c) < 7: continue
    sp = [] if c[5] == "-" else [int(x) for x in c[5].split(",")]
    per[c[0]] = dict(L=int(c[1]), nhsp=int(c[2]), ngene=int(c[3]),
                     pid=float(c[4]), sp=sp, mm=(c[6] == "Y"))
if not per: die("5S 히트 0")

trids = set()
if a.trnaids and os.path.exists(a.trnaids):
    for L in open(a.trnaids):
        if L.startswith("read_id"): continue
        trids.add(L.split("\t")[0])
sys.stderr.write("[1] 5S read %d · tRNA 히트 ID %d (기록만)%s" % (len(per), len(trids), NL))

# ---- tandem map PAF -> array span ----
# A truncated PAF is not passed over silently (if minimap2 dies from OOM the last line is cut).
span = {}
n_bad = 0
for L in open(a.paf):
    f_ = L.rstrip(NL).split("\t")
    if len(f_) < 12:
        n_bad += 1; continue
    try:
        s_ = int(f_[3]) - int(f_[2])
    except ValueError:
        n_bad += 1; continue
    if f_[0] not in span or s_ > span[f_[0]]: span[f_[0]] = s_
if n_bad:
    die("PAF 손상 (깨진 행 %d) — minimap2 미완주. 이 개체는 재처리해야 한다." % n_bad)
sys.stderr.write("[2] 탠덤 정렬 read %d%s" % (len(span), NL))

# ---- fasta ----
fa = a.fastq + ".fa"
with open(fa, "w") as fh, open(a.fastq) as f:
    i = 0
    for L in f:
        i += 1
        if i % 4 == 1: fh.write(">" + L[1:].split()[0] + NL)
        elif i % 4 == 2: fh.write(L)

# ---- NTS axis ----
db = "/tmp/nts_%s" % a.sample
subprocess.run(["makeblastdb", "-in", a.nts, "-dbtype", "nucl", "-out", db],
               capture_output=True, check=True)
o = "/tmp/nts_%s.tsv" % a.sample
r = subprocess.run(["blastn", "-task", "blastn", "-word_size", "11", "-dust", "no",
                    "-evalue", "1e-10", "-db", db, "-query", fa, "-num_threads", str(a.threads),
                    "-max_target_seqs", "10000", "-outfmt", "6 qseqid sseqid pident length",
                    "-out", o])
if r.returncode != 0: die("NTS blastn 실패")
nts = set()
for L in open(o):
    f_ = L.rstrip(NL).split("\t")
    if int(f_[3]) >= 200 and float(f_[2]) >= 85: nts.add(f_[0])

# ---- classification ----
lab = collections.Counter(); rows = []; arr_ids = set()
for rid, d5 in per.items():
    sp_ = span.get(rid, 0)
    periodic = any(2100 <= x <= 2400 for x in d5["sp"])
    if sp_ >= ARR_MIN:
        c5 = "ARRAY"; arr_ids.add(rid)
    elif d5["ngene"] >= 1:
        c5 = "DISPERSED"
    else:
        c5 = "AMBIGUOUS"
    lab[c5] += 1
    rows.append((rid, d5["L"], "Y", "Y" if rid in trids else "N", d5["ngene"],
                 "%.1f" % d5["pid"], "Y" if periodic else "N", "Y" if rid in nts else "N",
                 sp_, "%.2f" % (sp_ / 2240.0), "Y" if d5["mm"] else "N", c5))

with gzip.open(a.out + ".reads.tsv.gz", "wt") as fh:
    fh.write("read_id\tread_len\thas_5S\thas_tRNA\tn_gene\tbest_pid\tperiodic\thas_NTS"
             "\ttandem_span\ttandem_units\thas_MM\tlabel_5S" + NL)
    for r_ in rows: fh.write("\t".join(map(str, r_)) + NL)

# ---- keep ARRAY read sequences ----
n_arr = 0
with gzip.open(a.array, "wt") as gz, open(a.fastq) as f:
    keep = False; i = 0
    for L in f:
        i += 1
        if i % 4 == 1: keep = L[1:].split()[0] in arr_ids
        if keep:
            gz.write(L)
            if i % 4 == 0: n_arr += 1

with open(a.out + ".summary.tsv", "w") as fh:
    fh.write("sample\ttech\tn_read\tn_5S\tn_tRNA_id\tARRAY\tDISPERSED\tAMBIGUOUS"
             "\tn_NTS\tn_MM\tn_periodic\tn_array_kept" + NL)
    fh.write("\t".join(map(str, [a.sample, a.tech, len(per), len(per), len(trids),
                                 lab["ARRAY"], lab["DISPERSED"], lab["AMBIGUOUS"], len(nts),
                                 sum(1 for r_ in rows if r_[10] == "Y"),
                                 sum(1 for r_ in rows if r_[6] == "Y"), n_arr])) + NL)
for f_ in (fa, o):
    try: os.remove(f_)
    except OSError: pass
for e in (".ndb", ".nhr", ".nin", ".not", ".nsq", ".ntf", ".nto"):
    try: os.remove(db + e)
    except OSError: pass
print("%s: 5S read=%d ARRAY=%d(보존 %d) DISPERSED=%d NTS=%d tRNA_id=%d"
      % (a.sample, len(per), lab["ARRAY"], n_arr, lab["DISPERSED"], len(nts), len(trids)))
