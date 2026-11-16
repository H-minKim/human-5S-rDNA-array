# -*- coding: utf-8 -*-
# Stage 07, the merge: the reads are collected per unit and a consensus variant called at each position, for the other species only.
#
#   Why it is needed: the variant list of a single read is mostly ONT error (measured in gorilla: a median of 94 per 2,248 bp, 4.2%).
#     Only by collecting the reads covering each unit and dividing by the depth at each position do the real variants remain.
#
#   The method, the same depth-threshold consensus rule as the human pipeline:
#     1) a unit is (array_id, array_pos), and the read fragments covering it are collected, excluding NOALN;
#     2) the depth at a position is the number of fragments whose cov_start_in_ref..cov_end_in_ref spans it;
#     3) for each variant token, the support count over the depth at that position is the support fraction, and it is accepted at a support fraction of at least THR and a depth of at least MINDEP;
#     4) the accepted tokens are collected in positional order as that unit's consensus variant list.
#   No filling: an uncovered position is simply left empty, with cons_cov_start and cons_cov_end recording how far the coverage reached.
#   No mapq filter: mapq 0 is normal inside a tandem array.
#
#   The strict version (varlist_strict) keeps only substitutions and indels longer than 3 bp. ONT frequently produces homopolymer deletions,
#     so 1-3 bp deletions cross the threshold (measured in gorilla: 121 of the 204 accepted deletions are absent from the assembly,
#     while 1,514 of 1,593 substitutions, 95%, agree with it). Any per-unit conclusion should use the strict version.
#   Check: how far the result agrees with the assembly variants for the same unit, from the g8 map, is counted and reported, separately for the full and the strict version.
#   Sensitivity: the variant count at THR of 0.20, 0.35 and 0.50 is reported alongside.
#
#   Outputs: <out>/<species>__<sample>.unit_consensus.tsv, one row per unit,
#          and <out>/<species>__<sample>.site_support.tsv, one row per variant with a support fraction of at least REPORT.
#   Arguments: the read-unit table, the reference unit folder, the map folder, the output folder, [THR = 0.35], [MINDEP = 3], [REPORT = 0.10].
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


# The input is in read order and is not grouped by unit, so it is traversed once, accumulating only what each unit needs
#   (holding whole rows would come to several GB).
COV = collections.defaultdict(list)      # (array_id, array_pos) -> [(coverage start, coverage end), ...]
VAR = collections.defaultdict(collections.Counter)
NUMS = collections.defaultdict(lambda: [0, 0, 0])   # all, aligned, fully covered


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
