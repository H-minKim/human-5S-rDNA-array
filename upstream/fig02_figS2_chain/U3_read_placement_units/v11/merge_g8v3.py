# -*- coding: utf-8 -*-
# G8v3 placement merge (2026-09-11), by data source, read-stored, coding part included
#   Inputs (all stored placement reads, all passed the source check):
#     junction   : JX_G8_UNITS_main_* (610 slots whose old reference = latest canonical) + JX_G8_UNITS_c2_* (6 slots with consensus2 canonical)
#                  columns: sample read slot job source strand array_pos f0 f1 bnd_start bnd_end verify varlist
#     ARRAY_ONLY : AO_G8_UNITS_* columns: sample read read16 slot cfg source strand array_pos g6 f0 f1 bnd_s bnd_e verify varlist score inP
#   Rules (same as mergemap7):
#     - ARRAY_ONLY source is the single one v7 used for each individual (ao_g8cut already picks it that way), DUAL reads excluded per source
#     - a junction position is used only with 3 or more reads at the same (individual, source, slot, position) (original rjfreq2 rule)
#     - consensus = token read count >= 0.35 x depth, compared on UNIT, depth 3 or more, only units with map check OK/OK_FILLED
#     - read units only OK/OK_FILLED/OK_EDGE/OK_MIDDEL (input = .mid.tsv after recall_edge -> recall_mid). NO_FULL_READEND/NOALN/SEQ_MISMATCH are removed from the depth
#   Changes: data source added to the key, junction excludes HG00384 and HG01150 (their assignment table is an alignment self-reference readslot)
#   Output: READMAP_G8V3_BYSOURCE.tsv (per position read ids, consensus, map, spacer and gene part separately)
#   Arguments: collected folder (JX_G8_UNITS_*, AO_G8_UNITS_*), cmp folder (DUAL, G8v3 map)
import io
import os
import sys
import glob
import collections

NL = chr(10)
TAB = chr(9)
T = 0.35
RD_OK = tuple(os.environ.get("RDOK", "OK,OK_FILLED,OK_EDGE,OK_MIDDEL").split(","))   # read units: rule B (edge deletions <=10bp) + marking of deletions in middle read units
IN, W = sys.argv[1], sys.argv[2]
JX_EXCL = set(["HG00384", "HG01150"])


def num(t):
    n = ""
    for ch in t:
        if ch.isdigit():
            n += ch
        else:
            break
    return int(n) if n else -1


DU = {}
for c, fn in (("AOJV3", "DUAL20.tsv"), ("AOJV5", "DUAL5.tsv"), ("AOJV8", "DUAL8.tsv"), ("AOJVG", "DUALG.tsv")):
    DU[c] = set(tuple(L.rstrip(NL).split(TAB)[:2]) for L in io.open(os.path.join(W, fn), encoding="utf8") if L.strip())
M = {}
f = io.open(os.path.join(W, "g8v3", "CONSENSUS_MAP_G8_full.tsv"), encoding="utf8")
ix = dict((n, i) for i, n in enumerate(f.readline().rstrip(NL).split(TAB)))
for L in f:
    r = L.rstrip(NL).split(TAB)
    if r[ix["seg_type"]] != "UNIT":
        continue
    ok = r[ix["verify"]] in ("OK", "OK_FILLED")
    v = r[ix["varlist"]]
    M[(r[0], r[1], int(r[ix["array_pos"]]))] = (set() if v == "-" else set(v.split(","))) if ok else None
# junction read units
GENE_OK = ("OK", "OK_FILLED")           # state used for the gene part (2122-2240) consensus: intact units only
                                        #   units rescued by edge deletions get end fill that creates substitutions near 2235-2240 that are not there
JR = collections.defaultdict(dict)       # (smp, src, slot, pos) -> {read: (tokens, intact or not)}
JBAD = collections.Counter()
for p in sorted(glob.glob(os.path.join(IN, "JX_G8_UNITS_*.mid.tsv"))):
    for L in io.open(p, encoding="utf8"):
        r = L.rstrip(NL).split(TAB)
        if len(r) < 13:
            continue
        smp, rid, sl, job, src, pos, vd, vl = r[0], r[1], r[2], r[3], r[4], int(r[6]), r[11], r[12]
        if smp in JX_EXCL:
            JBAD["제외개체"] += 1
            continue
        if vd not in RD_OK:
            JBAD["변이못부름"] += 1
            continue
        JR[(smp, src, sl, pos)][rid] = ((set() if vl == "-" else set(vl.split(","))), vd in GENE_OK)
JR = dict((k, v) for k, v in JR.items() if len(v) >= 3)
# ARRAY_ONLY read units
AR = collections.defaultdict(dict)
ABAD = collections.Counter()
for p in sorted(glob.glob(os.path.join(IN, "AO_G8_UNITS_*.mid.tsv"))):
    for L in io.open(p, encoding="utf8"):
        r = L.rstrip(NL).split(TAB)
        if len(r) < 17:
            continue
        smp, rid, r16, sl, cfg, src, pos, vd, vl = r[0], r[1], r[2], r[3], r[4], r[5], int(r[7]), r[13], r[14]
        if (smp, r16) in DU.get(cfg, set()):
            ABAD["DUAL제외"] += 1
            continue
        if vd not in RD_OK:
            ABAD["변이못부름"] += 1
            continue
        AR[(smp, src, sl, pos)][rid] = ((set() if vl == "-" else set(vl.split(","))), vd in GENE_OK)
O = io.open(os.path.join(W, os.environ.get("OUTN", "READMAP_G8V3_BYSOURCE.tsv")), "w", encoding="utf8", newline=NL)
O.write(TAB.join(["sample", "data_source", "slot", "slot_space", "array_pos", "frame", "depth", "depth_gene", "depth_jx", "depth_ao",
                  "read_varlist", "map_varlist", "agree", "n_TP", "n_FN", "n_FP", "sp_TP", "sp_FN", "sp_FP",
                  "gene_TP", "gene_FN", "gene_FP", "reads_jx", "reads_ao"]) + NL)
G = collections.Counter()
PH = collections.defaultdict(lambda: collections.Counter())
for k in sorted(set(JR) | set(AR)):
    smp, src, sl, pos = k
    j = JR.get(k, {})
    a = AR.get(k, {})
    vals = list(j.values()) + list(a.values())
    dep = len(vals)                                   # spacer depth = all usable read units
    depg = sum(1 for v in vals if v[1])               # gene part depth = intact units only
    c = collections.Counter()
    cg = collections.Counter()
    for tk, full in vals:
        c.update(t for t in tk if num(t) <= 2121)
        if full:
            cg.update(t for t in tk if num(t) > 2121)
    R = set(x for x, n in c.items() if n >= T * dep) | set(x for x, n in cg.items() if depg and n >= T * depg)
    m = M.get((smp, sl, pos), "NOMAP")
    row = [smp, src, sl, "assembly", str(pos), "G8v3", str(dep), str(len(j)), str(len(a)), ",".join(sorted(R)) or "-"]
    row.insert(7, str(depg))
    if isinstance(m, set):
        sp = lambda s: set(t for t in s if num(t) <= 2121)
        ge = lambda s: set(t for t in s if num(t) > 2121)
        tp, fn_, fp = len(R & m), len(m - R), len(R - m)
        s3 = (len(sp(R) & sp(m)), len(sp(m) - sp(R)), len(sp(R) - sp(m)))
        g3 = (len(ge(R) & ge(m)), len(ge(m) - ge(R)), len(ge(R) - ge(m)))
        agree = "EXACT" if fn_ == 0 and fp == 0 else "DIFF"
        row += [",".join(sorted(m)) or "-", agree, str(tp), str(fn_), str(fp)] + [str(x) for x in s3 + g3]
        if dep >= 3:
            G["cmp"] += 1
            G["ex"] += (agree == "EXACT")
            for nm, t3 in (("all", (tp, fn_, fp)), ("sp", s3), ("ge", g3)):
                G[nm + "TP"] += t3[0]
                G[nm + "FN"] += t3[1]
                G[nm + "FP"] += t3[2]
            e = PH[(smp, src, sl)]
            e["n"] += 1
            e["TP"] += tp
            e["FN"] += fn_
            e["FP"] += fp
    else:
        row += ["NOMAP" if m == "NOMAP" else "NOCALL", "NOMAP" if m == "NOMAP" else "NOCALL"] + ["-1"] * 9
    row += [",".join(sorted(j)) or "-", ",".join(sorted(a)) or "-"]
    O.write(TAB.join(row) + NL)
    G["rows"] += 1
O.close()


def pct(p):
    t = G[p + "TP"] + G[p + "FN"] + G[p + "FP"]
    return 100.0 * G[p + "TP"] / t if t else float("nan")


print("junction read 유닛 제외 : %s · ARRAY_ONLY : %s" % (dict(JBAD), dict(ABAD)))
print("자리 %d · UNIT 깊이3+ 지도있음 %d : 변이염기 %.3f%% · 스페이서(1-2121) %.3f%% · 유전자부(2122-2240) %.3f%% · 전체염기 %.4f%% · 완전일치 %.2f%%"
      % (G["rows"], G["cmp"], pct("all"), pct("sp"), pct("ge"),
         100.0 * (1 - (G["allFN"] + G["allFP"]) / float(2240 * G["cmp"])) if G["cmp"] else 0, 100.0 * G["ex"] / G["cmp"] if G["cmp"] else 0))
print("유전자부 TP %d · FN %d · FP %d" % (G["geTP"], G["geFN"], G["geFP"]))
srcs = collections.defaultdict(set)
for (smp, src, sl) in PH:
    srcs[smp].add(src)
print("출처별 하플로타입 : %s" % dict(collections.Counter(k[1] for k in PH)))
print("출처가 둘 이상인 개체 %d : %s" % (sum(1 for v in srcs.values() if len(v) > 1), " ".join("%s[%s]" % (s, "+".join(sorted(v))) for s, v in sorted(srcs.items()) if len(v) > 1)))
