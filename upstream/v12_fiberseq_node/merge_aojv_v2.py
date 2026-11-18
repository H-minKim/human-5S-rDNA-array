# -*- coding: utf-8 -*-
"""Merge the array-only placement and the junction placement by rule

  A plain sort -u is not a merge: when the two sources give one read different slots, both lines survive
    and ao_g8cut takes the line it reads first (always slot "1" in lexical order), which is chance and not a rule.
    (measured : of 803 overlapping reads 6 conflict, and all 6 are cut to slot 1 in the output)

  Conflict resolution
    (1) votes with fewer than 10 discriminating sites are excluded: under a label-shuffling null, votes with <5 sites
       have 14.3% false positives and votes with 10 or more 0.0% (a margin of 1.000 on 4 sites arises by chance 1 time in 7)
    (2) remaining conflicts follow the array-only chain: every conflicting read has >=14 units,
       inside the range where the chain gate (K=12) was tested. The vote exists to fill the spans where the chain fails,
       not to overturn the chain.
    (3) the outcome and its basis are always logged (nothing is decided silently)

  usage : merge_aojv.py <array only tsv...> --jx <junction tsv> --vote <slotvote tsv> --out <merged> --log <log>
"""
import io
import os
import sys
import collections

NL = "\n"
TAB = "\t"
# the threshold is the environment variable MINSITE, default 10; the value used is 5 :
#   False positives in the null: 14.3% at fewer than 5 discriminating positions, 2.4% at 5-9 and 0.0% at 10 or more, so 5-9 is kept (about 3% of those recovered may be misassigned).
MINSITE = int(os.environ.get("MINSITE", "10"))
# the read key length: 16 for ONT and 0, the whole name, for PacBio, whose first 16 characters are the movie and would collide in every case
KEYLEN = int(os.environ.get("KEYLEN", "16"))


def rkey(x):
    return x[:KEYLEN] if KEYLEN else x

args = sys.argv[1:]
ao_files, jx, vote, out, log = [], None, None, None, None
i = 0
while i < len(args):
    a = args[i]
    if a == "--jx":
        i += 1
        jx = args[i]
    elif a == "--vote":
        i += 1
        vote = args[i]
    elif a == "--out":
        i += 1
        out = args[i]
    elif a == "--log":
        i += 1
        log = args[i]
    else:
        ao_files.append(a)
    i += 1

# the number of discriminating positions in the vote: (individual, read16) -> count
NSITE = {}
if vote and os.path.exists(vote):
    for k, L in enumerate(io.open(vote, encoding="utf8")):
        if k == 0:
            continue
        f = L.rstrip(NL).split(TAB)
        if len(f) < 8:
            continue
        NSITE[(f[0], rkey(f[1]))] = int(f[4])


def load(paths, tag):
    """(individual, read16) -> (slot, [rows...])"""
    D = collections.defaultdict(lambda: [None, []])
    for p in paths:
        if not os.path.exists(p):
            continue
        for L in io.open(p, encoding="utf8"):
            f = L.rstrip(NL).split(TAB)
            if len(f) < 5:
                continue
            k = (f[0], rkey(f[1]))
            if D[k][0] is None:
                D[k][0] = f[2]
            D[k][1].append(L.rstrip(NL))
    return D


AO = load(ao_files, "AO")
JX = load([jx], "JX") if jx else {}

lg = io.open(log, "w", encoding="utf8", newline=NL) if log else sys.stderr
lg.write("# merge log (conflict resolution rule : positions<%d excluded from the vote, remaining conflicts are array only)\n" % MINSITE)
lg.write("sample\tread16\tao_slot\tjx_slot\tn_site\t채택\t사유\n")

T = collections.Counter()
rows = []
for k in sorted(set(list(AO.keys()) + list(JX.keys()))):
    a = AO.get(k)
    b = JX.get(k)
    if a and not b:
        rows.extend(a[1])
        T["array only 만"] += 1
        continue
    if b and not a:
        ns = NSITE.get(k, 0)
        if ns < MINSITE:
            T["junction 만인데 자리부족 → 버림"] += 1
            lg.write("%s\t%s\t-\t%s\t%d\t버림\t자리<%d\n" % (k[0], k[1], b[0], ns, MINSITE))
            continue
        rows.extend(b[1])
        T["junction 만"] += 1
        continue
    # present in both
    if a[0] == b[0]:
        # the same slot: the array-only placement is taken as the base and the units the junction placement additionally found are added
        seen = set(x.split(TAB)[3] for x in a[1])
        extra = [x for x in b[1] if x.split(TAB)[3] not in seen]
        rows.extend(a[1])
        rows.extend(extra)
        T["겹침·같은슬롯"] += 1
        T["  junction 이 더한 유닛"] += len(extra)
        continue
    ns = NSITE.get(k, 0)
    rows.extend(a[1])                      # rule 2: follow the array-only placement
    T["슬롯충돌 → array only 채택"] += 1
    lg.write("%s\t%s\t%s\t%s\t%d\tAO\t%s\n"
             % (k[0], k[1], a[0], b[0], ns,
                ("자리<%d 라 투표 배제" % MINSITE) if ns < MINSITE else "사슬 범위(유닛>=12) 우선"))

with io.open(out, "w", encoding="utf8", newline=NL) as o:
    for r in rows:
        o.write(r + NL)

lg.write("\n# summary\n")
for kk, vv in sorted(T.items()):
    lg.write("#   %-34s %d\n" % (kk, vv))
lg.write("#   %-34s %d\n" % ("최종 행", len(rows)))
if log:
    lg.close()
print("  " + " · ".join("%s %d" % kv for kv in sorted(T.items())) + " · 행 %d" % len(rows))
