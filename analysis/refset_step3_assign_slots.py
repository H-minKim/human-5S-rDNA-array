# -*- coding: utf-8 -*-
# Stage 3: slot assignment, one file per haplotype.
#   Step 3: hap_final2.tsv plus seqidx_all.tsv -> hap_seq.tsv -> hap_seq_node_v8.tsv.
#   This script fixes the slot numbers.

# -*- coding: utf-8 -*-
# Assign an actual sequence to each fixed haplotype.
#
#  Inputs:
#    hap_final2.tsv, 314 individuals by slot with the unit count and the basis of the choice, in the unit frame;
#    seqidx_all.tsv, an exhaustive index of the 1,516 .fa files in cmp4d, in the unit frame.
#
#  The assignment rules:
#    1) Only a .fa whose unit count is within 1 of the fixed value is a candidate.
#       (The basis for choosing the haplotype was settled at stage 1; what remains here is only
#        choosing the most accurate sequence among those with that unit count.)
#    2) Choose in order of sequence accuracy: public assembly, then our hifiasm, then a spanning read.
#       An assembly is a polished sequence at about 99.9% or better, while a spanning read is a raw ONT read at about 99%.
#       Even where the spanning read won on unit count, if an assembled sequence gives the same count
#         that sequence is used: the spanning read is right about the count, the assembly about the bases.
#    3) Within a grade, the higher chain continuity wins, then the greater bp.
#    4) The same file is never assigned to both slots of one individual.
#
#  A slot must be backed by an actual sequence.
#    Where no .fa corresponds to the cluster stage 1 ranked first or second, the next cluster is taken.
#    This applies in two cases:
#      NA19120 takes 171 in place of 91, the spanning read's chain count: no 91-unit sequence exists at all,
#               the tier 1 sequence has 96 and hifiasm 171, so 91 is an undercount of the 96.
#      HG01028 takes 40 in place of 45: the tier 1 reference built from that same spanning read holds 40 units.
#    Where no sequence is found at all, the slot keeps its unit count and is marked as having no sequence
#    (a cluster supported only by a recombined-spanning read, whose sequence is synthetic and cannot be used).
#
#  A recombined-spanning sequence is synthetic, exactly N x 2240 bp, and is never used as a sequence.
import io
import os
import collections

S = os.environ["RDNA_S"]
NL = "\n"
TOL = 1
KPRI = {"pub": 0, "hifi": 1, "span": 2}
KNAME = {"pub": "공개조립본", "hifi": "우리hifiasm", "span": "관통"}

IDX = collections.defaultdict(list)
for i, L in enumerate(io.open(os.path.join(S, "seqidx_all.tsv"), encoding="utf8")):
    if i == 0:
        continue
    f = L.rstrip(NL).split("\t")
    if len(f) < 8:
        continue
    IDX[f[1]].append({"file": f[0], "kind": f[2], "dir": f[3],
                      "bp": int(f[4]), "u": int(f[5]), "cp": float(f[6])})

HAP = collections.defaultdict(list)
for i, L in enumerate(io.open(os.path.join(S, "hap_clusters.tsv"), encoding="utf8")):
    if i == 0:
        continue
    f = L.rstrip(NL).split("\t")
    HAP[f[0]].append({"rank": int(f[1]), "u": int(f[2]), "n": int(f[3]),
                      "strong": f[4], "nread": int(f[5]), "tier": f[6],
                      "m": f[7], "twin": f[8]})

OUT = io.open(os.path.join(S, "hap_seq.tsv"), "w", encoding="utf8", newline=NL)
OUT.write("sample\tslot\tunits_final\trank\ttier\tmethods\t"
          "seq_file\tseq_kind\tseq_dir\tseq_bp\tseq_units\tchain_pct\n")
ok = collections.Counter()
miss = []
kcnt = collections.Counter()
demoted = []
for s in sorted(HAP):
    used = set()
    got = []
    left = sorted(HAP[s], key=lambda x: x["rank"])
    # 1) fill the slots from the clusters that have a sequence, in rank order
    for h in left:
        if len(got) >= 2:
            break
        cand = [c for c in IDX.get(s, [])
                if abs(c["u"] - h["u"]) <= TOL and c["file"] not in used]
        if not cand:
            continue
        cand.sort(key=lambda c: (KPRI[c["kind"]], abs(c["u"] - h["u"]), -c["cp"], -c["bp"]))
        b = cand[0]
        used.add(b["file"])
        got.append((h, b))
        if h["rank"] > 2:
            demoted.append((s, h, b))
    # 2) if that leaves fewer than 2, a higher-ranked cluster with no sequence is kept on its unit count alone
    if len(got) < 2:
        for h in left:
            if len(got) >= 2:
                break
            if any(g[0]["u"] == h["u"] for g in got):
                continue
            got.append((h, None))
    for i, (h, b) in enumerate(got):
        if b is None:
            ok["서열없음"] += 1
            miss.append((s, h))
            OUT.write("%s\t%d\t%d\t%d\t%s\t%s\t-\t-\t-\t-\t-\t-%s"
                      % (s, i + 1, h["u"], h["rank"], h["tier"], h["m"], NL))
        else:
            ok["배정"] += 1
            kcnt[b["kind"]] += 1
            OUT.write("%s\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%.1f%s"
                      % (s, i + 1, b["u"], h["rank"], h["tier"], h["m"], b["file"],
                         KNAME[b["kind"]], b["dir"], b["bp"], b["u"], b["cp"], NL))
OUT.close()

tot = ok["배정"] + ok["서열없음"]
print("=" * 80)
print("  Stage 2. Sequence assignment: confirmed haplotypes %d" % tot)
print("=" * 80)
print("  sequence assigned: %d (%.1f%%), no sequence %d" %
      (ok["배정"], 100.0 * ok["배정"] / tot, ok["서열없음"]))
print()
print("  source of assigned sequences")
for k, v in kcnt.most_common():
    print("     %-12s %4d (%5.1f%%)" % (KNAME[k], v, 100.0 * v / ok["배정"]))
print()
print("  === cases that fell to the next rank for lack of a sequence ===")
for s, h, b in demoted:
    print("     %-11s rank%d adopted %3d units (%s): the higher-ranked cluster had no sequence" %
          (s, h["rank"], b["u"], KNAME[b["kind"]]))
print("  cases that fell: %d" % len(demoted))
print()
if miss:
    print("  slots without sequence %d (only the unit count kept: recombined spanning alone, etc.)" % len(miss))
    for s, h in miss:
        print("     %-12s %3d units  [%s]  sequence candidates of that individual: %s" %
              (s, h["u"], h["m"],
               ",".join("%s%d" % (c["kind"], c["u"]) for c in sorted(
                   IDX.get(s, []), key=lambda c: c["u"])) or "없음"))
