# The G8 consensus map, and how to rebuild it from an updated assembly (5S rDNA)


> **What this document is for.** The G8 in hand (`CONSENSUS_MAP_G8_full.tsv`, md5 `f1d6103c…`) was
> cut from an **older assembly**. The assembly has since been updated, so **do not use the G8 table
> as it stands**. Build it again **from the new assembly, from the beginning**, by the method below.
> The **method** of G8 (cut at the gene end, then check against the sequence) is kept as it is. Only
> the **input** is made anew.

---

## §0 The whole chain, five steps from assembly to map

```
[assembly]  public (HPRC) · the hifiasm assemblies of this work · spanning reads
   │
   ① cut out the array (flank to flank)      asm_units3.py / asm_arr_all.py / hiarr.py
   │   → ~/cmp4d/<source>__<node>__<name>.fa   (old set: 1,516 files)
   ② index the sequences (unit count per file) seqidx.py → seqidx_all.tsv
   ③ assign slots (one file per haplotype)     seqassign.py → hap_seq.tsv
   │   → hap_seq_node.tsv = the (sample, slot, seq_file) three columns of hap_seq.tsv
   ④ unify the orientation (gene on the + strand by majority)   inside consmapG8.py
   ⑤ cut the units, call the variants, check against the sequence   consmapG8.py
       → CONSENSUS_MAP_G8_<TAG>.tsv
```

**If the assembly changes, start again at ①.** When the output of ① changes, ② and ③ change with
it, and when ③ changes the **slot numbers** can change (§3-3, the most dangerous place).

| Step | Reusable? | Why |
|---|---|---|
| references (`unit_canon.fa`, `5S_gene_query.fa`, `probe.fa`) | **as they are** | independent of the assembly |
| the settled unit count (`units_final` of `hap_seq.tsv`) | **as it is** | read based (spanning reads and TUC), independent of the assembly |
| ① array sequences | **again** | the assembly itself changed |
| ② index, ③ assignment, `hap_seq_node.tsv` | **again** | they depend on ① |
| the revcomp column of `CONSENSUS_SUM_G6.tsv` | **discard** | it is keyed to the old files. §4 |
| ⑤ `consmapG8.py` | **same method, new input paths** | §5 |

---

## §1 ① Cutting out the array, flank to flank

The three extractors use **exactly the same method**. If they did not, no comparison between
sources would hold.

| Extractor | Target | Output directory |
|---|---|---|
| `asm_units3.py` | HPRC release2 public assemblies (unit count only) | `~/asmunits3/` |
| `asm_arr_all.py` | the same target, **storing every sequence and coordinate** | `~/asmarr/` |
| `hiarr.py` | the hifiasm assemblies of this work assemblies (the `S` lines of the GFA) | `~/hiasm/arr/` and others |

All three are in the compute node's home directory. A deposited copy of `asm_arr_all.py` is in
`upstream/fig02_figS2_chain/U1c_our_hifiasm/local_scripts_v6_hifiasm/`.

### 1-1 The probes, two 50 kb pieces **outside** the junctions

```
~/asmunits3/probe.fa  (made from CHM13v2.0 with samtools faidx)
  PROX50k = chr1:227,694,884-227,744,883     (PJ − 50000 .. PJ − 1,  PJ = 227,744,884)
  DIST50k = chr1:228,025,901-228,075,900     (DJ + 1 .. DJ + 50000,  DJ = 228,025,900)
```

**The array itself is not put in the probe.** It is tandem, so it would stick anywhere.

### 1-2 Find the contig that holds the array **by sequence**, never by guessing

Take four 25-mers from the 5S gene (at offsets 5, 30, 55 and 80, on both strands, eight in all),
sweep the assembly with `grep`, and take the **four contigs with the most motifs** as the
candidates. (Taking "the last 45 Mb of the largest contig" does not work: in HG01943
the array was on the seventy seventh contig.) For the hifiasm assemblies of this work, look at the `p_ctg` of `hap1`,
`hap2` **and** `primary`.

### 1-3 Align the probes, and the array is **between the inner ends** of the two

```python
paf = minimap2 -x asm20 -t TH --secondary=no  probe.fa  contig.fa
# keep only matches with >= 20,000 matching bases (PAF column 10), the longest one per probe
# both probes must align. If either is missing, record NO_BOTH_PROBE and discard
a = sorted([P.qs, P.qe, D.qs, D.qe])
array = contig[a[1] : a[2]]          # between the two inner coordinates = flank to flank
```

**The boundaries of this `.fa` are the junctions.** No later step touches the flank again.

### 1-4 The file naming rule, which must be kept

```
~/cmp4d/<source directory>__<node mark>__<base name>.fa
for example) asmunits3__58__arr_HG00097_1.fa      arr__58__hiarr_COLO829-BL_1.fa
             tier1ref__58__t1_HG03521_p2_d2.fa
```

`seqidx.py` splits the name on `__` and decides the tier from the **first piece, the source**:

```python
PUB_D  = {"asmunits3", "hgsvcunits"}                                  # public assemblies
SPAN_D = {"tier1ref"}                                                 # spanning reads
HIFI_D = {"arr", "resub", "resub2", "hicheck", "refix", "fill60", "sweep31", "nopurge"}
kind = "pub" if d in PUB_D else "span" if d in SPAN_D else "hifi" if d in HIFI_D else None
if not kind:
    continue          # a new source name is dropped here without a warning
```

> **A trap.** If the updated assembly goes in under a new directory name (say `asm0911`), it
> **drops out of the index in silence**. If a new name is to be used, **add it to `PUB_D` or
> `HIFI_D` first**, then count that the index has as many rows as there are files. A file with
> `len(sq) < 4000` is dropped in silence as well.

### 1-5 Do not mix the new files with the old

`~/cmp4d/` holds 1,516 `.fa` files from the old assembly. **If new files go into the same folder**,
the assignment of ③ can pick an old file. Collect them in a **new folder** (say `~/cmp4d_v2/`) and
do not delete the old one (deleting an existing backup is not allowed).

---

## §2 ② Indexing the sequences, measuring the unit count of each file

`seqidx.py` (split over three nodes, about 6 minutes for
the old set)

```bash
python3 ~/seqidx.py <number of parts> <this part>   # splits the file list by i % parts == this part
# -> ~/seqidxw/seqidx_<TAG>.tsv   columns: file sample kind dir bp units chain_pct n_hsp
```

- The unit count is **the number of gaps between gene HSP start points that fall in 2,100 to 2,400,
  plus one** (the `UNITS` frame, the same as `hap_final2.tsv`).
- **Never convert with `bp/2240`.** There is a precedent of being off by a factor of two on a
  sequence whose chain was broken.
- Sort with **Python `sorted()`** only. The `asorti` of awk sorts as text and scrambles the
  coordinates.
- Change `D = ~/cmp4d` inside the index code to the new folder of §1-5 before running it.

---

## §3 ③ Slot assignment, one `.fa` per haplotype

### 3-1 The rule

| Order | Rule |
|---|---|
| ① | only a `.fa` whose unit count is within **plus or minus one of the settled value (`units_final`)** is a candidate |
| ② | **by sequence accuracy**: public assembly, then the hifiasm assemblies of this work, then spanning reads |
| ③ | within one tier, **chain continuity %**, then **bp** |
| ④ | **the same file is never given to both slots** of one individual |
| ⑤ | a cluster with no candidate `.fa` is **demoted to the next cluster** (`hap_clusters.tsv`, 755 rows) |

- The reason for ②: the unit **count** is accurate from spanning reads, while the **bases** are
  accurate from an assembly. Where an assembled sequence gives the same unit count, that sequence
  is used.
- **A recombined spanning sequence (the old `refseq_*.fa`) is synthetic, exactly N times 2240, and
  is not used as a sequence.**
- Where the hifiasm assemblies of this work gives several contigs, **do not choose by the label (hap1 or hap2)**. For HG002,
  hap1 = 68 was right, hap2 = 77 was false and pri = 103 was right, and the false 77 was a partial
  duplicate overlapping the true 103 by 99.84%. **Go down by length, and discard anything covered
  90% or more by something already chosen. At most two per individual** (`build_refs.py`).

The old result: **616 of 621 slots assigned** (400 public, 198 the hifiasm assemblies of this work, 18 spanning). Against
the answer key (the 201 individuals with a public diploid assembly), both slots agreed within one
unit in **196 of 201, 97.5%**.

### 3-2 Rebuilding the slot assignment

The assignment rule above is the specification; the assignment is reproduced by applying it
rather than by rerunning a stored script. `hap_seq_node.tsv` is the `sample`, `slot` and
`seq_file` columns of `hap_seq.tsv` unchanged (616 of 616 identical). A reimplementation
should be given the original input (`seqidx_all.tsv`) and checked first against the 616 rows
of `hap_seq.tsv`, which it must reproduce exactly.

### 3-3 The slot number, the most dangerous place

`slot` is a number in **`assembly` space**. Every read placement table (`HAP_READ_FINAL.tsv`,
`RJF_*`, `AOJV*`, `READMAP_MERGED_*`) **joins on this number**. **If 1 and 2 swap within an
individual when the slots are assigned again on a new assembly, every read goes onto the wrong
haplotype.**

The precedent (2026-09-06 to 08): the pairing numbers were used without being checked against the
assembly numbers, **36.6% flipped on a coin toss and 83 individuals were discarded entirely from
the junction placement**. In the old set too, the slot and the file suffix **do not agree** (168
cases of slot 2 with `_1`). **Do not infer the slot from the suffix.**

**One of the two must be done:**

- **(a) Keep the old numbers.** Align the new array against each of the two old slot arrays of the
  same individual (`minimap2 -x asm10 -c`, and without `-c` the identity is underestimated) and
  give it the number of **the closer one**. Where the difference to the two old slots is small and
  the call is unclear, **record that individual separately** (no coin tosses).
- **(b) New numbers, and place the reads again.** Anchor the read groups **on the contigs of the
  new assembly** with `slotlib.anchor_slots()`, settle the slots again, and rebuild every
  downstream read table.

Either way, **write a `slot_space` column in the table**, and the reading side confirms it with
`slotlib.require_assembly()`.

---

## §4 ④ Unifying the orientation

**The orientation is decided by a majority vote of the gene HSP strands.** If more HSPs are on the
`-` strand, the sequence is reverse complemented. **Do not flip the coordinates by hand.** Flip the
sequence and **find the HSPs again**.

As it stands, `consmapG8.py` first applies the `revcomp` column (`RCF`) of the old
`CONSENSUS_SUM_G6.tsv` and then checks the strand and flips once more. The final orientation is
right, but the **input is tied to an old file**.

> **On a new assembly, take `RCF` out.** Its key is `(sample, slot)`, so the old value is applied
> even when the file has changed.
> **The more dangerous place is `slotlib.build_ref()`**, which applies `RCF` and **does not check
> the strand again**. If an old `RCF` catches a new file, the read placement reference is built
> **backwards**. Either change `build_ref()` to the same strand majority vote, or build a new
> orientation table and pass it in.

---

## §5 ⑤ G8, cutting the units, calling the variants, checking against the sequence

Code: `consmapG8.py`

### 5-1 Why the cut is at the **end** of the gene

```
canonical unit (unit_canon.fa, 2,240 bp, circular)
[AluY 1-318][NTS1 319-1178][GArich 1179-1551][NTS2 1552-2121][gene 2122-2240]
```

The canonical unit runs `[spacer][gene]`, so **it starts immediately after the gene ends**. G6 cut
at the **start** of the gene, which made it `[gene][spacer]` (a rotation by 119 bp), and because
minimap2 cannot wrap around a circle it **never saw 2122 to 2240** (zero variants in the gene
part).

### 5-2 The boundaries

```python
hits = gene HSPs (blastn -task blastn -word_size 9 -evalue 1e-5 -dust no)
       # where several HSPs fall in the same 1 kb bin (s // 1000), keep only the longest
ends = [e for (s, e, ln) in hits]
for i in range(len(ends) - 1):
    g = ends[i + 1] - ends[i]
    if LO <= g <= HI:                               # 2100 to 2400
        bounds.append((ends[i] + 1, ends[i + 1]))   # 1 based inclusive = the canonical phase
    else:
        nbreak += 1                                 # no unit is made for that interval
```

- From the junction to the first gene end is `PARTIAL_PROX`, and from the last gene end to the
  junction is `PARTIAL_DIST`. **No variant is called in either** (they are not whole units,
  `verify=SKIP`).
- If there are fewer than three gene HSPs, that haplotype is skipped and counted in `STAT` under
  the key for too few gene HSPs, which the code spells in Korean. **Always look at that count.**
- Note that the index of ② counts units by the gaps between gene **starts** while G8 cuts by the
  gaps between gene **ends**. They usually agree, but they can differ where an HSP was clipped, so
  compare `N` with the `units` of the index.

### 5-3 Calling the variants

```bash
minimap2 -x asm20 -c --cs=long -t TH --secondary=no  unit_canon.fa  units.fa
```

- **`asm20`**: on the same boundaries it differs from `asm5` in only 1.36% of cases (indel notation
  at 1455 to 1482) and passes the sequence check far better (`NO_FULL` 117 against 226).
- **`--cs=long`**: the sequence comes back as `=ACGT`. With the `:N` (a count) of the short form,
  **reconstruction is impossible**. The parser takes both forms with
  `([:*+\-=])([A-Za-z0-9]+)`.
- The tokens: `412G` (substitution), `1985D2` (a 2 bp deletion), `1290I1` (a 1 bp insertion,
  **length only**).

### 5-4 The sequence check, in three layers

**① Before starting**, doubt the canonical unit:
```python
if len(CANSEQ) != 2240 or len(GENESEQ) != 119: raise SystemExit(1)
if CANSEQ[2121:2240] != GENESEQ:               raise SystemExit(1)
```

**② Reconstruct and compare**: build the unit sequence again from the variant list and compare it
with the original **character by character**:
```python
if res["reb"] != seq[res["qs"]:res["qe"]]:  -> SEQ_MISMATCH
```

**③ Fill the ends directly**: minimap2 does not extend an alignment where the end mismatches, so it
**eats the variants there** (of 400 uncovered ends, **237 held a real variant**). Compare beyond the
end of the alignment against the sequence itself:
```python
while te + kt < 2240 and qe + kt < len(qseq):    # only within the unit's own sequence
    if qseq[qe + kt] != CANSEQ[te + kt]:
        toks.append("%d%s" % (te + kt + 1, qseq[qe + kt]))
    kt += 1
# the leading end the same way, going backwards
```
The fill cannot cross a unit boundary, so **nothing is counted twice**.

The `verify` values: `OK`, `OK_FILLED` (the ends were filled), `NO_FULL` (the canonical 2,240 was
not fully covered), `SEQ_MISMATCH` and `NOALN`. **Where it is not `OK` or `OK_FILLED`, `varlist` is
`-`**, which means the variants could not be called, not that there are none.

### 5-5 What to change for a new assembly

| Line | Now | Change to |
|---|---|---|
| reading `SEQF` | `~/hap_seq_node.tsv` | the **new** `hap_seq_node` |
| the `.fa` path | `~/cmp4d/` plus the file name | the **new** folder (§1-5) |
| `RCF` | `~/CONSENSUS_SUM_G6.tsv` | **delete it**, the strand majority vote alone (§4) |
| `W` | `~/g8` | a **new** folder (say `~/g8v2`), so the old G8 is not overwritten |

```bash
TAG=full PRESET=asm20 TH=8 LO=2100 HI=2400 python3 consmapG8.py <the new hap_seq_node.tsv>
```

> **Never run a test against a canonical path.** A run that writes where the canonical tables
> live overwrites them in place and nothing warns first. Put a new set in a **new folder**, and
> record the md5 of the old set before running.

---

## §6 The checks to pass before accepting it (all of them, or it is not canonical)

Use the check code of the old G8 as it is:

| # | Check | Code | Value for the old G8 | Pass if |
|---|---|---|---|---|
| 1 | the canonical unit, before starting | inside `consmapG8.py` | passed | it does not stop |
| 2 | the distribution of the check | the summary output | OK 99.15%, SEQ_MISMATCH 0, NOALN 0 | **SEQ_MISMATCH 0 and NOALN 0** |
| 3 | is the gene part visible | awk (the highest variant position) | highest 2240, gene part 12,405 | **highest 2240, gene part above 0** |
| 4 | double counting | an overlap check | overlap 0, the lengths sum to the array in 592 of 616 | **overlap 0** |
| 5 | what the holes are | `holes.py`, `holeid*.py` | 29 places, 80.6% are N gaps | **the N% of every hole, and what it is, recorded** |
| 6 | end artefacts | `chk2113.py` (an independent alignment of 800 bp of tandem) | position 2113 agrees in 99.89% | the top end positions **at or above 99%** |
| 7 | preset sensitivity | `cmpG.py` (asm20 against asm5) | 1.36% | confirm the difference is one of notation |
| 8 | against the index | `N` against `seqidx units` | not applicable | list the individuals that disagree |
| 9 | **against the old G8** | compare the same (sample, slot) | not applicable | **a list of the haplotypes that changed, with the reason** |

Number 9 matters. Only the haplotypes whose assembly changed should differ. **If the map changes
where the assembly did not, something is wrong.**

---

## §7 What was left open in the old G8 (look at it again in the new set)

- **The `LO=2100` threshold discards short units.** There are five holes of exactly `1,820 bp`
  (HG01175 twice, NA20799 twice, NA20806). Whether these are genuinely short units is undecided.
  Decide from the sequence before changing the threshold.
- **`NO_FULL` 117 (0.20%)**, cause unknown.
- **HG00116 slot 1** carries 30,360 bp of distal flank in the middle of the array (hg38
  chr1:228,646,039-228,676,039, MAPQ 60, 93.3% to the distal flank reference). **Whether it is
  there in the new assembly too** would separate a mis-join from the real thing.
- **Six N gaps** (HG01150 127 kb, NA20827 100 kb and others), whether they are filled in the new
  assembly.
- Six badly discordant units (n_var at or above 100), the two most extreme being HG01028 slot 2.

---

## §8 Rerunning the downstream work (once the new map is canonical)

`mergemap3.py`, `readmap_final.py`, `cmpmap3.py`, `cmpmap.py` and `haptest_ao.py` all read the map
from `CONSENSUS_MAP_G6.tsv`. Point them at the new G8, and filter on
**`seg_type != "UNIT"`** (the G6 `UNIT_LAST` and the G8 `PARTIAL_DIST` have different names).
**If the slot numbers of §3-3 were changed by route (b), rebuild the read placement tables first.**

---

## §9 The hard rules

- Do not estimate. Coordinates, lengths and identities are measured against the raw data.
- Never convert with `bp/2240`. The boundaries come from the gene HSPs.
- When building a map, **always check against the sequence**.
- Report both denominators: per variant base and per base of the whole 2,240 bp.
- Node addresses come from `$ANALYSIS_ROOT_INDEX\NODES.md` and nowhere else.
- Anything holding a backslash goes through a file writing tool, not a heredoc, because the path
  breaks in silence.
- Results are recorded in `LATEST_FILES.md` and synchronised with `$WORK_ROOT`. **An existing
  backup is never deleted.**
- The public code is MIT, copyright Hungmin Kim.
