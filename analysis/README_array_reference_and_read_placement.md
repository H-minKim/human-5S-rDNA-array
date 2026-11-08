# How an array reference is built and how reads are placed on it

This chain produces the array reference that every other figure measures against.

## The reference document

`analysis/METHOD_CONSENSUS_MAP_G8_REBUILD.md` is the method document that every script
below names in its header as its reference document. It states the rule each step applies,
and its md5 is recorded in `fig01/MANIFEST_Fig1.tsv`. This file is the index to it, and the
section numbers below are its section numbers.

It is worth reading for one thing in particular: it records the mistakes that set each threshold.
The contig is found by sequence because guessing "the largest contig" put the array on the wrong
contig for HG01943, where it is the 77th. The flank probes stop at the junction and do not include
any array, because the array is tandem and a probe containing it aligns anywhere.

## The chain, step by step

| § | step | script |
|---|---|---|
| — | pull 5S reads out of each BAM through a k-mer gate that preserves MM/ML | `readextract_5s_read_gate.py`, with `readextract_{hprc_ont,giab_ont_region,giab_hits_and_flank,cancer,hg008}_*.sh` per cohort |
| — | assemble, for individuals with no public assembly | `refset_assembly_queue.sh` (`hifiasm --ont`, one individual at a time, intermediates deleted immediately — a run needs 676 MB and the node had 43 GB) |
| §1 | cut the array flank to flank, from a public assembly | `refset_step1_cut_array_public_assembly.py` |
| §1 | the same cut, from the assemblies of this work (the GFA `S` lines) | `refset_step1_cut_array_our_assembly.py` |
| §1 | count units in a public assembly without storing sequence | `refset_step1_count_units_public_assembly.py` |
| §1 | check that the cut landed on the junction | `refset_validate_canon_flank.py` |
| §2 | index every cut sequence by its unit count | `refset_step2_sequence_index.py` |
| §3 | assign one sequence to each haplotype slot | `refset_step3_assign_slots.py`, with `refset_step0_*` for the pairing and fingerprint steps that feed it |
| §3 | repair the few references that needed it | `refset_step3b_canon_repair_*.py` |
| §4–5 | cut units, call variants, check the sequence back | `refset_step4_5_unit_map_and_variants.py` |
| §6 | build the registry and its tier flags | `refset_step6_registry_build.py`, `refset_step6_registry_add_flags.py` |
| — | place reads on the finished map | `readplacement_align_reads_to_map.py`, `readplacement_methylation_map.py` |

The array-only and junction placement routes themselves live under `upstream/`: `aoplace_v8fix.py`
(in `upstream/fig02_figS2_chain/U3_read_placement_units/v11/`) for reads lying wholly inside the
array, `upstream/v12_fiberseq_node/jxplace.py` and `slotvote.py` for reads that reach a junction, and
`upstream/v12_fiberseq_node/merge_aojv_v2.py` for the rule that merges the two.

## The numbers that matter, and where they come from

These are the values the Methods section quotes. Each is given here with the file it is set in, so
that a reader checking one does not have to search.

**Cutting the array** (`refset_step1_cut_array_*.py`, §1-1 to §1-3)

- The probes are two 50 kb windows cut from CHM13v2.0 immediately *outside* the junctions, never
  including array: `PROX50k = chr1:227,694,884-227,744,883` and
  `DIST50k = chr1:228,025,901-228,075,900`. They are aligned with `minimap2 -x asm20`, and the array
  is the span between their two **inner** ends. A contig where only one probe aligns is logged as
  `NO_BOTH_PROBE` and discarded.
- The contig is found by sequence: eight 25-mers from the 5S gene (offsets 5, 30, 55, 80, both
  strands), and the four contigs with the most hits are tested.
- The two routes differ only in how much probe alignment they demand, because those assemblies are
  short contigs rather than chromosomes: **≥ 20,000** matching bases (PAF column 10) for a public
  assembly, **≥ 5,000** for ours, where ours also keeps secondary hits (`-N 20 -p 0.1`), drops
  segments under 4 kb and treats contigs whose lengths agree to within 200 bp as one.
- For those assemblies the `hap1`, `hap2` and `primary` `p_ctg` sets are **all** read at this step;
  which haplotype a contig belongs to is settled later, at slot assignment.
- Units are counted from the 5S-gene blastn search (`-task blastn -dust no -word_size 9
  -evalue 1e-5`) as HSP start spacings of 2,100-2,400 bp. §2 of the reference document states the
  unit count as **that number of spacings plus one**, and warns against deriving it as `bp / 2240`,
  which was wrong by a factor of two on a sequence whose chain was broken.

**Placing reads** (`readplacement_align_reads_to_map.py`)

- The reference for placement is one sequence set holding **both haplotypes of the individual and
  both 50 kb flank probes**. Both haplotypes are present so a read cannot be pulled onto the wrong
  one by the absence of the right one; the flanks are present so a read running out of the array
  aligns to flank instead of being soft-clipped, which in an earlier version dropped the alignment
  rate to 31-89%.
- `minimap2 -a -x map-ont` for ONT and `-x map-hifi` for PacBio HiFi, with `--secondary=no` and `-y`.
  The platform is taken from the **read name**, not the directory: a cancer directory of HiFi reads
  was once run under `map-ont` and COLO829's alignment identity fell to 69% against 97-99% elsewhere.
- A read counts as reaching the flank when ≥ 1 kb aligns to a flank probe **and** that alignment
  comes within 50 bp of the junction (`MINFL = 1000`, `JXWIN = 50`). The 50 bp is the trough of a
  measured bimodal distribution; at 2 kb a 13-individual pilot called 27 flank-to-flank reads instead
  of 19 without gaining a true one, the extras sitting on a 5.9 kb repeat block inside the proximal
  probe.
- Mapping quality is kept as a column rather than used as a filter at this stage, because the two
  haplotypes of one individual are about 99.6% identical and a read that genuinely cannot be told
  apart gets MAPQ 0.

**Assigning a read to a haplotype** (`aoplace_v8fix.py`, `slotvote.py`, `merge_aojv_v2.py`)

- *Reads that reach a junction* are assigned by a vote at the positions where the two haplotypes'
  flanks differ; mismatches anywhere else are treated as sequencing error. Only flank is used,
  because the array is tandem and both haplotypes align equally well along it — including array
  sequence in the vote dropped accuracy from 100% to 84.7%. The vote is accurate to 99.3% at a margin
  of 0.90.
- *Reads lying wholly inside the array* are assigned by a **map-margin difference**. The read is
  placed as a chain of K consecutive units against the map of every individual in turn, and the
  assignment is accepted only when Δ = (margin on its own individual's map) − (best margin on any
  other individual's map) is at least 10. The difference, not the margin, is the axis that
  discriminates: unit indices run from one to a few hundred, so a high absolute margin arises by
  chance, and an absolute-margin gate raised the null rate from 8.2% to 16.2%. The median margin on
  another individual's map is 1, against 25 on the right one. On held-out data the gate keeps 95.3%
  of reads, is correct for 98.08%, and prefers another individual's map in 0.00%.
- K is fixed per dataset before the data are seen, from what the reads of that dataset can support:
  **K = 16** with a span of at least 36 kb for the whole-genome ONT reads, **K = 12** and 26,880 bp
  for Fiber-seq ONT, and **K = 4** with Δ ≥ 5 and 8,960 bp for the shorter PacBio reads. The
  "K = 12" named in `merge_aojv_v2.py` is the Fiber-seq value, not a different setting for the main
  cohort.
- Where the two routes disagree, votes resting on fewer than five discriminating positions are
  dropped (`MINSITE = 5`; the false-positive rate under permuted labels is 14.3% below five
  positions, 2.4% at five to nine and 0.0% at ten or more) and the remaining conflicts follow the
  chain, because every conflicting read is long enough to lie inside the range the chain was tested
  over. Each such case is written to a log rather than resolved silently — merging the two tables
  with `sort -u` instead had let six of 803 overlapping reads be decided by alphabetical order.