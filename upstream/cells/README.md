# The part of the K562 / GM12878 chain that comes before `cell_fhmm_units.py`

Run `` (the run that made the cell-line tables). Written 2026-10-22.
Covers the chain under Fig. 6 G-J and Suppl. Fig. S9 C-F.

Everything below was read from the compute node and from the `@PG` headers of
the BAMs, read-only. Nothing was recomputed. The scripts cannot be run again from
this deposit, because the reads and the node-made inputs are on the node.

## 1. The chain, step by step

```
(a) hits fastq $HOME/p11_v11/hits/K562r.hits.fastq.gz md5 1a1edaa2d7d1 (2,797,988 bytes)
 $HOME/p11_v11/hits/GM12878FS.hits.fastq.gz md5 19ec8e039999 (13,810,674 bytes)
 the header of each read carries ML (and MM) tags: the m6A calls already exist here
(b) analysis/fq2ubam.py fastq -> unaligned BAM on a fake reference `chrFake` -> all/<L>.u.bam
(c) samtools sort 1.19.2 -> all/<L>.bam (command line in fiberhmm_call_records.tsv)
(d) fiberhmm-call 2.16.8 -> out/<L>.calls.bam (command line in fiberhmm_call_records.tsv)
(e) analysis/cell_fhmm_units.py calls.bam -> units cut at the in-read anchor GTCTACGGC
(f) analysis/cell_units2.py gene-part methylation of each unit
(g) analysis/build_gene_all.py gene-part strata labels (CpG >= 5, beta 0.5)
(h) analysis/cell_sig_run_all.py cells/sig_<L>.json / .npz
 analysis/cell_strata_all.py cells/strata_<L>_gene.json / .npz
 analysis/cell_nucfree.py cells/nucfree_<L>.json
(i) fig06/fig06_G-J_cell_lines.py , figS9/figS9_C-F_cell_lines_onerow.py
```

`<L>` is `K562all` or `GM12878all`. Step (a) was made by another run (see section 4); (b)-(d) are ours.
The hits fastq in (a) are named `K562r` and `GM12878FS` on the node; `fq2ubam.py` was run on them with the
output names `all/K562all.u.bam` and `all/GM12878all.u.bam`.

## 2. FiberHMM, measured from the `calls.bam` headers

| | |
|---|---|
| program / version | `fiberhmm-call` 2.16.8 (`@PG VN`) |
| command (both whole-pool BAMs) | `fiberhmm-call -i all/<L>.bam -o out/<L>.calls.bam --enzyme hia5 --seq pacbio -c 8 --primary --no-qc` |
| model | bundled, none given with `-m`. The run log names it: `fiberhmm/models/hia5_pacbio.json` for both the apply model and the recall model |
| other recorded settings (`@PG DS`) | `mode=pacbio-fiber`, `enzyme=hia5`, `recall_nucs=True`, `nuc_recall_policy=conservative`, `chimera_filter=n/a`, `dedup=off`, `daf_snp_mask=off`; `phase_nrl=181` (K562all) and `186` (GM12878all) were estimated from the data by the program |
| log line | `min_llr=5.0 min_opps=3 unify_threshold=90 uplift=1.0` |
| samtools | 1.19.2 (`@PG` of the same files) |

The command lines of all four BAMs are in `fiberhmm_call_records.tsv`.

There is no saved script for the whole-pool calls. The command exists only as the `@PG CL` field of the
two BAMs and as a shell command typed on the node. The earlier set (`K562r`, `GM12878FS`, placed reads) is
not used by these panels.

## 5. The scripts of this stage, in `analysis/`

These are byte-for-byte copies of the files that ran on the compute node, paths
included, so they read the node's folders rather than paths on a local machine.

| File | md5 | Makes |
|---|---|---|
| `fq2ubam.py` | `59b721baa51e` | `all/<L>.u.bam` from the hits fastq |
| `build_gene_all.py` | `09573f21c40a` | gene-part strata labels for the whole pool |
| `cell_sig_run_all.py` | `3888eb6ab8a9` | `cells/sig_K562all.*`, `cells/sig_GM12878all.*` |
| `cell_strata_all.py` | `485d16c453c3` | `cells/strata_K562all_gene.*`, `cells/strata_GM12878all_gene.*` |
| `cell_nucfree.py` | `69d5c315c2db` | `cells/nucfree_K562all.json`, `cells/nucfree_GM12878all.json` |
