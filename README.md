# human-5S-rDNA-array

Analysis code and result tables for the manuscript on the human 5S rDNA array.
The code builds a flank-to-flank 5S rDNA array reference for each haplotype,
places long reads on it, and measures sequence variation, CpG methylation and
chromatin accessibility one repeat unit at a time.

## How this repository is laid out

| Folder | What is in it |
|---|---|
| `analysis/` | the scripts that produce the result tables, one file per step |
| `common/` | shared libraries and the per-figure chains that several scripts import |
| `upstream/` | the compute-node stages that produced the inputs of `analysis/` |
| `tables/` | the result tables the figures read, and the Source Data workbooks |
| `checks/` | two check tables on sequencing chemistry, written by `analysis/chem_flip_2026-10-24.py` and `analysis/chem_strata_2026-10-24.py` |
| `fig01/` … `figS13/` | one manifest per figure |
| `array_acquisition_other_species/` | the cross-species array acquisition code |

**Start from the manifests.** Every `figNN/MANIFEST_*.tsv` has one row per panel
and gives the script that made the panel's data, the input tables it read, and
the md5 of each, written as a path inside this repository. That is the route from
a number in the paper to the file it came from: the Index sheet of the paper's Source
Data file lists every panel of every figure with the sheet that holds its values.

Four result tables are too large for a repository and are in the Zenodo record,
named with their sizes in `tables/README.md`.

## Running the code

To run without writing over existing tables, point every output at an empty
directory:

```
export OUT_DIR=/tmp/check SOURCEDATA_DIR=/tmp/check RDNA_OUT_DIR=/tmp/check
export DC_OUT=/tmp/check/DC_TABLE.tsv TRIM_OUT_DIR=/tmp/check
export UNITCMP_OUT=/tmp/check/FIG2_data_unitcmp.tsv
export HIFI_OUT=/tmp/check/FIG2_data_hifi.tsv
```

Run with `PYTHONDONTWRITEBYTECODE=1` (or `python -B`) so that no `.pyc` is
written next to the source. All analyses were carried out in Python with NumPy
and SciPy; 28 scripts read BAMs through `pysam`.

## Environment variables

Every root is taken from the environment, with no default. A script stops with
a message naming the variable it needs.

| Variable | What it points at |
|---|---|
| `DATA_ROOT` | the result tables (`02_tables/*.parquet`, `*.tsv`) |
| `TABLE_DIR` | the per-figure input tables, and where outputs are written |
| `V_ROOT` | consensus maps and the raw methylation calls |
| `CANON_FA` | the canonical 2,240 bp unit sequence (`canon_v2.fa`) |
| `WORK_ROOT` | the working tree holding edition stores and result tables |
| `RDNA_RESULTS` | the result tables under that working tree |
| `REGISTRY_DIR` | the registry: the tier and slot tables shared between figures |
| `RDNA_FIBERSEQ_DIR` | the Fiber-seq accessibility tree |
| `RDNA_FHMM_DIR` | the FiberHMM and transcription-factor footprint campaign |
| `DISPERSED_ROOT` | the dispersed-versus-array methylation comparison |
| `RDNA_HPRC_R2_ONT_INDEX` | the HPRC release 2 ONT sequencing index, read by three scripts in `analysis/` |
| `OUT_DIR` | where the figures are written |
| `FIG1_SCHEMATIC_PDF` | the vector schematic that panels A and B are drawn on |
| `RDNA_DATA_ROOT` | root of the deposited input data; every input is `<root>/<relative path>`, the relative path being the one the script writes |
| `RDNA_OUT_DIR` | where a script writes when `FIGOUT` / `OUTDIR` is not given; must be outside this deposit |
| `V11_CODE` | folder from which the hooks and producers import `v11_io` |
| `V11_0913_CODE` | folder `S24_NOBETA.py` copies `S24_optimal_zoning.py` from |
| `RDNA_FONT_DIR` | font folder |
| `CHROME` | headless Chrome used for the composite png and pdf |
| `P04_SRC` | folder holding `OUT_P04/` for `COSHARE_BSTRAT.py` |
| `RDNA_WORK_ROOT` | the working data root: the per-array and per-CpG tables of the 1kGP-ONT (P10) analysis, the dispersed-site tables, and the extracted region BAMs under `P10/` |
| `RDNA_OUT_ROOT` | the output root: `02_output_tables/` (the result tables the figures read) and `05_requests/Atlas_reply_2026-09-30/` (the transcript annotation) |
| `RDNA_V11_ROOT` | the v11 read-variant root: `02_tables/` with the read-unit tables and the slot QC |
| `RDNA_VERSION_ROOT` | the version root holding the consensus maps under `08_versions/v11/03_maps/` and the dispersed-site source under `05_checks/` |
| `RDNA_HUMAN_L1_PARQUET` | the file `L1_READ_UNIT.parquet` (57 MB), in the Zenodo record |
| `RDNA_POP_PANEL` | the file `pop.panel`: the 1000 Genomes population panel, a public reference |
| `RDNA_ACCESS_DIR` | the folder holding `a13_unit_state_gmm.tsv.gz` and `a01_tracks/` |
| `RDNA_S6_DIR` | the s6 call folder holding its own `a13_unit_state_gmm.tsv.gz` |
| `RDNA_V12_CANON_DIR` | the canonical layer: `02_canonical/`, `08_methylation_s5/`, `08_methylation_s6/` |
| `RDNA_RESULTMAP_DIR` | the result-map folder: `r20_*`, `s01_*` tables |
| `RDNA_A01_TRACKS_OUT_DIR` | where `a01_unit_tracks.py` writes its per-sample tracks |
| `RDNA_PREV_OUT_DIR` | an output root that a few scripts read for inputs |
| `RDNA_SVG_DIR` | where figure SVGs are copied |
| `RDNA_V11_READVAR_DIR` | the read-variant folder the chemistry-stratified checks read |
| `RDNA_METH_BLOCKS_DIR` | the methylation-block folder, including `layer_v11_fix4_HPRC/out` |
| `RDNA_ARRAY_METH_DIR` | the array-methylation folder holding `02_output_tables/` and `03_figures/` |
| `RDNA_DEPOSIT_DIR` | this repository, when a script needs to read a file from it by absolute path |
| `RDNA_FIGURE_CODE_DIR` | the figure-drawing code |
| `RDNA_SUBMISSION_FIG_DIR` | the folder holding the submitted figure files |
| `RDNA_FIG4_LEGEND_DIR` | the folder holding the Figure 4 and Supplementary Figure 4 legend files |
| `RDNA_OTHERSP_DIR` | the cross-species folder, whose `data_v2/` the acquisition scripts read |
| `RDNA_V11_FIX4_DIR` | the v11-fix4 working tree the unit-methylation run used |
| `NODE_HOME` | the home directory of the compute node, for the `upstream/` stages that ran there. |
| `NODE_WORK_DIR` | the working directory of the compute node, for the cross-species scripts in `array_acquisition_other_species/code/` that read and write their tables there |
| `RDNA_LIB_DIR` | the folder a script imports its helper modules from: `analysis/` for `a_lib.py`, and `common/fig07_figS7/` for `betabins.py`, `prof_smooth.py`, `landsc.py`, `figorder.py` and `tab_en.py`, which `make_sourcedata_fig07.py` imports |
| `RDNA_SOURCEDATA_DIR` | where `make_sourcedata_fig07.py` writes its workbook |
| `RDNA_CAPTION_DIR` | where `zzab_nucfree_bin.py` writes `TableS7_nucleosome_removed_occupancy.tsv` |
| `RDNA_CELLS_DIR` | the cell-line Fiber-seq folder for K562 and GM12878: `cell_fhmm_units.py` writes `UNIT_<label>.tsv`, `FHMM_UNITS_<label>.tsv`, `FHMM_READCOV_<label>.tsv` and `SITE_M6A_<label>.tsv` there and the other `cell_*.py` scripts read them |
| `RDNA_HEPG2_DIR` | the same folder for HepG2 |
| `RDNA_FS_OUT_DIR` | the folder of extracted reads `<sample>.hit.fq` that `tf_untreated.py` reads |
| `RDNA_GENE_QUERY_FA` | the 5S gene query sequence that `tf_untreated.py` gives to `blastn` |
| `RDNA_TF_TMP_DIR` | the scratch folder of `tf_untreated.py` |
| `P11_DIR` | the folder holding `RU_GM12878FS.tsv.gz`, `RU_GM12878.tsv.gz` and `RU_K562r.tsv.gz`, read by `cell_recompute.py` |
| `P11_TABLE_DIR` | the input tables of the Figure 6 zone and marker scripts (`analysis/fig06_*.py`) |
| `P11_TABLE_OUT_DIR` | where those scripts write |
| `P11_ASSEMBLY_DIR` | the folder holding `MAP_hapB_v4_110units.tsv` |
| `V11PIPE_DIR` | the folder holding `MAP_G8_K562r.tsv` and the read-unit tables of the Figure 6 zone scripts |
| `ARRAYMETH_ROOT` | the root whose `02_output_tables/P10/` is read by `p10_19c_steps_unitlevel.py` and written by `p10_20_array_bounds.py` |
| `S42SET29_SITES_TSV` | a table whose first column lists the 29 sites that `s42_dedup_tables_2026-11-01.py` reports separately |
| `S42_KEEP_NODE` | the node label of the NA21133 run that `s42_dedup_tables_2026-11-01.py` keeps when the sample is counted once |
| `RDNA_S` | the working folder of the reference-set steps `refset_step0_*.py` and `refset_step3_assign_slots.py`: `units_all.tsv`, `seqidx_all.tsv` and `hap_clusters.tsv` are read from it, and `HAP_READ_GLOBAL.tsv`, `HAP_READ_FINAL.tsv` and `hap_seq.tsv` are written to it |
| `RDNA_ASM_COMPARE` | the folder holding `asm_units_all.tsv` |
| `RDNA_UNIFIED_REFMAP` | the folder holding the `refmap_hap_ip-*.tsv` tables |
| `RDNA_V11` | the v11 tree: `02_canonical/`, `03_maps/`, `06_units_final/`, `08_methylation/raw/` and `12_flank/` |
| `RDNA_V11_FIX1` | the v11-fix1 tree holding `03_maps/g8canonv11f1/CONSENSUS_SUM_G8_full.tsv` |
| `RDNA_V11_FIX4` | the v11-fix4 tree; `refset_step6_registry_build.py` stops without it |
| `RDNA_READVAR` | the read-variant tree holding `02_tables/SLOT_QC.tsv` and `09_inventory/hprc_r2_phasing_meta.csv` |
| `RDNA_WORK_DIR` | the working root whose `cmp/g8place/fix1/` is read and written by `build_units_f1.py` and `build_meth_f1.py` |
| `X12_WINDOWS` | `nine` or `zone`: the set of variant windows `figx12_zone.py` uses |
| `X12_OUT` | where `figx12_zone.py` and `x12_poptoken_perm.py` write |
| `X12_PUB` | the folder holding `X12_POP_VAR_DOWNSAMPLE.tsv`, read by `x12_poptoken_perm.py`, and `X12_TYPE_GROUPS.tsv`, `X12_SIG_TOKENS.tsv` and `X12_HAP_GROUPS.tsv`, which `figx12_zone.py` compares its output with when `X12_WINDOWS` is `nine` |
| `STAGE_OUT` | the folder to which `hooks/hook_stage/sitecustomize.py` redirects the output of `v11_io` |

The scripts below take the inputs of one run from the environment. They read
`os.environ` directly, so a missing variable stops the script with a `KeyError`
that names it.

| Variable | Read by | What it is |
|---|---|---|
| `HSN` | `consmapG8v4_fixdedup.py`, `refset_step4_5_unit_map_and_variants.py`, `rjf_repro.py`, `jxplace.py`, `slotvote.py`, `slotvote_qc.py` | the `hap_seq_node` list |
| `HSNF` | `fsplace_v12.py`, `pbfsplace_v11.py` | the `hap_seq_node` list |
| `CMPD` | the eight scripts of the two rows above | the folder of array `.fa` files |
| `WDIR` | `consmapG8v4_fixdedup.py`, `refset_step4_5_unit_map_and_variants.py`, `aoplace_v3fix.py`, `aoplace_v5fix.py`, `aoplace_v8fix.py` | the work folder |
| `MAPF` | `aoplace_v8fix.py`, `rjf_repro.py`, `fsplace_v12.py`, `jxplace.py`, `pbfsplace_v11.py` | the `CONSENSUS_MAP` table |
| `SUMF` | `fsplace_v12.py`, `jxplace.py`, `pbfsplace_v11.py` | the `CONSENSUS_SUM` table, whose `revcomp` column gives the orientation of each array |
| `RCFF`, `RCFCOL` | `rjf_repro.py` | the orientation table, and the index of its column that holds `Y` for a reverse-complemented array |
| `G8MAP` | `rjf_repro.py` | the consensus map that supplies the G8 unit boundaries |
| `RJW` | `rjf_repro.py` | the work folder |
| `READS` | `jxplace.py`, `slotvote.py`, `slotvote_qc.py` | the folder of extracted reads `<sample>.hit.fq` |
| `BOUND` | `slotvote.py`, `slotvote_qc.py` | the folder of `<sample>.bound.tsv` |
| `READD` | `pbfsplace_v11.py` | the folder of `<label>.hit.fq` and `<label>.bound.tsv` |
| `SAMD` | `pbfsplace_v11.py` | the folder its SAM files are written to |
| `SRCD`, `FTW` | `ftprep_v12b.py` | the folder of `<sample>.AOJV_PB.sam`, and the work folder |
| `CFG`, `LABEL`, `OUT` | `fsmeth_v12c.py` | `AOJV_FS` or `AOJV_PB`, selecting the SAM files that are read; `UCSC_ONT_FIBERSEQ` or `PACBIO_FIBERSEQ`; the output folder |
| `ASM`, `OUT` | `flank_reasm.py` | the folder of assembled contigs `<sample>.bp.hap?.p_ctg.fa`, and the output folder |
| `OUT` | `flank_cons.py`, `flank_reads.py` | the output folder |
| `REPRO_IN`, `AOREADS` | `extract_ao_reads_env.py` | the folder of accepted-read lists (`AOJV_GIAB.tsv`, `AOJV8.tsv`), and the output folder |
| `EXP6` | `run_fix2_ao.sh` | the md5 that `MAP_G6C2_f2u.tsv` must have; the script stops when it differs |

## Environment

`requirements.txt` lists the Python packages the code imports. `environment.yml`
is a conda environment that adds the aligners.

The releases of the external tools are in Supplementary Table 1 of the paper,
and `environment.yml` repeats the compute-node ones. The code runs on Python 3.

## Tools the code relies on

The releases are those the paper records in its Methods and in Supplementary
Table 1.

| Tool | Release | Used for | Reference |
|---|---|---|---|
| BBDuk (BBTools) | 40.02 | k-mer gate that extracts the 5S reads from each BAM | https://escholarship.org/uc/item/1h3515gn |
| hifiasm | 0.25.0-r726 (`--ont`) | the references assembled here | doi:10.1038/s41592-020-01056-5, doi:10.1038/s41586-026-10105-6 |
| hifiasm, verkko | 0.19.7-0.19.9, 2.2.1 | assemblers of the public references | doi:10.1038/s41587-023-01662-6, doi:10.1101/gr.280383.124 |
| minimap2 | 2.26-r1175, 2.28-r1209 | unit map and read placement, upstream junction scan | doi:10.1093/bioinformatics/bty191, doi:10.1093/bioinformatics/btab705 |
| SAMtools | 1.19.2, 1.9 | unit map and read placement, upstream junction scan | doi:10.1093/gigascience/giab008 |
| BLAST+ | 2.12.0+, 2.16.0 | 5S gene search | doi:10.1186/1471-2105-10-421 |
| parasail | 1.3.4 | pairwise global alignment of unit sequences | doi:10.1186/s12859-016-0930-z |
| dorado | 0.6.0, 0.9.0 | basecalls and methylation calls of the R10.4.1 reads, made by the data providers | https://github.com/nanoporetech/dorado |
| guppy | 6.5.7 | basecalls and methylation calls of the R9.4.1 reads, made by the data providers | https://community.nanoporetech.com |
| jasmine | 2.3.0, 2.7.99, 26.1.3 | m6A calls of the HPRC PacBio Fiber-seq, made by the submitters | https://github.com/PacificBiosciences/jasmine |
| fibertools-rs | 0.13.0, 0.1.4, 0.3.8 | fiber segments (`ft add-nucleosomes`) and the GM12878 FIRE release (`ft predict`, `ft fire`) | doi:10.1101/gr.279095.124 |
| FiberHMM | 2.16.8 | nucleosome and MSP segments | https://github.com/fiberseq/FiberHMM |
| NumPy, SciPy | | all analyses | doi:10.1038/s41586-020-2649-2, doi:10.1038/s41592-019-0686-2 |

## Data sources

All input data are public. The Data availability statement and the Methods of
the paper describe them in full; Supplementary Tables 1, S6 and S7 give the
per-sample detail.

| Data | Accession or location |
|---|---|
| HPRC assemblies | release 2 index `assemblies_release2_v1.0.index.csv` |
| HPRC ONT reads | BioProject PRJNA701308; `s3://human-pangenomics/working/HPRC/` |
| 1kGP ONT release, with its GIAB HG003 assembly | `s3://1000g-ont/PROCESSED_DATA/` |
| GIAB ONT reads | `s3://ont-open-data/giab_2023.05/` |
| GIAB HG001 PacBio HiFi reads | `s3://giab/data/NA12878/PacBio_CCS_15kb_20kb_chemistry2/` |
| GM12878 ONT reads | ENA ERR17530411 (BioProject PRJEB114757); `s3://ont-open-data/giab_2025.01/analysis/wf-human-variation/sup/HG001/` |
| K562 ONT reads | BioProject PRJNA1371648 |
| Tumor and matched normal, five cell-line pairs | BioProject PRJNA1086849 |
| Tumor and matched normal, COLO829 | `s3://ont-open-data/colo829_2024.03/` |
| Tumor and matched normal, HG008 | `ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/` |
| HPRC Fiber-seq | `s3://human-pangenomics/submissions/`, submissions `UCSC_HPRC_ONT_YEAR1_FIBERSEQ` and `WUSTL_HPRC_PacBio_FiberSeq` |
| K562 Fiber-seq | ENA PRJEB121523, run ERR17608832 |
| GM12878 Fiber-seq | `s3://stergachis-manuscript-data/2024/Vollger_et_al/FIRE/GM12878/fire/GM12878.fire.bam` |
| ENCODE ChIP-seq, ATAC-seq and whole-genome sequence of K562 and GM12878 | file accessions in Supplementary Table 6 |
| GM12878 whole-genome denominator | 1000 Genomes run ERR3239334 |
| GM12878 GRO-cap | GEO GSE60453 |
| Annotation sequences | EMBL X12811; RefSeq NR_002716, NM_001306068, NR_147172, NR_147173 |
| Cross-species ONT reads and reference assemblies | per run in Supplementary Table 7 |

## Paths of the `upstream/` stages

The node stages build their paths under the node's home directory with
`os.path.expanduser("~/…")`. `expanduser` reads `HOME`, so pointing `HOME` at a
tree that reproduces that layout relocates them. A node script that takes a
folder or a table from a variable of its own (`WDIR`, `CMPD`, `HSN`, `MAPF`,
`SUMF`, `READS`, `OUT` and the like) reads it with no default.

A few node shell scripts under `array_acquisition_other_species/node_scripts/`
read their region table (`zf_regions_parts.tsv`) from the current directory.

## Operating system

The scripts were run on Linux and assume a POSIX system. Some write scratch
files to `/tmp`, and the shell scripts are written for `sh`/`bash`.

## Terms used in the code

Names used in the code and the terms the paper uses for them.

| In the code | What it means | In the paper |
|---|---|---|
| `slot` | the array position assigned to one haplotype of one person | haplotype |
| `unit` | one repeat of the array | repeat unit |
| `depth` | a unit's number counted from a junction | depth |
| `zone` | a stretch of the array set off by methylation or by variant content | zone |
| `canon`, `canonical unit` | `unit_canon.fa`, the 2,240 bp reference unit cut from CHM13 | the canonical unit |
| `tier` | the cohort subsets H0 to H3 | cohort tier (H0-H3) |
| `ruler` | the definition a quantity is measured with: its cut-off, denominator and aggregation unit | the ruler / the definition |
| `hook` | a `sitecustomize.py` that replaces `v11_io.load` for one run | not used in the paper |
| `producer` | the script that makes a result table | the script that made the panel's data |
| `v11`, `v12`, `v13`, `fix1` to `fix4` | internal generations of the pipeline | not used in the paper |
| `P04`, `P08` | internal stage numbers of the variant and methylation analyses | not used in the paper |
| `P10` | the 1kGP-ONT analysis | not named as such; its results are the 1kGP comparisons |
| `P11` | the cell-line marker analysis | not named as such; its results are the ENCODE zone measurements |
| `NEWZONE` | table names that carry the zoning the paper reports | zone |
| `readvar`, `edition`, `campaign`, `atlas`, `handover`, `S24` | working names that appear only in comments and file names | not used in the paper |

## Scripts that produce no figure panel

| Script | What it does |
|---|---|
| `analysis/supp_tables_unit_vs_gene_beta_2026-11-01.py` | produces Supplementary Tables 2 and 3 |
| `analysis/n_vs_chem_2026-10-25.py` | a check of whether the 1kGP curve depends on mixed sequencing chemistry or on the smaller number of units |
| `analysis/profile_by_chem_2026-10-25.py` | the positional methylation profile computed separately per sequencing chemistry |

## Randomness

Every randomized procedure is seeded in the source and the seed is given where
it is used.

## Licence and citation

MIT, see `LICENSE`. `CITATION.cff` gives the citation for this repository; the
archived release has its own DOI, which the Data Availability statement of the
manuscript names.

## AI assistance

The analysis code, scripts, and documentation in this repository were
developed with the assistance of Claude (Anthropic), used as a coding and
editing assistant (via Claude Code). The author designed the analyses,
reviewed the code and outputs, and takes full responsibility for the
content of this repository.
