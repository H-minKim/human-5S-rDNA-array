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
| `checks/` | `rulers.py`, which several scripts import so that a measurement and the thing it measures use the same ruler, and the two tables that read it |
| `fig01/` … `figS13/` | one manifest per figure, and nothing else |
| `array_acquisition_other_species/` | the cross-species array acquisition code |

**Start from the manifests.** Every `figNN/MANIFEST_*.tsv` has one row per panel
and gives the script that made the panel's data, the input tables it read, and
the md5 of each, written as a path inside this repository. That is the route from
a number in the paper to the file it came from: Supplementary Table S11 lists
every panel of every figure with the manifest that holds its row.

## What is here and what is not

This repository holds **the code that produced the data**, and the result tables
it produced. The code that takes a table and draws a figure is not here, so a
manifest row may name a drawing script that this repository does not contain;
the table that panel was drawn from is here, with its md5 and its row count.

Four result tables are too large for a repository and are in the Zenodo record
instead, named with their sizes in `tables/README.md`. Image files
(`*.png`, `*.pdf`, `*.svg`) are excluded by `.gitignore`, as the figures
themselves are.

## Running this code: read this before the first run

**Every default output path points at the live working tree, not at a scratch
folder.** That is deliberate, because it is what makes the scripts run unchanged
on the machine the analysis was done on, but it means that running this deposit
as it stands overwrites result tables in place. Nothing warns you first.

To exercise the code without touching anything, point every output at an empty
directory:

```
export OUT_DIR=/tmp/check SOURCEDATA_DIR=/tmp/check RDNA_OUT_DIR=/tmp/check
export DC_OUT=/tmp/check/DC_TABLE.tsv TRIM_OUT_DIR=/tmp/check
export UNITCMP_OUT=/tmp/check/FIG2_data_unitcmp.tsv
export HIFI_OUT=/tmp/check/FIG2_data_hifi.tsv
```

Run with `PYTHONDONTWRITEBYTECODE=1` (or `python -B`) so that no `.pyc` is
written next to the source. All analyses were carried out in Python with NumPy
and SciPy; three scripts read BAMs through `pysam`.

## Environment variables

Every root is taken from the environment, with no default, so that no absolute
path of the machines that ran the analysis is baked into the code. A script
stops with a message naming the variable it needs.

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
| `RDNA_CELLS_ALLPOOL_DIR` | the pooled cell-line campaign |
| `ENCODE_UNIT_ROOT` | the ENCODE unit-placement tree |
| `DISPERSED_ROOT` | the dispersed-versus-array methylation comparison |
| `ANALYSIS_ROOT` | the analysis disk, one dated folder per campaign |
| `RDNA_HPRC_R2_ONT_INDEX` | the HPRC release 2 ONT sequencing index, read by three scripts in `analysis/` |
| `P10_DIR` | the 1kGP P10 frame |
| `OUT_DIR` | where the figures are written |
| `FIG1_SCHEMATIC_PDF` | the vector schematic that panels A and B are drawn on |
| `RDNA_DATA_ROOT` | root of the deposited input data; every input is `<root>/<relative path>` as listed under "Raw data read by the producers" |
| `RDNA_OUT_DIR` | where a script writes when `FIGOUT` / `OUTDIR` is not given; must be outside this deposit |
| `V11_CODE` | folder from which the hooks and producers import `v11_io` |
| `V11_0913_CODE` | folder `S24_NOBETA.py` copies `S24_optimal_zoning.py` from |
| `RDNA_FONT_DIR` | font folder (Korean figures only) |
| `CHROME` | headless Chrome used for the composite png and pdf |
| `P04_SRC` | folder holding `OUT_P04/` for `COSHARE_BSTRAT.py` |
| `RDNA_WORK_ROOT` | the working data root: the per-array and per-CpG tables of the 1kGP-ONT (P10) analysis, the dispersed-site tables, and the extracted region BAMs under `P10/` |
| `RDNA_OUT_ROOT` | the output root: `02_output_tables/` (the result tables the figures read) and `05_requests/Atlas_reply_2026-09-30/` (the transcript annotation) |
| `RDNA_V11_ROOT` | the v11 read-variant root: `02_tables/` with the read-unit tables and the slot QC |
| `RDNA_VERSION_ROOT` | the version root holding the consensus maps under `08_versions/v11/03_maps/` and the dispersed-site source under `05_checks/` |
| `RDNA_HUMAN_L1_PARQUET` | the file `L1_READ_UNIT.parquet` (57 MB). **This archive cannot rebuild it** - it is an output of the `v11_readvar` analysis. The same variable name is used by the cross-species part of this work, so one copy serves both |
| `RDNA_POP_PANEL` | the file `pop.panel`: the 1000 Genomes population panel, a public reference. It points at the file and not at its folder, because nothing else in that folder is read |
| `RDNA_ACCESS_DIR` | the folder holding `a13_unit_state_gmm.tsv.gz` and `a01_tracks/` |
| `RDNA_S6_DIR` | the s6 call folder holding its own `a13_unit_state_gmm.tsv.gz` |
| `RDNA_V12_CANON_DIR` | the canonical layer: `02_canonical/`, `08_methylation_s5/`, `08_methylation_s6/` |
| `RDNA_RESULTMAP_DIR` | the result-map folder: `r20_*`, `s01_*` tables |
| `RDNA_A01_TRACKS_OUT_DIR` | where `a01_unit_tracks.py` writes its per-sample tracks |
| `RDNA_PREV_OUT_DIR` | an earlier output root that a few scripts still read for inputs |
| `RDNA_SVG_DIR` | where figure SVGs are copied |
| `RDNA_V11_READVAR_DIR` | the read-variant folder the chemistry-stratified checks read |
| `RDNA_METH_BLOCKS_DIR` | the methylation-block folder, including `layer_v11_fix4_HPRC/out` |
| `RDNA_ARRAY_METH_DIR` | the array-methylation folder holding `02_output_tables/` and `03_figures/` |
| `RDNA_DEPOSIT_DIR` | this repository, when a script needs to read a file from it by absolute path |
| `RDNA_FIGURE_CODE_DIR` | the figure-drawing code, which is not part of this repository |
| `RDNA_SUBMISSION_FIG_DIR` | the folder holding the submitted figure files |
| `RDNA_FIG4_LEGEND_DIR` | the folder holding the Figure 4 and Supplementary Figure S4 legend files |
| `RDNA_OTHERSP_DIR` | the cross-species folder, whose `data_v2/` the acquisition scripts read |
| `RDNA_V11_FIX4_DIR` | the v11-fix4 working tree the unit-methylation run used |

## Randomness

Every randomized procedure is seeded in the source and the seed is given where
it is used. Where a quantity was computed in more than one way, both values are
reported rather than reconciled.

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
