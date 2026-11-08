# Upstream code behind Figure 2 and Supplementary Figure 2

The compute-node and local code behind the two figures:
the code that made the data the panels draw, followed back to the public data. The code that runs inside
this deposit (`common/fig02_figS2/tables/`) starts from intermediate tables; the scripts here are the steps
before it. `MANIFEST_upstream_fig02_figS2.tsv` lists every file with the chain and step it belongs to, where
it was found, the md5 of the file found and of the file deposited, what it does and the line that shows it.
The full account (chains U1–U10, every panel, the public data and how to cite them) is kept with the working
files (`00_Fig2_FigS2_code_chain_raw_data_references_2026-10-23.md`, "code chain, raw data and references") and summarised below.

## Chains

| folder | what it made |
|---|---|
| `U0_reference_sequences` | `mkref_aws.sh`: the canonical unit `unit_canon.fa` (reverse complement of CHM13 chr1:227,887,698-227,889,937, 2,240 bp), the 5S gene query `5S_gene_query.fa` (unit positions 2122-2240), the discriminating segment, the tandem reference and the 50 kb flank probes, all cut from CHM13 v2.0 |
| `U1e_reference_slot_list` | `hap_seq_node_v8_chain.sh`: the reference list `hap_seq_node_v8.tsv` (616 slots) from the slot assignment `seqassign.py` (deposited as `analysis/refset_step3_assign_slots.py`) through v7 to v8 |
| `U5_cohort_tiers` | `tiers_H0H3_build.sh`: the cohort tier list `TIERS_H0H3_SLOTS_2026-10-05.tsv` (H0 604, H1 557, H2 510, H3 490 slots) |
| `U1a_5S_read_extraction` | 5S reads from each ONT submission: BBDuk (k = 25, hdist = 2, `mkh=2 rcomp=t`) on a discriminating 5S segment, MM/ML tags kept; array call from a minimap2 map-ont alignment to a tandem reference; the same for GIAB ONT; the 1kGP-ONT array region cut from the CHM13-aligned R10 BAMs (`kgp_ont/manifest.tsv` lists the BAMs read) |
| `U1b_HPRC_array_cut` | 400 array references cut from HPRC release 2 assemblies with CHM13 proximal and distal 50 kb probes (minimap2 `-x asm20`) |
| `U1c_our_hifiasm` | 197 array references assembled here with `hifiasm --ont` 0.25.0-r726 from the reads of U1a and cut with the same probes; `tier1ref.py`, the 12 single spanning reads |
| `U1d_reference_repair` | repair of reference windows (read fragments, in-house backbone polishing on the racon principle, splice); `v11` made `cmp4d_v10`, `v11-fix1` replaced two haplotypes (`cmp4d_v11_fix1`); `readmap_repo` holds the read-consensus reference scripts of the earlier pipeline repository |
| `U2_unit_map_variant_calling` | **the per-base variant calls of Figure 2A and of every panel that reads units**: units cut between consecutive 5S genes found by blastn, each aligned with minimap2 `-x asm20 -c --cs=long` to the canonical unit (CHM13 chr1:227,887,698-227,889,937, reverse complement) with 200 bp of context, the `cs` string written as tokens; `rung8v10.sh` made the map read by the read-level chain (U7), `run_map_fix1.sh` the map read through `v11_io.py` |
| `U3_read_placement_units` | ONT reads placed on the references (ARRAY_ONLY and junction), units cut, variants re-called, merged by data source into `READMAP_V11F2_BYSOURCE.tsv`; `stage_all4.py` is the assembly-stage rule that `tables/stage_v11.py` copies |
| `U4_unit_methylation` | per read × unit methylation parsed from MM/ML (ML ≥ 128), merged as Σmeth/Σcalls; fix1, fix2 and fix4 versions of `READMAP_METH` |
| `U7_read_level_tables` | read-level tables behind Figure 2 M, N and Supplementary I (varlist from the U2 map `g8canonv10`) |
| `U8_flank` | `flank_final2.tar`: 50 kb flanks from HPRC contigs, junction-read consensus and `hifiasm --ont` re-assembly |
| `U9_trio_references` | the trio map of Supplementary Figure S2A (children from HPRC release 2; parents assembled here from 1kGP-ONT reads) |
| `U10_cell_line_arrays` | the K562 and GM12878 arrays of Supplementary Figure S2O (K562: ERR17608832 HiFi + PRJNA1371648 ONT, `hifiasm --ul`; GM12878: GIAB HG001 HiFi, `hifiasm`) and the distal 20 units |

## Paths

Code that ran on the compute nodes is deposited as it ran, so it refers to the
node's own home folder. Scripts that ran on a local machine take every data root
from an environment variable with **no default**: `RDNA_DATA_ROOT` for the data
and `RDNA_WORK_DIR` for scratch, and a script stops with `KeyError` if one is
unset. No computing line differs from the version that produced the tables.

`U1a_5S_read_extraction/kgp_ont/manifest.tsv` is a data list, the 1kGP-ONT BAMs
that were read, not an archive manifest; it keeps the name the scripts read.

Nothing here runs without the node data. These files are deposited so that the
procedure can be read and checked.