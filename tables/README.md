# tables/, the result tables the figure code reads

120 tables, 45.8 MB. Each one is the exact file a figure was drawn from: it was
selected by matching the `input_table_md5` recorded in the `MANIFEST_*.tsv` files,
not by file name, because several of these names exist in more than one version
folder in the working tree with different content.

Point the table roots at this folder when re-running the code, for example
`TABLE_DIR`, `DATA_ROOT`, `P11_TABLE_DIR` and `REGISTRY_DIR`
(the full list is in the repository README, section *Environment variables*).

**4 tables are not here because of their size**; they are in the Zenodo record
that archives this repository:

| MB | table |
|--:|---|
| 68.7 | `P10_FRAME_ALL.tsv` |
| 65.2 | `P10_CPG_ADD3.tsv.gz` |
| 57.0 | `L1_READ_UNIT.parquet` |
| 54.0 | `P10_CPG_ALL.tsv.gz` |

Of those four, `L1_READ_UNIT.parquet` is the one that **this repository cannot
rebuild**: it is an output of the `v11_readvar` analysis, which runs on the
compute-node side from the placed SAM files of `upstream/v11_node/`, and those
SAM files are intermediates of the public reads and are not deposited. Take the
file from the Zenodo record and point `RDNA_HUMAN_L1_PARQUET` at it. The other
three are rebuildable from what is here; they are out only for their size.

**9 tables carry English names.** The content is
unchanged, so the md5 recorded in the manifests still applies:

| name here | name in the working tree |
|---|---|
| `GM12878_interaction_NEWZONE_GMWGS_2026-10-14.tsv` | `GM12878_interaction_NEWZONE_GMWGS_2026-10-14.tsv` |
| `K562_interaction_NEWZONE_withBDP1BRF1_2026-10-14.tsv` | `K562_interaction_NEWZONE_with_BDP1BRF1_2026-10-14.tsv` |
| `K562_zoneratio_Q20_NEWZONE_withBDP1BRF1_2026-10-14.tsv` | `K562_Q20_zone_ratio_NEWZONE_with_BDP1BRF1_2026-10-14.tsv` |
| `MAP_hapB_v4_110units.tsv` | `MAP_hapB_v4_110units.tsv` |
| `P11_percellline_NEWZONE_GMWGS_bootstrap_2026-10-14.tsv` | `P11_by_cell_line_summary_NEWZONE_GMWGS_bootstrap_2026-10-14.tsv` |
| `TableS1_Hia5_untreated_control_by_position.tsv` | `TableS1_Hia5_untreated_control_by_position.tsv` |
| `TableS3_alternative_explanations_predicted_vs_observed.tsv` | `TableS3_alternative_explanation_predicted_vs_observed.tsv` |
| `TableS7_nucleosome_removed_occupancy.tsv` | `TableS7_nucleosome_removed_occupancy.tsv` |
| `lowmeth_changepoint_4haplotypes_2026-10-12.tsv` | `low_methyl_segment_changepoint_4haplotypes_2026-10-12.tsv` |

| table | md5 | bytes |
|---|---|--:|
| `1kGP_3202_ped_population_2026-09-23.txt` | `bac176a443f9` | 97483 |
| `ALIGNED_POS_RATE.tsv` | `a0f3a2ce577d` | 1151660 |
| `ARRAYPROF_noHG008T_2026-10-16.tsv` | `a80b86bc5741` | 9461 |
| `BOUNDARY_WIDE_H2_508.tsv` | `963331f95ce8` | 882317 |
| `CONSENSUS_MAP_G8_full.tsv` | `34dc9870ee54` | 8864831 |
| `CONSENSUS_SUM_G8_full.tsv` | `4c8252462824` | 21186 |
| `COSHARE_PERHAP.tsv` | `398d1259d27a` | 844653 |
| `COV_AFTER_FILTER.tsv` | `55cec51e3838` | 205808 |
| `DC_TABLE.tsv` | `6029cf9c721e` | 21483 |
| `DISP_PROFILES.tsv` | `0599c32be6b0` | 32869 |
| `DISP_SITES_TABLE.tsv` | `ae1bf0bcc790` | 11801 |
| `EXAMPLE_PROFILE.tsv` | `ddf56e9d56c0` | 3675 |
| `F8X2_pairs.tsv` | `2ce60ee20897` | 621 |
| `FIG1_2_MID.tsv` | `a7638c89acf3` | 678 |
| `FIG1_2_UNITS.tsv` | `5c98df7d82d3` | 1949 |
| `FIG1_V5_A_bp_HC.tsv` | `38e6e4aaa87f` | 29546 |
| `FIG2_data_hifi.tsv` | `1ff07c00c41a` | 624 |
| `FIG2_data_unitcmp.tsv` | `978bda6f0176` | 6666 |
| `FIG3_A_UNITPROFILE_ZONE.tsv` | `a315cdac6f4a` | 31581 |
| `FIG3_B_POSBETA.tsv` | `e6b0b57f9d6b` | 11495 |
| `FIG3_C_LOWCOUNT.tsv` | `580823d3003c` | 699 |
| `FIG3_G_TYPEDIST.tsv` | `77f00dfd554a` | 17489 |
| `FIG3_SUPP_A_SEQ1K.tsv` | `072cb3bc1c42` | 11004 |
| `FIG5V2_PAIRZONE5.tsv` | `d2a865760e06` | 3706 |
| `FIG8_PAIRS.tsv` | `355fb486a09f` | 1945 |
| `FIGX25_DENSITY.tsv` | `92b829838fe4` | 27495 |
| `FIGX25_PROFILE.tsv` | `1b027fc9e331` | 4624 |
| `FIG_data_table.tsv` | `2c0c009285c9` | 1233 |
| `FILTER_DECOMPOSE.tsv` | `fdfa9107be64` | 241 |
| `FILTER_DECOMPOSE_BDP1_2026-11-01.tsv` | `4d24ab436b8a` | 116 |
| `FILTER_DECOMPOSE_GMATAC_2026-11-01.tsv` | `13f08b85afcd` | 76 |
| `GM12878_interaction_NEWZONE_GMWGS_2026-10-14.tsv` | `b2355947b209` | 139 |
| `GM_QUALGATE_2026-11-01.tsv` | `0a596a6b423b` | 32158 |
| `INTEGRATED_13SITES.tsv` | `c4882945692a` | 1634 |
| `K562_QUALGATE_MIX.tsv` | `1414cf754e97` | 37510 |
| `K562_QUALGATE_MIX_BDP1BRF1.tsv` | `c627d845531e` | 9931 |
| `K562_QUALGATE_NEWZONE_2026-10-13.tsv` | `c8a0b6b322da` | 22756 |
| `K562_interaction_NEWZONE_withBDP1BRF1_2026-10-14.tsv` | `270dc585ce15` | 411 |
| `K562_zoneratio_Q20_NEWZONE_withBDP1BRF1_2026-10-14.tsv` | `79afcf0bb43e` | 725 |
| `L2_HAP_UNIT.parquet` | `2d0217ca81b8` | 9477685 |
| `MAP_G8_GM12878FS.tsv` | `c5fa06b45060` | 9205700 |
| `MAP_G8_K562r.tsv` | `a91b74d6df1b` | 9201404 |
| `MAP_hapB_v4_110units.tsv` | `68f5caa084c2` | 18897 |
| `METH_TERM_units.tsv` | `fad6e6a66b4a` | 232924 |
| `MODEL_TEST_ALL6.tsv` | `ae81acc0f287` | 5368 |
| `NORM_VS_CANCER_SAMPLE.tsv` | `d889dbccca0f` | 2297 |
| `P08_X2_zone_region.tsv` | `e6900744f779` | 3273 |
| `P10_ARRAYS.tsv` | `e6eb7058310e` | 32072 |
| `P10_CPG_SEGS_ADD3.tsv` | `a95f2e09f8dd` | 9186 |
| `P10_DXZ4_POLARITY_FINAL.tsv` | `a6a120f498ee` | 1689 |
| `P10_FRAME_SEGS_ALL.tsv` | `4108ba22cc74` | 15631 |
| `P10_ONOFF_TESTS.tsv` | `74af1e5ed562` | 869 |
| `P10_UNIT_STEPS_ABS.tsv` | `f95b9e89a86c` | 248 |
| `P11_GM12878_MARKERS_NEWZONE_2026-10-13.tsv` | `4fec23b0c0a3` | 3367 |
| `P11_K562_MARKERS_NEWZONE_2026-10-13.tsv` | `7a7ddea532af` | 2875 |
| `P11_percellline_NEWZONE_GMWGS_bootstrap_2026-10-14.tsv` | `a04e49a9d8e7` | 2623 |
| `POP_COHORT.tsv` | `c68b2063eb7a` | 51974 |
| `PROB_A4_DIST_SITE_TONE.tsv` | `01770c17a2fc` | 582858 |
| `PROB_A4_SITE.tsv` | `68ec37019501` | 5145 |
| `PROB_A4_SITE_TONE.tsv` | `0a43fa70958b` | 21873 |
| `Q1_results.tsv` | `5da4978927ef` | 1676 |
| `Q2b_results.tsv` | `a2dc296981e3` | 7311 |
| `REGCHAR_axis.tsv` | `d2da5973f84a` | 46409 |
| `REGCHAR_relcount.tsv` | `97d494a4469b` | 2955 |
| `RU_GM12878FS.tsv.gz` | `1fffbcea765f` | 1709527 |
| `RU_K562r.tsv.gz` | `f8b13b7421bb` | 380352 |
| `RU_K562u.tsv.gz` | `5e0799454a25` | 397158 |
| `S14_length.tsv` | `1de63221dd12` | 484 |
| `S18_changepoint.tsv` | `fb7c0859f602` | 107612 |
| `S1_per_position_profile.tsv` | `98996f44644b` | 169950 |
| `S24_gain_curve.tsv` | `5adce291cd1f` | 20837 |
| `S24_optimal_zoning.tsv` | `bbe1c8dab70d` | 1347 |
| `S26b_boundary_drivers.tsv` | `dfc5b7df3b36` | 1451 |
| `S26b_position_profile.tsv` | `097ff70985d7` | 10299 |
| `STAGE_ALL_V11F2.tsv` | `16c79bac4378` | 59729 |
| `SUBUNIT_BIMODAL2.tsv` | `24c5cc0b9ea9` | 14980 |
| `SUPP3_2_POSSD.tsv` | `5b1c646147a2` | 2108 |
| `SUPP3_B_HAP3.tsv` | `4650d5dd5fa0` | 267349 |
| `T1_pair_summary.tsv` | `e4f8375e8199` | 15115 |
| `TIERS_H0H3_SLOTS_2026-10-05.tsv` | `64725a33454f` | 94741 |
| `TRANSCRIBED_FOR_PANEL_A.tsv` | `0eab96ae032d` | 1806 |
| `TableS1_Hia5_untreated_control_by_position.tsv` | `07b25442b25a` | 943 |
| `TableS3_alternative_explanations_predicted_vs_observed.tsv` | `7350a4eea9a7` | 1878 |
| `TableS7_nucleosome_removed_occupancy.tsv` | `4644d7f1856d` | 775 |
| `V11F4_DATASET_hap.tsv` | `b6c5011279b5` | 143783 |
| `V11F4_DATASET_hap_flags.tsv` | `d1a8f1484fb9` | 177229 |
| `V11F4_DATASET_hap_flags_v2.tsv` | `24316816ab1d` | 179719 |
| `VTSTRAT_cpg_profile.tsv` | `ddf3d6a8c495` | 7341 |
| `VTZONE4_byzone.tsv` | `9bb1ba01cb84` | 16031 |
| `VTZONE_zones.tsv` | `217a2aa6a5d2` | 1055 |
| `W15_minority_group_gaps.tsv` | `aa4d816b01ac` | 632 |
| `W2_marker_by_depth.tsv` | `1b9e8079cfc8` | 693 |
| `X12_HAP_GROUPS.tsv` | `3a041c252d26` | 14649 |
| `X12_POP_VAR_DOWNSAMPLE.tsv` | `f3713ec0043b` | 5025 |
| `X12_POP_VAR_PERM.tsv` | `3b5bd853de47` | 322 |
| `X12_SIG_TOKENS.tsv` | `48580f870520` | 28420 |
| `X12_TYPE_GROUPS.tsv` | `6b85e0089e01` | 2527 |
| `X12_WINDOWS.tsv` | `849da655f1dd` | 24 |
| `X1_population_control.tsv` | `4aa0f1bfdbca` | 1553 |
| `X4_lag2_and_patch.tsv` | `1721d7dff79b` | 1209 |
| `X5_lag2_noncircular.tsv` | `4ce156fb8f62` | 1056 |
| `X7_homogenisation_evidence.tsv` | `6b902e916555` | 896 |
| `X9_patch_floor.tsv` | `2f5368ed3d30` | 610 |
| `XA_decay_constant_robustness.tsv` | `e51df126bb52` | 591 |
| `Y4_cellline_summary.tsv` | `4b3233519253` | 399 |
| `Y5_terminal_defining_variants.tsv` | `b013af3755fa` | 3210 |
| `ZB_neighbour_divergence.tsv` | `be24394880b3` | 398 |
| `ZM_population_carrier.tsv` | `f9d80cbf0aad` | 646 |
| `ZM_terminal_similarity.tsv` | `5398b28a3f9c` | 551 |
| `ZN_depthset_census.tsv` | `b3343a3840c1` | 2094 |
| `ZP_gap_shape.tsv` | `17ddf3048255` | 320 |
| `cohort_tiers_Fig1_2026-10-08.tsv` | `f3e56d90ad6c` | 728 |
| `fig_VTORG_single_region_type_composition_source_values.tsv` | `64ffa4f303e1` | 27296 |
| `fig_VTPD_region_type_composition_source_values.tsv` | `914783b86972` | 32514 |
| `hap_seq_node_v8.tsv` | `bc11cdc3aec7` | 25342 |
| `hor_anchor_per_read_unit.tsv` | `7ac0f3445d31` | 325340 |
| `lowmeth_changepoint_4haplotypes_2026-10-12.tsv` | `af278c7324de` | 595 |
| `peaks_both_v2.tsv` | `6f3ab4d07673` | 4249 |
| `privvar_per_read_unit.tsv` | `e6daf9ef2296` | 6098 |
| `unit_marker_meth.tsv` | `52fcf2405a86` | 172689 |

**74 more files are in this folder**, 44.0 MB, which brings it to 194 files and
89.8 MB: the 16 Source Data workbooks (`.xlsx`), and 58 data files that the scripts
or the manifests of this repository read or name.

| file | md5 | bytes |
|---|---|--:|
| `ATH_CONTEXT_PERARRAY.tsv` | `a4e9ff17195a` | 9543 |
| `ATH_READUNIT_CONTEXT.tsv.gz` | `0c38f2edab05` | 5672358 |
| `CENTROMERE.tsv` | `d852cada38d4` | 3523 |
| `CHROM_NAMES.tsv` | `073216067cba` | 1165 |
| `FIG5V2_ZONE5.tsv` | `7a6f851682cf` | 237 |
| `FIGDATA_A1_chrom_map.tsv` | `b88266084ec4` | 9750 |
| `FIGDATA_SUPP4_ATH_CONTEXT.tsv` | `288d77754a80` | 1901 |
| `FIGX19_E_HAP_cut0508.tsv` | `5029d62eb4a5` | 426 |
| `FIGX19_E_PAIR.tsv` | `2c44ea4b3f14` | 358 |
| `FIGX19_F_SEG.tsv` | `ea37d896bbd8` | 1132 |
| `FIGX25_BIMODAL.tsv` | `aea480f5560d` | 418 |
| `Fig2_SourceData.xlsx` | `4baaf42e4e07` | 161648 |
| `Fig3_main_SourceData.xlsx` | `b468e2c9fded` | 173998 |
| `Fig3_supp_SourceData.xlsx` | `11d18e79bcb8` | 206823 |
| `FigS2_SourceData.xlsx` | `27ec7bead7c9` | 57592 |
| `HMM_20261014_by_zone.tsv` | `d0cc563613cb` | 452 |
| `HMM_20261014_largeMSP_pair_comparison.tsv` | `9ee647a67a1c` | 338 |
| `HMM_20261014_merged_profile.tsv` | `7aace8ebf929` | 7012 |
| `HMM_20261014_occupancy_pair_points.tsv` | `3bdd95d1411a` | 193453 |
| `HMM_20261014_raw_vs_paired.tsv` | `21aa4c1f2f17` | 166 |
| `HMM_20261029_alu_by_zone.tsv` | `b0b4d814b369` | 557 |
| `HMM_20261029_label_rate_ratio.tsv` | `6506b547d0b1` | 136 |
| `HMM_20261029_paired_molecule_fraction.tsv` | `e8ae7c7ff3c1` | 82 |
| `ORIENTATION.tsv` | `cc5e5780f128` | 13460 |
| `PERARRAY_LANDSCAPE.tsv.gz` | `c5289d8dc99a` | 111585 |
| `PERARRAY_META.tsv` | `8298a482d3e8` | 9274 |
| `PICKED_CALLS.tsv` | `c1907a894396` | 299120 |
| `PICKED_MOLECULES.tsv` | `37500a163200` | 627 |
| `PROFILE_BP_v13.tsv.gz` | `ec247298584c` | 310483 |
| `PROFILE_META_v13.tsv` | `9e87223cad1d` | 25388 |
| `REF2_20261010c_merged_profile.tsv` | `faa1b22513b3` | 7355 |
| `REF2_20261010d_m6a_depth_denominator_by_position.tsv` | `ef38425678a7` | 4190 |
| `REF2_20261011_m6a_in_unit_by_state_cut508.tsv` | `02b100ddae69` | 173477 |
| `REF2_20261012_by_cell_positions_and_depth.tsv` | `0fbd2fb242c8` | 3524 |
| `SourceData_Fig1_FigS1.xlsx` | `ae322e3f3ec8` | 75293 |
| `SourceData_Fig4.xlsx` | `23b22a87f1b0` | 99909 |
| `SourceData_Fig5_FigS7.xlsx` | `fa60ca0f025f` | 478510 |
| `SourceData_Fig6_A-H.xlsx` | `dd310666d04d` | 11649 |
| `SourceData_Fig6_I-L.xlsx` | `3e586b26fc82` | 16657 |
| `SourceData_Fig7_FigS11_FigS12_FigS13.xlsx` | `efb88d00ded9` | 2960647 |
| `SourceData_FigS10.xlsx` | `54d0f6eb0492` | 241866 |
| `SourceData_FigS4.xlsx` | `3473269d0f58` | 222450 |
| `SourceData_FigS5.xlsx` | `4a90ca14e0ad` | 102396 |
| `SourceData_FigS6.xlsx` | `83e89a7fc0a0` | 46800 |
| `SourceData_FigS8.xlsx` | `4baa203b6405` | 10372 |
| `SourceData_FigS9_A-B.xlsx` | `25f4953bb476` | 18265 |
| `UNITBETA_PERARRAY_v13.tsv.gz` | `5d405e2a01b7` | 3341511 |
| `XB2_band_all.npz` | `b40a3757dfb9` | 1877 |
| `XB2_band_noodd.npz` | `e06b180d07a9` | 1888 |
| `cells/nucfree_GM12878all.json` | `5e342d4eb87d` | 5239 |
| `cells/nucfree_K562all.json` | `f3466e657600` | 4987 |
| `cells/sig_GM12878all.json` | `785e81f3cc2c` | 4349 |
| `cells/sig_GM12878all.npz` | `91dd00c0a7ba` | 26226 |
| `cells/sig_K562all.json` | `44e1a70e4c02` | 3600 |
| `cells/sig_K562all.npz` | `e0b42e9b4d59` | 11750 |
| `cells/strata_GM12878all_gene.json` | `8cd1c97e007f` | 3636 |
| `cells/strata_GM12878all_gene.npz` | `bf448a695782` | 21356 |
| `cells/strata_K562all_gene.json` | `cbf3567bf030` | 3300 |
| `cells/strata_K562all_gene.npz` | `66fe0dfedad6` | 16358 |
| `fhmm_peaks_5bp_bin.npz` | `770e363293f6` | 27191 |
| `fhmm_peaks_v2.npz` | `0408bc5a43a5` | 205594 |
| `gm_v5_contigsig_v2.json` | `bd718546936a` | 1868368 |
| `hprc_genebeta.npz` | `73612952a226` | 32446 |
| `k562_v5_contigsig_v2.json` | `fb2971da132f` | 572804 |
| `nucfree_bin.npz` | `d25e96217b56` | 4674 |
| `nucfree_flankcore_span.npz` | `33daa5dbe840` | 9097 |
| `occlusion_bin.npz` | `631b18a6aaf6` | 3449 |
| `r20_actuation_profile_REF.tsv` | `0c09bff49e3e` | 230917 |
| `r20_actuation_unified.tsv` | `74e195d6694d` | 869 |
| `ru_v11.parquet` | `1995b1901373` | 25732589 |
| `s01_element_bounds.tsv` | `b8ad392d397d` | 234 |
| `scan_maps_v2.npz` | `d8f6ffe89e40` | 62634 |
| `singlemol_fhmm_bin.npz` | `fb6af6b5e357` | 13374 |
| `site_confirm_callerfree.npz` | `c41a70751cb7` | 30584 |
