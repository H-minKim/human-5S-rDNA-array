# tables/ — the result tables the figure code reads

121 tables, 43.7 MB. Each one is the exact file a figure was drawn from: it was
selected by matching the `input_table_md5` recorded in the `MANIFEST_*.tsv` files,
not by file name, because several of these names exist in more than one version
folder in the working tree with different content.

Point the table roots at this folder when re-running the code, for example
`TABLE_DIR`, `DATA_ROOT`, `P11_TABLE_DIR`, `P10_DIR` and `REGISTRY_DIR`
(the full list is in the repository README, section *Running this code*).

**4 tables are not here because of their size**; they are in the Zenodo record
that archives this repository:

| MB | table |
|--:|---|
| 68.7 | `P10_FRAME_ALL.tsv` |
| 65.2 | `P10_CPG_ADD3.tsv.gz` |
| 57.0 | `L1_READ_UNIT.parquet` |
| 54.0 | `P10_CPG_ALL.tsv.gz` |

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
| `BOUNDARY_WIDE_H2_508.tsv` | `963331f95ce8` | 882512 |
| `CONSENSUS_MAP_G8_full.tsv` | `bb68fe3a2a23` | 8864818 |
| `CONSENSUS_SUM_G8_full.tsv` | `4c8252462824` | 21186 |
| `COSHARE_PERHAP.tsv` | `398d1259d27a` | 844653 |
| `COV_AFTER_FILTER.tsv` | `55cec51e3838` | 205808 |
| `DC_TABLE.tsv` | `6029cf9c721e` | 21483 |
| `DISP_PROFILES.tsv` | `0599c32be6b0` | 33231 |
| `DISP_SITES_TABLE.tsv` | `4edb4f1dad70` | 11740 |
| `EXAMPLE_PROFILE.tsv` | `ddf56e9d56c0` | 3675 |
| `F8X2_pairs.tsv` | `2ce60ee20897` | 621 |
| `FIG1_2_MID.tsv` | `a7638c89acf3` | 678 |
| `FIG1_2_UNITS.tsv` | `fc4eb1dba1c8` | 1950 |
| `FIG1_V5_A_bp_HC.tsv` | `c30a31686bc7` | 29542 |
| `FIG2_data_hifi.tsv` | `1ff07c00c41a` | 624 |
| `FIG2_data_unitcmp.tsv` | `978bda6f0176` | 6666 |
| `FIG3_A_UNITPROFILE_ZONE.tsv` | `a315cdac6f4a` | 31581 |
| `FIG3_B_POSBETA.tsv` | `e6b0b57f9d6b` | 11495 |
| `FIG3_C_LOWCOUNT.tsv` | `580823d3003c` | 751 |
| `FIG3_G_TYPEDIST.tsv` | `f893f27bea9a` | 17433 |
| `FIG3_SUPP_A_SEQ1K.tsv` | `590470d4af17` | 11005 |
| `FIG5V2_PAIRZONE5.tsv` | `f1cfef234fb3` | 3554 |
| `FIG8_PAIRS.tsv` | `d83d18590dc7` | 1895 |
| `FIGX25_DENSITY.tsv` | `92b829838fe4` | 27495 |
| `FIGX25_PROFILE.tsv` | `cbb8972f3b91` | 4628 |
| `FIG_data_table.tsv` | `2c0c009285c9` | 1233 |
| `FILTER_DECOMPOSE.tsv` | `fdfa9107be64` | 241 |
| `FILTER_DECOMPOSE_BDP1_2026-11-01.tsv` | `4d24ab436b8a` | 116 |
| `FILTER_DECOMPOSE_GMATAC_2026-11-01.tsv` | `13f08b85afcd` | 76 |
| `GM12878_interaction_NEWZONE_GMWGS_2026-10-14.tsv` | `b2355947b209` | 139 |
| `GM_QUALGATE_2026-11-01.tsv` | `0a596a6b423b` | 32158 |
| `INTEGRATED_14SITES.tsv` | `0d518807d74e` | 1591 |
| `K562_QUALGATE_MIX.tsv` | `1414cf754e97` | 37510 |
| `K562_QUALGATE_MIX_BDP1BRF1.tsv` | `c627d845531e` | 9931 |
| `K562_QUALGATE_NEWZONE_2026-10-13.tsv` | `c8a0b6b322da` | 22756 |
| `K562_interaction_NEWZONE_withBDP1BRF1_2026-10-14.tsv` | `270dc585ce15` | 411 |
| `K562_zoneratio_Q20_NEWZONE_withBDP1BRF1_2026-10-14.tsv` | `31c86607c38e` | 740 |
| `L2_HAP_UNIT.parquet` | `2d0217ca81b8` | 9477685 |
| `MAP_G8_GM12878FS.tsv` | `c5fa06b45060` | 9205700 |
| `MAP_G8_K562r.tsv` | `a91b74d6df1b` | 9201404 |
| `MAP_hapB_v4_110units.tsv` | `68f5caa084c2` | 18897 |
| `METH_TERM_units.tsv` | `fad6e6a66b4a` | 232924 |
| `MODEL_TEST_ALL6.tsv` | `ae81acc0f287` | 5368 |
| `NORM_VS_CANCER_SAMPLE.tsv` | `977c4a3e441a` | 2299 |
| `P08_X2_zone_region.tsv` | `584cbac13842` | 3287 |
| `P10_ARRAYS.tsv` | `4299508a8ca8` | 32060 |
| `P10_CPG_SEGS_ADD3.tsv` | `a95f2e09f8dd` | 9186 |
| `P10_DXZ4_POLARITY_FINAL.tsv` | `a6a120f498ee` | 1689 |
| `P10_FRAME_SEGS_ALL.tsv` | `4108ba22cc74` | 15631 |
| `P10_ONOFF_TESTS.tsv` | `b1add5f9ad93` | 854 |
| `P10_UNIT_STEPS_ABS.tsv` | `f95b9e89a86c` | 248 |
| `P11_GM12878_MARKERS_NEWZONE_2026-10-13.tsv` | `4fec23b0c0a3` | 3367 |
| `P11_K562_MARKERS_NEWZONE_2026-10-13.tsv` | `7a7ddea532af` | 2875 |
| `P11_percellline_NEWZONE_GMWGS_bootstrap_2026-10-14.tsv` | `c4c303b1b5d6` | 2749 |
| `POP_COHORT.tsv` | `82edeaab20ea` | 52499 |
| `PROB_A4_DIST_SITE_TONE.tsv` | `aab77f6cd3d7` | 579566 |
| `PROB_A4_SITE.tsv` | `8e0915d1dda2` | 4883 |
| `PROB_A4_SITE_TONE.tsv` | `8a06d5f3f884` | 21807 |
| `Q1_results.tsv` | `b1a29eaec02f` | 1697 |
| `Q2b_results.tsv` | `8638d02a23db` | 7323 |
| `REGCHAR_axis.tsv` | `d2da5973f84a` | 46409 |
| `REGCHAR_relcount.tsv` | `97d494a4469b` | 2955 |
| `RU_GM12878FS.tsv.gz` | `1fffbcea765f` | 1709527 |
| `RU_K562r.tsv.gz` | `f8b13b7421bb` | 380352 |
| `RU_K562u.tsv.gz` | `5e0799454a25` | 397158 |
| `S14_length.tsv` | `1de63221dd12` | 484 |
| `S18_changepoint.tsv` | `fb7c0859f602` | 108993 |
| `S1_per_position_profile.tsv` | `98996f44644b` | 169950 |
| `S24_gain_curve.tsv` | `5adce291cd1f` | 20837 |
| `S24_optimal_zoning.tsv` | `bbe1c8dab70d` | 1347 |
| `S26b_boundary_drivers.tsv` | `dfc5b7df3b36` | 1451 |
| `S26b_position_profile.tsv` | `097ff70985d7` | 10299 |
| `STAGE_ALL_V11F2.tsv` | `683992af188f` | 59189 |
| `SUBUNIT_BIMODAL2.tsv` | `a6ab3256ba53` | 14980 |
| `SUPP3_2_POSSD.tsv` | `5b1c646147a2` | 2108 |
| `SUPP3_B_HAP3.tsv` | `4650d5dd5fa0` | 267349 |
| `SourceData_SuppFig1.xlsx` | `63c58b4b03e9` | 18338 |
| `T1_pair_summary.tsv` | `960664742fc1` | 14597 |
| `TIERS_H0H3_SLOTS_2026-10-05.tsv` | `fc34c6326111` | 94793 |
| `TRANSCRIBED_FOR_PANEL_A.tsv` | `0eab96ae032d` | 1806 |
| `TableS1_Hia5_untreated_control_by_position.tsv` | `07b25442b25a` | 949 |
| `TableS3_alternative_explanations_predicted_vs_observed.tsv` | `7350a4eea9a7` | 1888 |
| `TableS7_nucleosome_removed_occupancy.tsv` | `4644d7f1856d` | 781 |
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
| `X1_population_control.tsv` | `a19205eabc85` | 1565 |
| `X4_lag2_and_patch.tsv` | `18cdbbe1c221` | 1223 |
| `X5_lag2_noncircular.tsv` | `4ce156fb8f62` | 1056 |
| `X7_homogenisation_evidence.tsv` | `0db41b93c163` | 910 |
| `X9_patch_floor.tsv` | `2f5368ed3d30` | 610 |
| `XA_decay_constant_robustness.tsv` | `462ec3f54a59` | 606 |
| `Y4_cellline_summary.tsv` | `3e59d7e8fc71` | 405 |
| `Y5_terminal_defining_variants.tsv` | `b013af3755fa` | 3210 |
| `ZB_neighbour_divergence.tsv` | `be24394880b3` | 398 |
| `ZM_population_carrier.tsv` | `f9d80cbf0aad` | 670 |
| `ZM_terminal_similarity.tsv` | `5398b28a3f9c` | 551 |
| `ZN_depthset_census.tsv` | `b3343a3840c1` | 2094 |
| `ZP_gap_shape.tsv` | `17ddf3048255` | 320 |
| `cohort_tiers_Fig1_2026-10-08.tsv` | `f3e56d90ad6c` | 728 |
| `fig_VTORG_single_region_type_composition_source_values.tsv` | `64ffa4f303e1` | 27296 |
| `fig_VTPD_region_type_composition_source_values.tsv` | `914783b86972` | 32514 |
| `hap_seq_node_v8.tsv` | `bc11cdc3aec7` | 25342 |
| `hor_anchor_per_read_unit.tsv` | `7ac0f3445d31` | 325340 |
| `lowmeth_changepoint_4haplotypes_2026-10-12.tsv` | `af278c7324de` | 607 |
| `peaks_both_v2.tsv` | `61dac6219562` | 4427 |
| `privvar_per_read_unit.tsv` | `e6daf9ef2296` | 6098 |
| `unit_marker_meth.tsv` | `52fcf2405a86` | 174271 |
