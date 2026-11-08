# Upstream — v12 Fiber-seq node pipeline

Archived on 2026-10-22 at the 's request.
**26 files** (23 at first archiving; three producers added the same day — see the correction in §2).
**Nothing here is edited in this repository.** These files are copies of code that ran on the three
AWS compute nodes; they are kept so that the chain behind the published numbers is visible, not so
that it can be re-run from this folder.

---

## 1. Which figures hang on this code

| Figure | Through |
|---|---|
| **Fig. 5 A (lower profile), D, E** | `m6a_site/`, `m6a_readunit/` → `zg2_m6a_frac_profile.py`, `zk1b_m6a_by_state_cut508.py` |
| **Fig. 5 B, C, F, G, H** · **Suppl. Fig. S7 A, B, C** | the placed unit tables → `a13_beta_mixture.py` → `zy7_make_hmm_tables.py`, `zg1_*`, `zt2_*` |

Other figures in this paper rest on the **v11** stage, which is archived separately in
`upstream/v11_node/` and `upstream/v11_local/`. This folder is the **v12** stage only.

---

## 2. The run scripts

Four run scripts drove the pipeline, and they are deposited with the Python they
call, because they are the record of which settings ran:

```
04_placement_s5/<node>/run_s5.sh        the merge at MINSITE=5 and the re-call passes
04_placement_s5/<node>/run_meth_s5.sh   the m6A/5mC extraction
01_code/<node>/run_aojv.sh              read placement
01_code/<node>/setup_v12fs.sh           the working tree the placement runs in
```

Two near-identical names are worth separating, because the wrong one would give
a different result:

| Not used here | Why |
|---|---|
| `merge_aojv.py` | `MINSITE = 10` is hardcoded. The data was produced by **`merge_aojv_v2.py`**, which takes `MINSITE` from the environment and was called as `MINSITE=5` in `run_s5.sh`. |
| `run_fsmeth.sh` | writes `fsmeth_v12b`. The data files are `fsmeth_v12c_*` and `pbmeth_v12d_*`, which `run_meth_s5.sh` produced. |
| `merge_aojv_v3.py` | a later `_s6` series used by other analyses. Figure 5 and Supplementary Figure S7 do not read it. |

There is no job scheduler in this pipeline: `#SBATCH`, `srun`, `qsub`, `bsub`,
`module load` and `conda activate` do not appear anywhere in it. The only
parallelism is `nohup python3 ... &` in `run_aojv.sh`.
## 3. Node copies

The pipeline ran on three nodes (`iA`, `iB`, `w3`, internal hostnames `ip-172-31-32-12`,
`ip-172-31-57-143`, `ip-172-31-58-58`). **All fourteen `.py` files are byte-identical across the
three nodes**, so one copy is archived here; `MD5_node_s5.txt` is the node-side checksum list.
`MANIFEST_upstream_v12.tsv` records, for every file, the md5 and the folder it was taken from.

The three `recall_*_v10.py` files are the **sources**, not themselves run: `run_s5.sh` derives
`recall_*_fs3.py` and `recall_*_pb2.py` from them with a single `sed` substitution of the SAM
directory. Both the sources and the derived scripts are archived, so the substitution can be
checked. The sources live in the v11 archive (`08_versions/v11/01_code/node/`).

---

## 4. Randomness

This stage contains the only random draws anywhere in the Fig. 5 chain, and every one is seeded.

| File | Draw | Seed |
|---|---|---|
| `fsplace_v12.py`, `pbfsplace_v11.py` | `random.Random(rid).choice(...)` when a read ties between array slots | `random.seed(20260909)`; the per-read generator is seeded with the **read id**, so the choice does not depend on processing order or on how the work is split across nodes |
| `recall_edge_fs3.py`, `recall_edge_pb2.py` (and their v10 sources) | `random.sample(oks, 20)` | `random.seed(20260911)`; used for a spot check written to the run log, not for any tabulated value |
| `fsmeth_v12c.py` | — | no random draw at all |

Counting randomness only in the local analysis scripts would report "none", and that would be
wrong. This table is the reason the Methods can say the pipeline is reproducible.

---
