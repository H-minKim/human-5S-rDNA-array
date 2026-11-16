# `upstream/` — code from the shared v11 pipeline, copied unchanged

These scripts belong to the shared v11 reference pipeline that several figures
depend on. The numbers of published panels come out of them.

The paths in them are the ones they ran with.

| File here | md5 | Lines | Original folder | What it does |
|---|---|---|---|---|
| `v11_node/meth_units_jx_v11.py` | `80f36bf5ec5b` | 196 | the pipeline's `node` folder | reads the `MM`/`ML` tags of the batch BAMs and writes one read-by-unit methylation table per node |
| `v11_node/meth_units_ao_v11.py` | `781bf20a345a` | 196 | the pipeline's `node` folder | the same for the second placement set |
| `v11_local/merge_meth.py` | `55064038fe1d` | 98 | the pipeline's `local` folder | merges the per-node tables into `READMAP_METH_G8V4.tsv` and `METH_SITE_G8V4.tsv.gz` — **those two names and no others** |

**The folder names carry `node` and `local`.** The originals live in two
folders, and the same filename exists in both with different content
(`meth_units.py`, `merge_meth.py`), so each copy keeps the side it came from.

**The node scripts need a compute node to run again**, because the reads they
read are there. The output filename ends in the node's own hostname, which the
script takes from `hostname` at run time, as in
`METH_READUNIT_<something>_ip-172-31-…​.tsv.gz`.

**What these scripts do not do.** They do not call a basecaller. They read the
`MM`/`ML` base-modification tags that the public submissions already carry, count
a site as methylated at `ML >= 128`, and report sites with no tag separately.
The 5mC calls are taken from the submitted tags rather than produced here.

**Do not attribute a table to a script here without matching the written filename.**
`merge_meth.py` writes `READMAP_METH_G8V4.tsv` and no other name. A similar name in
a similar pipeline is not a provenance.
