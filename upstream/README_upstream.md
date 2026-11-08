# `upstream/` — code from the shared v11 pipeline, copied unchanged

These scripts belong to the shared v11 reference pipeline that several figures
depend on. They are deposited because the numbers of published panels come out of
them.

They are byte-for-byte copies of the originals that ran. Nothing was edited, not
even a path, so that the code here is demonstrably the code that produced the
tables.

| File here | md5 | Lines | Original folder | What it does |
|---|---|---|---|---|
| `v11_node/meth_units_jx_v11.py` | `369925607294` | 195 | the pipeline's `node` folder | reads the `MM`/`ML` tags of the batch BAMs and writes one read-by-unit methylation table per node |
| `v11_node/meth_units_ao_v11.py` | `3f55989d7bac` | 195 | the pipeline's `node` folder | the same for the second placement set |
| `v11_local/merge_meth.py` | `0ad651bd8314` | 98 | the pipeline's `local` folder | merges the per-node tables into `READMAP_METH_G8V4.tsv` and `METH_SITE_G8V4.tsv.gz` — **those two names and no others** |

**Why the folder names carry `node` and `local`.** The originals live in two
folders, and the same filename exists in both with different content
(`meth_units.py`, `merge_meth.py`). Naming the copies after the figure they serve
would lose that, so each keeps the side it came from.

**The node scripts need a compute node to run again**, because the reads they
read are there. The output filename ends in the node's own hostname, which the
script takes from `hostname` at run time — so a file called
`METH_READUNIT_<something>_ip-172-31-…​.tsv.gz` was named by the machine, not by a
person, and the hostname is not a network address anyone needs.

**What these scripts do not do.** They do not call a basecaller. They read the
`MM`/`ML` base-modification tags that the public submissions already carry, count
a site as methylated at `ML >= 128`, and report sites with no tag separately.
This is the direct evidence for the Methods statement that 5mC calls are taken
from the submitted tags rather than produced here.

**Do not attribute a table to a script here without matching the written filename.**
`merge_meth.py` writes `READMAP_METH_G8V4.tsv` and no other name. A similar name in
a similar pipeline is not a provenance.
