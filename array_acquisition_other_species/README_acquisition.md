# How the 5S rDNA arrays were obtained in the non-human species

This folder holds the code that produced the array catalogue, the canonical array
references and the read placements for the eight non-human species of Figure 7 and
Supplementary Figures S11 to S13.

Unlike the scripts under `upstream/`, these do not carry the absolute paths of the
machines that ran them; section 4 says what was replaced.

**No assembly was performed.** Every genome used for these species is a published assembly,
downloaded by accession through the NCBI Datasets API. A sweep of all 114 scripts here
and the 113 recovered from the compute nodes for
`hifiasm|verkko|flye|canu|shasta|wtdbg|raven|miniasm|spades` returns no hit. The one file
file built here is the zebra finch diploid reference, which is not an assembly but
the concatenation of two published assemblies (`cat` of the maternal and paternal FASTA).
This is a point of difference
from the human part of the study, where haplotypes assembled with `hifiasm` are used
alongside HPRC assemblies. The two must not be described together.

## 1. The six steps

| # | Step | Code | The rule that was applied | Output |
|---|---|---|---|---|
| 1 | Download the published assembly | `ape_prep.sh`, `zfinch_v03.sh`, `zfB_setup.sh`, `ath_*` | NCBI Datasets API by accession, `samtools faidx`; no assembly step | `asm/<species>.fna` |
| 2 | Locate the 5S genes | `cms_all.sh`, `ape_arrays.sh`, `zfinch_v03.sh` | `cmsearch` (Infernal 1.1.5) against Rfam `RF00001.cm` (CM length 120, gathering threshold 38.00 bits, checksum 242171328) with `--noali --cut_ga`, that is at the default gathering threshold with no parameter changed | `cms_tbl/<species>.tbl` |
| 3 | Group the genes into arrays | `array_catalog.py` | **Hits no more than 10 kb apart belong to one array**; an array is catalogued when it holds **three or more copies**; **copy number is the number of `cmsearch` hits**, never a length divided by a period; assembly gaps within 100 kb are checked | `5S_array_catalog.tsv` |
| 4 | Measure structure and period | `array_struct.py`, `ladder_all.py`, `period_direct2.py`, `selfdot.py` | **The repeat period is measured directly, from the offset ladder of a self dot plot of the array sequence; the repeat length given in a reference annotation is not used**; gene and spacer GC and CpG are measured separately; 5S genes and gaps within 20 kb of flank are checked | `5S_array_structure.tsv`, `LADDER_ALL.tsv` |
| 5 | Build the canonical array reference | `ape_arrays.sh`, `refunit_build.py` | The reference spans the **array plus 20 kb of flank on each side**, and is **stored in the direction of transcription of the gene**; one repeat unit runs from one base past the end of a gene to the end of the next gene, so a unit is `[spacer][gene]`; a unit boundary is accepted only on an **intact gene, covered at 0.80 or more** | `01_canonical/<species>/*.fa`, `02_reference_units/` |
| 6 | Place the reads | per species, see section 2; `read_units.py`, `units_matrix_g8.py` | Reads are aligned to the **whole assembly** and the array interval is then extracted, never aligned to an array-only reference; `samtools view -F 0x900` removes supplementary and secondary records; **reads are not filtered on mapping quality**, because inside a tandem array MAPQ 0 is the normal case and a `-q 20` filter empties the window | `04_bam/<species>/*.5S.bam`, `03_reads/` |

The 10 kb and three-copy rules of step 3 were applied identically to every species; no
threshold was set per species.

## 2. Read preparation differs between species

| Species | Read source | What was done | Re-basecalled |
|---|---|---|---|
| Great apes and siamang | GenomeArk, per-basecaller BAM with `MM`/`ML` preserved | candidate reads selected by 5S k-mer while streaming the BAM from object storage, then realigned to that species' own assembly | no (as released: `dorado-0.9.1`, `guppy-6.3.7`, `guppy-6.3.8`) |
| Zebra finch | GenomeArk raw FASTQ, 13 flowcells | candidate reads selected by four 22-mers of the 5S gene and their reverse complements while streaming, aligned to the diploid reference, then re-basecalled | yes (`dorado 2.1.2`, 9 of the 13 flowcells) |
| Zebrafish | ENA, ONT BAM aligned to GRCz11 | reads over the 5S intervals extracted and realigned to GRCz12tu | no (as deposited, `Guppy 6.4.2`) |
| *A. thaliana* | ENA raw signal, 164 GB of pod5 | 5S reads selected while streaming, then re-basecalled on CPU in all three contexts | yes (`dorado 2.1.2`, CG, CHG and CHH) |

Because the basecaller version shifts methylation, cross-species comparisons are made only
between read sets called with the same version; the per-flowcell versions are recorded in
the data catalogue of section 3.

## 3. Per-species values

| Species | Assembly read | `cmsearch` hits | Arrays catalogued | Arrays with structure | Arrays with reads placed | Reference unit (bp) |
|---|---|--:|--:|--:|--:|---|
| *Pan troglodytes* | `GCA_028858775.3` | 485 | 1 | 1 | 1 | 2,218 |
| *Pan paniscus* | `GCA_029289425.3` | 588 | 1 | 1 | 1 | 2,242 |
| *Gorilla gorilla* | `GCA_029281585.3` | 539 | 1 | 1 | 1 | 2,248 |
| *Pongo pygmaeus* | `GCA_028885625.3` | 514 | 1 | 1 | 1 | 2,313 |
| *Pongo abelii* | `GCA_028885655.3` | 537 | 1 | 1 | 1 | 2,297 |
| *Symphalangus syndactylus* | `GCA_028878055.3` | 421 | 1 | 1 | 1 | 3,296 |
| *Taeniopygia guttata* | `GCA_048771995.1` + `GCA_048772025.1` | 2,238 | 33 | 8 | 8 | 368-372 |
| *Danio rerio* | `GCA_049306965.1` | 8,638 | 46 | 46 | 31 | 151-9,990 |
| *Arabidopsis thaliana* | `GCA_028009825.2` | 3,935 | 10 | 10 | 9 | 475-749 |
| *Homo sapiens* (pipeline check only) | `GCA_009914755.4` | 489 | 1 | 1 | - | 2,240 |

**Three array counts, three definitions.** The columns above narrow: everything `cmsearch`
proposed, then what a structure was measured on, then what reads were placed on. A single
number quoted without its definition will disagree with another table. For the zebra finch,
33 catalogued rows become 8 because only four loci on two haplotypes were kept, and 25 of
those rows come from two other assemblies that were screened and not used
(`GCA_051427915.1`, 19 rows; `GCA_051428105.1`, 6 rows). For zebrafish, 46 become 31 because
only those arrays allowed a unit frame to be built. **The three-to-five-copy arrays that
Figure 7A draws in grey are not excluded by copy number**: copy number is what they have in
common, not the rule that excluded them.

The siamang is in this folder but not in the analysis. Its assembly-to-read agreement is
0.094 on the strict measure against 0.696 to 0.982 in the other species, so it was not used
for unit-level work; its single array of 3,837 units is excluded from every figure.

The human row is a check that the same pipeline run on human recovers the known array; no
reads were placed. The human results of this study do not use a single assembly as their
frame: they use a flank-to-flank array reference built per haplotype.

## 4. Paths are given by root variables

The scripts were written to run on one workstation and on three compute nodes. The paths of
those machines appear here as variables, so no machine-specific location is carried in the
code.

| Variable | What it points at |
|---|---|
| `$RDNA_OTHERSP_DIR` | the folder holding this campaign's data and outputs |
| `$SCRATCH_DIR` | a scratch folder for intermediate files |
| `$NODE_KEY` | the private key used to reach the compute nodes |
| `$NODE_HOST_1`, `$NODE_HOST_2`, `$NODE_HOST_3` | the three compute nodes |

No value is recorded here for the key path or the three host addresses: they are
credentials and infrastructure, not inputs a reader would follow. Nothing in this folder
needs them in order to be read; they are needed only to re-run the work, which requires
new nodes in any case.

`node_scripts/node1/zfwgbs/zf_prep.sh` names `~/.ssh/id_rsa`. That is the default key name
under a node's own home directory, the same string on any machine; it identifies nothing and
no file in this deposit.

## 5. Which scripts are current

`CODE_INDEX.tsv` lists all 110 scripts in `code/` with an `mtime`, a `status` and a one-line
description. 87 are marked `current` plainly, and the rest carry a status with a reason in
parentheses. Superseded and withdrawn files are kept rather than deleted so that
a number traced to an older table can still be matched to the code that made it. **Read the
status before running anything in `code/`.**

The five statuses are `current`, `superseded` (the parenthesis names the file that replaced
it), `withdrawn` (with why it was not used), `decision record`, and `reference only`. The two
rows whose description is empty, the two `.bed` files, have no description.

The file names in section 1 are the entry points.

## 6. The compute-node scripts are here too

The upstream extraction, alignment and methylation calling ran on three Linux compute nodes.
The scripts that ran there are in `node_scripts/`, one folder per node: `node3` 44 files,
`node2` 41, `node1` 28, **113 files in all**.

They are deposited as plain files rather than as archives, so that a text search of this
folder reaches them.

Where a node copy and a local copy of the same script differ, **the node copy is the one that
ran**.
