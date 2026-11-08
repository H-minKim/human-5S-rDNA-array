#!/bin/bash
# 위 글을 돌린 SLURM 배치 (12 작업 배열)
#   단계 : S1 B 상류
#   Method document : analysis/METHOD_CONSENSUS_MAP_G8_REBUILD.md
#SBATCH --job-name=s5s_jxhifi
#SBATCH --partition=cpu-short
#SBATCH --array=0-11%6
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=11:00:00
# NOTE (24 Oct 2026) the two paths below were absolute as run; they were made relative to take the
# user account out of the deposit. #SBATCH directives are read by sbatch before any shell variable
# exists, so RDNA_JX_WORK cannot be used here and there is no --chdir in this file: a relative path
# resolves against the SUBMIT directory. Submit from the junction-scan working directory and make
# sure logs/ exists there, or the array job will not start. Log destination only; no result changes.
#SBATCH --output=logs/jh_%A_%a.out
#SBATCH --error=logs/jh_%A_%a.err
# 검정 2 — HiFi junction 배치 (5S rDNA)
# %6 limits it to 6 at a time, because the cluster is shared
set -u
W="${RDNA_JX_WORK:?set RDNA_JX_WORK to the junction-scan working directory on the cluster}"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${RDNA_CONDA_ENV:?set RDNA_CONDA_ENV to the conda environment holding minimap2 2.28-r1209, samtools 1.9 and BLAST+ 2.16.0}"
command -v minimap2 >/dev/null || { echo "FATAL minimap2 없음"; exit 9; }
echo "host=$(hostname) task=$SLURM_ARRAY_TASK_ID minimap2=$(command -v minimap2)"
LIST=$W/src/part_${SLURM_ARRAY_TASK_ID}.txt
[ -s "$LIST" ] || { echo "목록 없음 $LIST"; exit 0; }
TAG=a${SLURM_ARRAY_TASK_ID} TH=8 TMPD=$W/tmp/$SLURM_ARRAY_TASK_ID \
  python3 $W/src/jxhifi_snu.py "$LIST"
