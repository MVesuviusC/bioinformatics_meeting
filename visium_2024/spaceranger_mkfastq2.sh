#!/bin/bash
#SBATCH --job-name=spaceranger
#SBATCH --output=spaceranger_%j.txt
#SBATCH --error=spaceranger_%j.txt
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=30
#SBATCH --wait
#SBATCH --partition=general,himem

set -e

## set path info
export PATH=/gpfs0/home2/gdrobertslab/lab/Tools/10x/spaceranger-2.1.1:$PATH
export PATH=/gpfs0/home2/gdrobertslab/lab/Tools/bcl2fastq/2.20/bin:$PATH

## demultiplex bcl and create fasta files
spaceranger mkfastq \
  --localcores=30 \
  --localmem=200 \
  --samplesheet=240227_VH01085_30_AAFGFVJM5_Roberts/SampleSheet.csv \
  --run=240227_VH01085_30_AAFGFVJM5_Roberts \
  --output-dir /home/gdrobertslab/lab/FASTQs_2/spatial_0001
