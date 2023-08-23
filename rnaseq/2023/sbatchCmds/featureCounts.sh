#!/bin/sh
#SBATCH --error=slurmOut/featureCount-%j.txt
#SBATCH --output=slurmOut/featureCount-%j.txt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --job-name featureCount
#SBATCH --wait
#SBATCH --time=0-12:00:00

set -e ### stops bash script if line ends with error

ml purge
ml load GCC/8.3.0 \
        Subread/2.0.0

# Make variable with gtf file here since it's long
gene_gtf=/reference/mus_musculus/GRCm38/ensembl/release-86/Annotation/Genes/gtf/Mus_musculus.GRCm38.86.gtf

featureCounts \
  -T 10 \
  -a ${gene_gtf} \
  -o output/geneCounts/geneCountTable.txt \
  -s 2 \
  output/align/*.bam
