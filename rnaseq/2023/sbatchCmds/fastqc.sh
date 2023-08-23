#!/bin/sh
#SBATCH --array=0-23
#SBATCH --error=slurmOut/fastqc-%j.txt
#SBATCH --output=slurmOut/fastqc-%j.txt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --job-name fastqc
#SBATCH --wait
#SBATCH --time=0-12:00:00

set -e ### stops bash script if line ends with error

ml purge
ml load FastQC/0.11.9-Java-11.0.2

# make array with all the input R1 files
r1_array=(input/fastqs/*.fastq.gz)

# Use SLURM_ARRAY_TASK_ID to select the correct R1 file for this index
r1_file=${r1_array[${SLURM_ARRAY_TASK_ID}]}

# Get the basename of the R1 file by trimming off the path and .fastq.gz
r1_base=${r1_file##*/}
r1_base=${r1_base%.fastq.gz}

echo ${HOSTNAME} ${SLURM_ARRAY_TASK_ID} ${r1_file}

# Run fastqc
fastqc \
    -o output/fastqc \
    -t 4 \
    --extract \
    ${r1_file}

# We don't need this since it was already extracted
rm output/fastqc/${r1_base}_fastqc.zip