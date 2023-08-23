#!/bin/sh
#SBATCH --array=0-23
#SBATCH --error=slurmOut/align-%j.txt
#SBATCH --output=slurmOut/align-%j.txt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --job-name align
#SBATCH --wait
#SBATCH --time=0-12:00:00

set -e ### stops bash script if line ends with error

ml purge
ml load HISAT2/2.2.1 \
        SAMtools/1.15

# make array with all the input R1 files
r1_array=(input/fastqs/*.fastq.gz)
# Use SLURM_ARRAY_TASK_ID to select the correct R1 file for this index
r1_file=${r1_array[${SLURM_ARRAY_TASK_ID}]}

# Get the basename of the R1 file by trimming off the path and .fastq.gz
r1_base=${r1_file##*/}
r1_base=${r1_base%.fastq.gz}

echo ${HOSTNAME} ${SLURM_ARRAY_TASK_ID} ${r1_file}

# Make variable with ref here since it's long
ref=/reference/mus_musculus/GRCm38/ensembl/release-86/Sequence/Hisat2Index/genome

# Run hisat2
hisat2 \
        -x ${ref} \
        -U ${r1_file} \
        -k 1 \
        -p 10 \
        --summary-file output/align/${r1_base}Summary.txt \
    | samtools sort \
        -@ 5 \
        -m 5G \
        -O BAM \
    | samtools markdup \
        -@ 5 \
        -r - - \
    | samtools view \
        -F 3 \
        -b \
        - \
    > output/align/${r1_base}.bam

samtools index -@ 5 output/align/${r1_base}.bam
