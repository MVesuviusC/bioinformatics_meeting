#!/bin/sh
#SBATCH --account=gdrobertslab
#SBATCH --array=0-23
#SBATCH --error=slurmOut/getSRA.txt
#SBATCH --output=slurmOut/getSRA.txt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --job-name getSRA
#SBATCH --wait
#SBATCH --time=24:00:00

module load module load SRAToolkit/3.0.1

# Pull the SRA accession number from the csv file and make a bash array
# The outer parenthesis make it an array
# The inner parenthesis with a $ in front of it make it a command that is run
sra_array=($(cut -f 1 -d "," misc/MiiiceIiiiinSpaaaaaaaaaaaaace.csv | grep -v "Run"))

# Pick the SRA number to download with this array index
sra_number=${sra_array[${SLURM_ARRAY_TASK_ID}]}

# download the SRA file (this is not the fastq file)
# -O means to put the SRA file in the specified directory
prefetch \
    -O input/sra \
    ${sra_number}

# convert the SRA file to fastq
# -e 4 means to use 4 threads
# -S means to split the fastq into two files (R1 and R2)
#    In this case, all of the reads are single end, so there is only one file
# -O means to put the fastq files in the specified directory
fasterq-dump \
    -e 4 \
    -S \
    -O input/fastqs \
    input/sra/${sra_number}/${sra_number}.sra

pigz input/fastqs/${sra_number}.fastq

# Get rid of the SRA file folder that was created by prefetch
rm input/sra/${sra_number}/${sra_number}.sra
rmdir input/sra/${sra_number}
