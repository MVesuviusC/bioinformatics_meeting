#!/bin/bash
#SBATCH --job-name=spaceranger_count
#SBATCH --output=spaceranger_count_%j.txt
#SBATCH --error=spaceranger_count_%j.txt
#SBATCH --time=3-00:00:00
#SBATCH --array=0-1
#SBATCH --cpus-per-task=20
#SBATCH --wait
#SBATCH --partition=general,himem

set -e

## set path info
export PATH=/gpfs0/home2/gdrobertslab/lab/Tools/10x/spaceranger-2.1.1:$PATH

sample_array=(ROB491D-1 ROB495D-2)
sample_name=${sample_array[$SLURM_ARRAY_TASK_ID]}

our_name_array=(SP0001 SP0002)
our_name=${our_name_array[$SLURM_ARRAY_TASK_ID]}

echo ${sample_name}
echo ${our_name}

base_path=/home/gdrobertslab/lab

cd ${base_path}/Counts_2/

spaceranger count \
    --id=${our_name} \
    --transcriptome=${base_path}/GenRef/10x-mm10_spaceranger \
    --probe-set=${base_path}/GenRef/10x-mm10_spaceranger/Visium_Mouse_Transcriptome_Probe_Set_v1.0_mm10-2020-A.csv \
    --fastqs=${base_path}/FASTQs_2/spatial_0001/240221_Miller-GSL-KM-3778/ \
    --sample=${sample_name} \
    --cytaimage=${base_path}/BCLs/spatial_0001/${sample_name}_framed.TIF \
    --image=${base_path}/BCLs/spatial_0001/${sample_name}.tif \
    --unknown-slide=visium-2-large \
    --localcores=16 \
    --localmem=140

