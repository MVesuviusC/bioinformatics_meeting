/*
 * Define "align_hisat" process
 *
 */
process align_hisat {
    publishDir params.align,  mode: "copy", pattern: "*.{bam,bam.bai,txt}"
    publishDir params.slurmOut, mode: "copy", pattern: "*.log"

    input:
    path (r1_file)

    output:
    path ("*.bam"),         emit: bam
    path ("*.bam.bai"),     emit: bai
    path ("*.txt"),         emit: txt
    path ("*.log"),         emit: log


    shell:
    '''
    set -e ### stops bash script if line ends with error
    
    ml purge
    ml load HISAT2/2.2.1 \
            SAMtools/1.15
    
    # Get the basename of the R1 file by trimming off the path and .fastq.gz
    r1_base=!{r1_file}
    r1_base=${r1_base##*/}
    r1_base=${r1_base%.fastq.gz}
    
    echo ${HOSTNAME} !{r1_file}
    
    # Make variable with ref here since it's long
    ref=/reference/mus_musculus/GRCm38/ensembl/release-86/Sequence/Hisat2Index/genome
    
    # Run hisat2
    hisat2 \
            -x ${ref} \
            -U !{r1_file} \
            -k 1 \
            -p 10 \
            --summary-file ./${r1_base}.txt \
        | samtools sort \
            -@ 5 \
            -m 5G \
            -O BAM \
        | samtools view \
            -F 3 \
            -b \
            - \
        > ./${r1_base}.bam
    
    samtools index -@ 5 ./${r1_base}.bam

    # save the slurm output
    cat .command.log > ${SLURM_JOB_ID}_align_${r1_base}.log
    '''
}

