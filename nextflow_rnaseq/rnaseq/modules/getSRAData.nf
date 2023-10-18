/*
 * Define "getSRAData" process
 *
 */
process getSRAData {
    publishDir params.fastqs,  mode: "copy", pattern: "*.{gz}"
    publishDir params.slurmOut, mode: "copy", pattern: "*.log"

    input:
    val (sra_number)

    output:
    path("*.gz"),        emit: fastq
    path("*.log"),          emit: log


    shell:
    '''
    module purge
    module load SRAToolkit/3.0.1

    # download the SRA file (this is not the fastq file)
    # -O means to put the SRA file in the specified directory
    prefetch \
        -O ./ \
        !{sra_number}
    
    # convert the SRA file to fastq
    # -e 4 means to use 4 threads
    # -S means to split the fastq into two files (R1 and R2)
    #    In this case, all of the reads are single end, so there is only one file
    # -O means to put the fastq files in the specified directory
    fasterq-dump \
        -e 4 \
        -S \
        -O ./ \
        ./!{sra_number}/!{sra_number}.sra
    
    pigz ./!{sra_number}.fastq

    # save the slurm output
    cat .command.log > ${SLURM_JOB_ID}_fastqs_!{sra_number}.log
    '''
}

