/*
 * Define "fastqc" process
 *
 */
process fastqc {

    publishDir params.fastqc,  mode: "copy", pattern: "SRR*"
    publishDir params.fastqc,  mode: "copy", pattern: "*.html"
    publishDir params.slurmOut, mode: "copy", pattern: "*.log"

    input:
    path (r1_file)

    output:
    path ("SRR*"),          emit: fastqc_out
    path ("*.html"),        emit: html
    path ("*.log"),         emit: log

    shell:
    '''
    set -e  ### stops bash script if line ends with error 

    ml purge
    ml load FastQC/0.11.9-Java-11.0.2
    
    # Get the basename of the R1 file by trimming off the path and .fastq.gz
    r1_base=!{r1_file}
    r1_base=${r1_base##*/}
    r1_base=${r1_base%.fastq.gz}

    echo ${HOSTNAME} !{r1_file}
    
    # Run fastqc
    fastqc \
        -o ./ \
        -t 4 \
        --extract \
        !{r1_file}

    rm ./${r1_base}_fastqc.zip
    
    cat .command.log > ${SLURM_JOB_ID}_fastqc_${r1_base}.log

    '''
}
