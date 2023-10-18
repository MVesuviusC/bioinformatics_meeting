/*
 * Define "featureCounts" process
 *
 */
process featureCounts {
    publishDir params.geneCounts,  mode: "copy", pattern: "*.{txt,summary}"
    publishDir params.slurmOut, mode: "copy", pattern: "*.log"

    input:
    path (bam_files)

    output:
    path ("*.summary"),                 emit: summary
    path ("geneCountTable.txt"),        emit: counts
    path ("*.log"),                     emit: log


    shell:
    '''
    set -e ### stops bash script if line ends with error
    
    ml purge
    ml load GCC/8.3.0 \
            Subread/2.0.0
    
    # Make variable with gtf file here since it's long
    gene_gtf=/reference/mus_musculus/GRCm38/ensembl/release-86/Annotation/Genes/gtf/Mus_musculus.GRCm38.86.gtf
    
    featureCounts \
      -T 10 \
      -a ${gene_gtf} \
      -g gene_name \
      -o ./geneCountTable.txt \
      -s 2 \
      !{bam_files}

    # save the slurm output
    cat .command.log > ${SLURM_JOB_ID}_featureCounts.log
    '''
}

