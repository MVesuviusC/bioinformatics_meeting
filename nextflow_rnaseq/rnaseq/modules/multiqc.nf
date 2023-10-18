/*
 * Define "multiqc" process
 *
 *
 */
process multiqc {
    publishDir params.slurmOut, mode: "copy", pattern: "*_multiqc.log"

    input:
    path("*.log")

    output:
    path ("*_multiqc.log"), emit: log

    shell:
    '''
    ml purge
    ml load GCC/9.3.0 \
            OpenMPI/4.0.3 \
            MultiQC/1.9-Python-3.8.2
    
    workDir=$(pwd)
    cd !{params.projDir} 
    multiqc \
        -o output/multiqc \
        output/

    cd ${workDir}
    cat .command.log > ${SLURM_JOB_ID}_multiqc.log
    '''
}

