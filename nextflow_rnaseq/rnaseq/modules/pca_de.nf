/*
 * Define "pca_de" process
 *
 */
process pca_de {
    publishDir params.slurmOut, mode: "copy", pattern: "*pca_de.log"

    input:
    path ("*multiqc.log")

    output:
    path ("*pca_de.log"),         emit: log


    shell:
    '''
    # Create the directories
    mkdir -p !{params.figures}
    mkdir -p !{params.rdata}
    mkdir -p !{params.de}
    mkdir -p !{params.gsea}

    ml purge
    ml GCC/9.3.0 OpenMPI/4.0.3 R/4.2.2

    workDir=$(pwd)
    cd !{params.projDir}

    Rscript ./MainRNAseqAnalysis.R 

    # save the slurm output
    cd ${workDir}
    cat .command.log > ${SLURM_JOB_ID}_pca_de.log
    '''
}

