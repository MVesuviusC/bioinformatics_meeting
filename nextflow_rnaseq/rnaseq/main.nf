/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NCH RNASEQ NEXTFLOW CONFIG FILE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INPUT INFORMATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
params.projDir = "/home/gdhpcgroup/yxz006/training/2023/rnaseq/"
params.sampledir = "$params.projDir/misc"
params.samples = "$params.sampledir/MiiiceIiiiinSpaaaaaaaaaaaaace.csv"
params.inputdir = "$params.projDir/input"
params.fastqs = "$params.inputdir/fastqs"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    OUTPUT DIRECTORIES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
params.outdir = "$params.projDir/output"
params.fastqc = "$params.outdir/fastqc"
params.multiqc = "$params.outdir/multiqc"
params.align = "$params.outdir/align"
params.geneCounts = "$params.outdir/geneCounts"
params.figures = "$params.outdir/figures"
params.rdata = "$params.outdir/rdata"
params.de = "$params.outdir/de"
params.gsea = "$params.outdir/gsea"
params.slurmOut = "$params.projDir/slurmOut"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    REFERENCE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
params.ref = "/reference/mus_musculus/mm10/ucsc_assmebly/illumina_download/Sequence/Hisat2Index/genome"
params.genome = "hg38"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    OTHER PARAMETERs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
params.email = ""


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN RNASEQ PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { RNASEQ } from './workflows/rnaseq'

workflow {
    RNASEQ ()
}

