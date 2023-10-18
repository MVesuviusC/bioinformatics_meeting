/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// MODULE: Loaded from modules/
include { getSRAData      } from '../modules/getSRAData'
include { multiqc         } from '../modules/multiqc'
include { align_hisat     } from '../modules/align_hisat'
include { fastqc          } from '../modules/fastqc'
include { featureCounts   } from '../modules/featureCounts'
include { pca_de          } from '../modules/pca_de'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow RNASEQ {

    // Extract SRA numbers from the sample file 
    //   and slpit them into 24 elements in a queue channel.
    ch_samples = Channel.fromPath(params.samples)
    ch_samples.splitCsv ( header:true, sep:',' )
	.map { it.Run }
	.set {sra_array}

    // Download SRA Data in 24 tasks
    getSRAData (sra_array) 

    // Run fastqc on 24 downloaded fastq.gz files
    ch_fastq = getSRAData.out.fastq
    fastqc ( ch_fastq )

    // Align the reads to genome in 24 tasks
    align_hisat (ch_fastq)    

    // Combine the bam files from align_hisat process and counts the reads 
    ch_bam = align_hisat.out.bam.collect()
    featureCounts (ch_bam)

    // Run multiqc on the current output folder
    ch_log = featureCounts.out.log
    multiqc (ch_log)

    // Run PCA analysis and DE analysis
    ch_log_2 = multiqc.out.log
    pca_de (ch_log_2)
}
