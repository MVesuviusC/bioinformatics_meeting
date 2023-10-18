#!/usr/bin/env nextflow 
nextflow.enable.dsl=2
params.input = "/home/gdhpcgroup/yxz006/data/fastqs/SRR12005075.fastq.gz"

workflow { 
	input_ch = Channel.fromPath(params.input) 
	NUM_LINES(input_ch) 
	NUM_LINES.out.view() 
}

process NUM_LINES { 
	input: path read 
	output: stdout 
	script: 
	""" 
	printf '${read}: ' 
	gunzip -c ${read} | wc -l 
	""" 
}
