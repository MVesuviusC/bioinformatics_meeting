# Code to teach RNAseq analysis

2023 version

This is a demo of an RNAseq analysis.

We are analyzing data from this paper:
"Nrf2 contributes to the weight gain of mice during space travel"
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7479603/

The data are available on GEO and SRA:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152382
https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA639142

This is intended as a demonstration on one method to analyze RNAseq data, but there are many other ways to do it. The code is set up to work on the Franklin cluster here at NCH, so if you want to use it in other contexts you will need to modify the shell scripts.

## Getting started

To run this analysis, you can either copy/paste the commands into the bash or R terminal, or you can run the entire thing at once by starting R, loading in the rmarkdown library with `library(rmarkdown)`, then running `render("MainRNAseqAnalysis.Rmd")`. This will run the entire analysis (including downloading the data) and produce a report in the form of an html file. This will take an hour or two.

After it runs, check out all the intermediate files that were created so you can get a feel for the different outputs.

Good luck!

Matt
