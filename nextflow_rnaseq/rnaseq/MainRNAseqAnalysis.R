## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
    echo = TRUE,
    error = TRUE,
    cache = TRUE
)


## ----libraries, eval=TRUE-----------------------------------------------------
library(edgeR)
library(tidyverse)
theme_set(theme_bw())
theme_update(plot.title = element_text(hjust = 0.5))
library(qvalue)


## for directoryName in \

##     misc \

##     slurmOut \

##     input/fastqs \

##     output/fastqc \

##     output/multiqc \

##     output/align \

##     output/geneCounts \

##     output/figures \

##     output/rdata \

##     output/de \

##     output/gsea

## do

##   if [ ! -d ${directoryName} ]

##   then

##     mkdir -p ${directoryName}

##   fi

## done


## sbatch sbatchCmds/getSRAData.sh


## zcat input/fastqs/SRR12005052.fastq.gz | head -n 40


## fastqc \

##   -o output/fastqc \

##   -t 5 \

##   --extract \

##   input/fastqs/*.fastq.gz


## sbatch sbatchCmds/fastqc.sh


## hisat2-build \

##   -p 8 \

##   /folder/path/to/reference/pandaGenome.fa \

##   /folder/path/to/reference/pandaGenome


## hisat2 --help

## 

## hisat2 \

##         -x /reference/homo_sapiens/hg38/ucsc_assembly/illumina_download/Sequence/HiSat2Index/genome \

##         -1 input/inputFile_R1.fastq.gz \

##         -2 input/inputFile_R2.fastq.gz \

##         -p 8 \

##         --summary-file output/aligned/inputFile_Summary.txt \

##     | samtools fixmate \

##         -m - - \

##     | samtools sort \

##         - \

##     | samtools markdup \

##         -r - - \

##     | samtools view \

##         -f 2 \

##         -b \

##         - \

##     > output/aligned/inputFile.bam

## 

## samtools index \

##   output/aligned/inputFile.bam


## tophat --help

## 

## tophat \

##     -G /data/Indexes/Annotation/mm9.gtf \

##     -p 8 \

##     -o output/aligned \

##     /reference/homo_sapiens/hg38/ucsc_assembly/illumina_download/Sequence/HiSat2Index/genome \

##     input/inputFile_R1.fastq.gz \

##     input/inputFile_R2.fastq.gz


## sbatch sbatchCmds/align_hisat.sh


## ml purge

## ml load IGV/2.8.0-Java-11.0.2

## 

## igv.sh


## featureCounts --help

## 

## featureCounts \

##   -T 10 \

##   -a /reference/homo_sapiens/hg38/ucsc_assembly/illumina_download/Annotation/Genes/genes.gtf \

##   -o output/geneCounts/combined.txt \

##   -p \

##   --countReadPairs \

##   -s 0 \

##   output/aligned/*.bam


## sbatch sbatchCmds/featureCounts.sh


## ml purge

## ml load GCC/9.3.0 \

##         OpenMPI/4.0.3 \

##         MultiQC/1.9-Python-3.8.2

## 

## multiqc \

##     -o output/multiqc \

##     output/


## ----countsData, dependson='actualFeatureCounts', eval=TRUE-------------------
gene_expr <-
    read_tsv("output/geneCounts/geneCountTable.txt",
        show_col_types = FALSE,
        comment = "#"
    ) %>%
    rename_with(
        ~ str_remove(.x, ".+\\/") %>%
            str_remove(".bam"),
        starts_with("SRR")
    ) %>%
    column_to_rownames("Geneid") %>%
    select(starts_with("SRR"))


## ----metadata, eval=TRUE------------------------------------------------------
metadata <-
    read_csv("misc/MiiiceIiiiinSpaaaaaaaaaaaaace.csv",
        show_col_types = FALSE
    ) %>%
    select(
        "Run",
        "CONDITION",
        "Genotype"
    ) %>%
    mutate(group = paste(CONDITION, Genotype, sep = "_")) %>%
    arrange(match(
        colnames(gene_expr),
        Run
    )) %>%
    rename(sample = Run)

# Make sure the metadata and gene expression data are in the same order
# Should be TRUE
all(metadata$sample == colnames(gene_expr))


## ----normAllData, dependson=c('countsData', 'metadata'), eval=TRUE------------
groups <- factor(metadata$group)

dataDGE <- DGEList(counts = gene_expr, group = groups)

cpm(dataDGE) %>%
    apply(1, function(x) {
        sum(x != 0)
    }) %>%
    as_tibble() %>%
    ggplot(aes(x = value)) +
    geom_histogram()

# Check how many genes are in our dataset before filtering
dim(dataDGE$counts)

keptDf <- filterByExpr(dataDGE, min.count = 1)

dataDGE <- dataDGE[keptDf, keep.lib.sizes = FALSE]

# Check how many genes are in our dataset after filtering
# This step gets rid of really lowly expressed genes that can have highly inflated fold changes
dim(dataDGE$counts)

cpm(dataDGE) %>%
    apply(1, function(x) {
        sum(x != 0)
    }) %>%
    as_tibble() %>%
    ggplot(aes(x = value)) +
    geom_histogram()

# Generate normalization factors to scale the data across samples by read depth
dataDGE <- calcNormFactors(dataDGE)

dataDGE <- estimateDisp(dataDGE, model.matrix(~groups))

# Get log transformed normalized counts per million
norm_counts <-
    as.data.frame(cpm(dataDGE, log = TRUE))

# Save the data out to a file
write.table(norm_counts,
    file = "output/de/normCounts.txt",
    col.names = TRUE,
    quote = FALSE,
    sep = "\t"
)


## ----pca, dependson='normAllData', eval=TRUE----------------------------------
norm_counts <-
    read.delim("output/de/normCounts.txt",
        sep = "\t"
    )
message("test")
# Calculate the mean and standard deviation of the log transformed normalized
# counts per million so we can assess heteroscadasticity
mean_sd_df <-
    apply(norm_counts, 1, function(x) {
        c(
            mean = mean(x),
            sd = sd(x)
        )
    }) %>%
    t() %>%
    as.data.frame()

ggplot(
    mean_sd_df,
    aes(
        x = mean,
        y = sd
    )
) +
    geom_point(alpha = 0.1) +
    ggtitle("Mean vs. SD of log transformed normalized counts per million")

# Run PCA, scaling the data so lowly expressed genes can contribute
pca_out <-
    prcomp(t(norm_counts))

# Add metadata to pca results
pc_plus_metadata <-
    pca_out$x[, 1:9] %>%
    as.data.frame() %>%
    rownames_to_column(var = "sample") %>%
    as_tibble() %>%
    left_join(metadata)

# Calculate the percent variance explianed by each PC
variance <-
    (pca_out$sdev)^2 / sum(pca_out$sdev^2) * 100
names(variance) <- colnames(pca_out$x)

plot(as.numeric(variance[1:20]))

# Plot of PC 1 and 2
pca_all <-
    ggplot(
        pc_plus_metadata,
        aes(
            x = PC1,
            y = PC2,
            color = group
        )
    ) +
    geom_point(size = 2) +
    ggrepel::geom_text_repel(aes(label = sample),
        size = 2,
        color = "gray"
    ) +
    labs(
        x = paste0(
            "PC1 (%var. = ",
            round(variance[1], 1),
            "%)"
        ),
        y = paste0(
            "PC2 (%var. = ",
            round(variance[2], 1),
            "%)"
        ),
        title = "PCA of Mutant Mice From Outer Spaaaaaaace"
    ) +
    scale_color_brewer(
        palette = "Paired",
        name = "Group"
    )

ggsave("output/figures/allSamplePCA_1_2.png",
    plot = pca_all,
    width = 8,
    height = 6
)

# Plot of PC 1, 2 and 3
pca_all <-
    ggplot(
        pc_plus_metadata,
        aes(
            x = PC1,
            y = PC2,
            color = group
        )
    ) +
    geom_point(aes(size = PC3)) +
    ggrepel::geom_text_repel(aes(label = sample),
        size = 2,
        color = "gray"
    ) +
    labs(
        x = paste0(
            "PC1 (%var. = ",
            round(variance[1], 1),
            "%)"
        ),
        y = paste0(
            "PC2 (%var. = ",
            round(variance[2], 1),
            "%)"
        ),
        size = paste0(
            "PC3 (%var. = ",
            round(variance[3], 1),
            "%)"
        ),
        title = "PCA of Mutant Mice From Outer Spaaaaaaace"
    ) +
    scale_color_brewer(
        palette = "Paired",
        name = "Group"
    )

ggsave("output/figures/allSamplePCA_1_2_3.png",
    plot = pca_all,
    width = 8,
    height = 6
)

top_9_pcs <-
    pc_plus_metadata %>%
    pivot_longer(
        cols = starts_with("PC"),
        names_to = "pc",
        values_to = "value"
    ) %>%
    ggplot(aes(
        x = pc,
        y = value,
        color = group,
        fill = group
    )) +
    geom_boxplot(color = "gray") +
    scale_fill_brewer(
        palette = "Paired",
        name = "Group"
    ) +
    labs(
        x = "",
        y = "PC value",
        title = "Distribution of top 9 PCs"
    )

ggsave("output/figures/top9pcs.png",
    plot = top_9_pcs,
    width = 8,
    height = 4
)


## ----runDE, dependson=c('countsData', 'metadata'), eval=TRUE------------------
# Filter down the metadata to just the non-mutant mice
sub_metadata <-
    metadata %>%
    filter(Genotype == "WT")

# Use the filtered metadata to subset the gene expression data
sub_gene_expr <-
    gene_expr[, match(
        sub_metadata$sample,
        colnames(gene_expr)
    )]

# Check that the metadata and gene expression data are in the same order
# Should be TRUE
all(sub_metadata$sample == colnames(sub_gene_expr))

groups <- factor(sub_metadata$group)

dataDGE <- DGEList(counts = sub_gene_expr, group = groups)

dim(dataDGE)

keptDf <- filterByExpr(dataDGE, min.count = 1)

dataDGE <- dataDGE[keptDf, keep.lib.sizes = FALSE]

dim(dataDGE$counts)

dataDGE <- calcNormFactors(dataDGE)

dataDGE <- estimateDisp(dataDGE, model.matrix(~groups))

# Run the actual DE analysis
tested <- exactTest(dataDGE)

# Get all the results
results <- topTags(tested, n = Inf)
results$table$Gene <- rownames(results$table)

# Gather the results,metadata, group average expression values and individual values
results_df <-
    results$table %>%
    as_tibble() %>%
    mutate(signif = FDR <= 0.1 & abs(logFC) > log2(1.5)) %>%
    left_join(cpmByGroup(dataDGE) %>%
        as.data.frame() %>%
        rownames_to_column("Gene")) %>%
    left_join(cpm(dataDGE) %>%
        as.data.frame() %>%
        rownames_to_column("Gene"))

write.table(results_df,
    file = "output/de/edgeROut.txt",
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
)


## ----dePlots, dependson='runDE', eval=TRUE------------------------------------
# Histogram of logFC to see if the values are centered around 0
logfc_hist <-
    ggplot(results_df, aes(x = logFC)) +
    geom_vline(xintercept = 0, color = "red") +
    geom_histogram(bins = 100) +
    ggtitle("Histogram of logFC")

ggsave("output/figures/logfc_hist.png",
    plot = logfc_hist,
    width = 8,
    height = 6
)

# Histogram of p-values to get a sense of the distribution
pval_hist <-
    ggplot(results_df, aes(x = PValue)) +
    geom_histogram(bins = 200) +
    ggtitle("Histogram of p-values")

ggsave("output/figures/pval_hist.png",
    plot = pval_hist,
    width = 8,
    height = 6
)

# Histogram of FDR values which should be similar to p-values but shifted left
fdr_hist <-
    ggplot(results_df, aes(x = FDR)) +
    geom_histogram(bins = 200) +
    ggtitle("Histogram of FDR values")

ggsave("output/figures/fdr_hist.png",
    plot = fdr_hist,
    width = 8,
    height = 6
)

# Volcano plot
top_de <-
    results_df %>%
    filter(signif == TRUE) %>%
    arrange(desc(abs(logFC))) %>%
    group_by(logFC > 0) %>%
    slice_head(n = 50) %>%
    ungroup()

volcano <-
    ggplot(
        results_df,
        aes(
            x = logFC,
            y = -log10(PValue),
            color = signif
        )
    ) +
    geom_point(alpha = 0.5) +
    ggrepel::geom_text_repel(
        data = top_de,
        aes(label = Gene),
        size = 1.5,
        color = "gray"
    ) +
    ggtitle("Volcano plot\nFDR <= 0.1 denoted in blue") +
    scale_color_brewer(palette = "Paired")

ggsave("output/figures/volcano.png",
    plot = volcano,
    width = 10,
    height = 8
)

# Scatterplot of group mean expression values (log10 transformed)
sample_scatterplot <-
    ggplot(
        results_df,
        aes(
            x = Flight_WT,
            y = Ground_control_WT,
            color = signif
        )
    ) +
    geom_point(alpha = 0.2) +
    geom_abline(slope = 1, intercept = 0) +
    ggtitle("Scatterplot of mean CPM values") +
    scale_y_log10() +
    scale_x_log10() +
    labs(
        color = "Significant",
        x = "Space mouse",
        y = "Ground control to Major Tom"
    ) +
    scale_color_brewer(palette = "Paired")

ggsave("output/figures/sample_scatterplot.png",
    plot = sample_scatterplot,
    width = 8,
    height = 6
)


## ----gsea, dependson='runDE', eval=TRUE---------------------------------------
kegg_ref <-
    msigdbr::msigdbr(
        species = "Mus musculus",
        category = "C2",
        subcategory = "CP:KEGG"
    ) %>%
    split(x = .$gene_symbol, f = .$gs_name)

gsea_input <-
    results_df %>%
    arrange(desc(logFC)) %>%
    pull(logFC, name = Gene)

gsea_out <-
    fgsea::fgseaMultilevel(kegg_ref,
        gsea_input,
        minSize = 15,
        maxSize = 500,
        nPerm = 1000
    ) %>%
    arrange(desc(NES)) %>%
    filter(padj < 0.05)

write_tsv(gsea_out,
    file = "output/gsea/gsea_out.tsv"
)

head(gsea_out)

my_fav_path <- "KEGG_MISMATCH_REPAIR"
fav_path_index <- grep(my_fav_path, gsea_out$pathway)
fav_path_genes <- gsea_out$leadingEdge[[fav_path_index]]

# Lets look at the expression of genes in this pathway
fav_path_expr <-
    norm_counts[rownames(norm_counts) %in% fav_path_genes, ] %>%
    rownames_to_column("gene") %>%
    pivot_longer(
        cols = starts_with("SRR"),
        names_to = "sample",
        values_to = "cpm"
    ) %>%
    left_join(metadata) %>%
    ggplot(aes(
        x = gene,
        y = cpm,
        fill = group
    )) +
    geom_boxplot() +
    labs(
        fill = "Group",
        x = "",
        y = "CPM",
        title = paste(
            "Expression of genes in",
            my_fav_path,
            "pathway"
        )
    ) +
    scale_fill_brewer(palette = "Paired")

ggsave("output/figures/fav_path_expr.png",
    plot = fav_path_expr,
    width = 10,
    height = 6
)

