#!/usr/bin/env /home/sychen9584/miniconda3/envs/cardio/bin/Rscript
library(zellkonverter)
library(SingleCellExperiment)
library(scDblFinder)

# Read AnnData directly into SingleCellExperiment
sce <- readH5AD("/home/sychen9584/projects/cardio_paper/data/raw/scATAC_m3_s1/raw.h5ad")
counts_matrix <- as(assays(sce)[["X"]], "dgCMatrix")
sce_new <- SingleCellExperiment(list(counts=counts_matrix))

sce <- scDblFinder(sce_new, clusters=TRUE, aggregateFeatures=TRUE, nfeatures=25, processing="normFeatures")
dbl_score <- sce$scDblFinder.score

dbl_results <- data.frame(barcode = colnames(sce), doublet_score = sce$scDblFinder.score, doublet_class = sce$scDblFinder.class)
head(dbl_results)

write.table(dbl_results, file="/home/sychen9584/projects/cardio_paper/data/raw/scATAC_m3_s1/doublet_scores.csv", sep=",", row.names=FALSE, col.names=TRUE, quote=FALSE)
