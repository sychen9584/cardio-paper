# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Ensure an input file is provided
if (length(args) == 0) {
  stop("Usage: Rscript scDblFinder_script.R <input_h5ad_file>", call. = FALSE)
}

# Load necessary libraries
suppressPackageStartupMessages(library(zellkonverter))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(scDblFinder))

# Read input file path from command line
input_file <- args[1]
cat(input_file, " loaded... \n")

# Extract directory and filename for output
output_dir <- dirname(input_file)
output_file <- file.path(output_dir, "doublet_scores.csv")

# Read AnnData and convert to SingleCellExperiment
sce <- readH5AD(input_file)
counts_matrix <- as(assays(sce)[["X"]], "dgCMatrix")

# Run scDblFinder
sce <- scDblFinder(SingleCellExperiment(list(counts=counts_matrix)), clusters=TRUE, aggregateFeatures=TRUE, nfeatures=25, processing="normFeatures")

# Extract doublet scores and classifications
dbl_results <- data.frame(
  barcode = colnames(sce),
  doublet_score = sce$scDblFinder.score,
  doublet_class = sce$scDblFinder.class
)

# Save results in the same directory as the input file
write.table(dbl_results, file=output_file, sep=",", row.names=FALSE, col.names=TRUE, quote=FALSE)

# Print success message
cat("Doublet scores saved to:", output_file, "\n")
