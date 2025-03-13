#!/bin/bash

# Define base and target directories
BED_DIR="/home/sychen9584/projects/cardio_paper/data/homer_motifs/input"
OUT_DIR="/home/sychen9584/projects/cardio_paper/data/homer_motifs"

INPUT_FILE="$BED_DIR/input.bed"
BACKGROUND_FILE="$BED_DIR/background.bed"
GENOME="mm10"

# Run HOMER motif enrichment analysis
findMotifsGenome.pl "$INPUT_FILE" "$GENOME" "$OUT_DIR" -bg "$BACKGROUND_FILE" -p 4

echo "Finished motif enrichment analysis"