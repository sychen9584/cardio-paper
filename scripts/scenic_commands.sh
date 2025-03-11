#!/bin/bash
nohup pyscenic grn input.loom mm_mgi_tfs.txt -o adj.tsv --num_workers 8 > pyscenic_grn.log 2>&1 &

wget https://your-database-url/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
wget https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl

nohup pyscenic ctx adj.tsv mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather --annotations_fname motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl --mask_dropouts --expression_mtx_fname input.loom --output reg.csv --num_workers 8 > pyscenic_ctx.log 2>&1 &

nohup pyscenic aucell input.loom reg.csv --output output.loom --num_workers 8 > pyscenic_aucell.log 2>&1 &