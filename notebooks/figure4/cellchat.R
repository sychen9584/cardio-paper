library(anndata)
library(CellChat)

DATA_PATH <- "/home/sychen9584/projects/cardio_paper/data/"
FIGURE_PATH <- "/home/sychen9584/projects/cardio_paper/figures"


ad <- read_h5ad(file.path(DATA_PATH, "processed/scRNA_m12.h5ad"))
# access count data matrix 
data.input <- t((ad$X))
# access meta data
meta <- ad$obs 
meta$labels <- meta[["cell_type_fine"]] 

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")

CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)


# subset the expression data of signaling genes for saving computation cost
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
options(future.globals.maxSize = 5 * 1024^3)  # Increase to 5 GB
future::plan("multisession", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#> The number of highly variable ligand-receptor pairs used for signaling inference is 692

cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

pathways.show <- c("IGF") 

par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

saveRDS(cellchat, file = file.path(DATA_PATH, "cellchat/m12_cellchat.rds"))
