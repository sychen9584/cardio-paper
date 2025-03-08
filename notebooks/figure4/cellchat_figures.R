library(CellChat)
library(grid)
library(gridGraphics)

DATA_PATH <- "/home/sychen9584/projects/cardio_paper/data/"
FIGURE_PATH <- "/home/sychen9584/projects/cardio_paper/figures"

cellchat_m3 <- readRDS(file.path(DATA_PATH, "cellchat/m3_cellchat.rds"))
cellchat_m12 <- readRDS(file.path(DATA_PATH, "cellchat/m12_cellchat.rds"))
cellchat_m24 <- readRDS(file.path(DATA_PATH, "cellchat/m24_cellchat.rds"))

cellchat_m3

object.list <- list(m3=cellchat_m3, m12=cellchat_m12, m24=cellchat_m24)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

cell_types <- c(
  "Endo.1", "Endo.3", "Endo.4", "Endo.6", "Endo.7", 
  "Fib.1", "Fib.2", "Fib.3", "Fib.4", "Fib.5", "Fib.6", 
  "Neutrophil", 
  "MC.1", "MC.2", "MC.3", 
  "Smooth Muscle", "T-Cell"
)

# remove cell types with less than 10 cells in at least one sample
cellchat <- subsetCellChat(cellchat, idents.use = cell_types)

sample_colors <- c("m3" = "#98df8a", "m12" = "#FFED6F", "m24" = "#ff9896")
# figure 4a
gg1 <- compareInteractions(cellchat, show.legend = F, group = c("m3", "m12", "m24"), color.use = sample_colors, measure = "count", width = 0.9)
gg2 <- compareInteractions(cellchat, show.legend = F, group = c("m3", "m12", "m24"), measure = "weight", color.use = sample_colors, width = 0.9)
gg1 + gg2
ggsave(file.path(FIGURE_PATH, "figure4/figure4a.png"), width = 6.5, height = 4.2, dpi = 150)


# figure 4b
cell_type_colors <- c(
  
  # Endo types in distinct blue shades
  "Endo.1" = "#2171b5", 
  "Endo.3" = "#6baed6", 
  "Endo.4" = "#9ecae1", 
  "Endo.6" = "#deebf7", 
  "Endo.7" = "#08306b", 
  
  # Fib types in distinct orange shades
  "Fib.1" = "#e6550d", 
  "Fib.2" = "#fd8d3c", 
  "Fib.3" = "#fdae6b", 
  "Fib.4" = "#fdd0a2", 
  "Fib.5" = "#feedde", 
  "Fib.6" = "#a63603", 
  
  "Neutrophil" = "#9467bd",
  # MC types in similar colors
  "MC.1" = "#31a354", 
  "MC.2" = "#74c476", 
  "MC.3" = "#a1d99b", 
  
  # Others
 
  "Smooth Muscle" = "#17becf",
  "T-Cell" = "#aec7e8"
)


png(file.path(FIGURE_PATH, "figure4/figure4b_1.png"), width = 8, height = 8, units = "in", res = 150)
par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0))
netVisual_diffInteraction(cellchat, comparison = c(1, 2), weight.scale = T, measure = "weight", color.use = cell_type_colors, title.name = "m3 vs m12")
dev.off()

png(file.path(FIGURE_PATH, "figure4/figure4b_2.png"), width = 8, height = 8, units = "in", res = 150)
par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0))
netVisual_diffInteraction(cellchat, comparison = c(2, 3), weight.scale = T, measure = "weight", color.use = cell_type_colors, title.name = "m12 vs m24")
dev.off()

# figure 4c
png(file.path(FIGURE_PATH, "figure4/figure4c_1.png"), width = 9, height = 8, units = "in", res = 150)
netVisual_heatmap(cellchat, comparison = c(1, 2), title.name = "3 months vs 12 months" , color.use = cell_type_colors, font.size.title = 14, font.size = 10)
dev.off()

png(file.path(FIGURE_PATH, "figure4/figure4c_2.png"), width = 9, height = 8, units = "in", res = 150)
netVisual_heatmap(cellchat, comparison = c(2, 3), title.name = "12 months vs 24 months",  color.use = cell_type_colors, font.size.title = 14, font.size = 10)
dev.off()

# figure 4d
cell_type_colors["MC/B"] <- "#c7e9c0" 
cell_type_colors["B-Cell"] <- "#1f77b4"

object.list<- list()
object.list[["3 months"]] <- cellchat_m3
object.list[["12 months"]] <- cellchat_m12
object.list[["24 months"]] <- cellchat_m24

generate_cc_pathway_heatmap <- function(object.list, pathway, file_name, figure_path, 
                             heatmap_color = "Reds", font_size_title = 14, font_size = 10, 
                             image_width = 13, image_height = 5, resolution = 150) {
  library(ComplexHeatmap)
  library(grid)
  
  # Generate heatmaps for each object in the list
  ht <- list()
  for (i in seq_along(object.list)) {
    ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = c(pathway),  
                                 title.name = names(object.list)[i], 
                                 font.size.title = font_size_title, font.size = font_size,
                                 color.heatmap = heatmap_color, color.use = cell_type_colors)
  }
  
  # Define full file path
  file_path <- file.path(figure_path, file_name)

  # Save as PNG
  png(file_path, width = image_width, height = image_height, units = "in", res = resolution)
  
  # Draw the combined heatmaps with a clear title
  ComplexHeatmap::draw(ht[[1]] + ht[[2]] + ht[[3]], 
                        ht_gap = unit(0.5, "in"), 
                        column_title = paste(pathway, "Signaling"), 
                        column_title_gp = grid::gpar(fontsize = 16, fontface = "bold"))
  
  dev.off()
  
  message("Heatmap saved to: ", file_path)
}

# IGF signaling
generate_cc_pathway_heatmap(object.list, pathway = "IGF", 
                 file_name = "figure4/figure4d_1.png", 
                 figure_path = FIGURE_PATH)

# TGFb signaling
generate_cc_pathway_heatmap(object.list, pathway = "TGFb", 
                 file_name = "figure4/figure4d_2.png", 
                 figure_path = FIGURE_PATH)

# CCL signaling
generate_cc_pathway_heatmap(object.list, pathway = "CCL", 
                 file_name = "figure4/figure4e_1.png", 
                 figure_path = FIGURE_PATH)

# IL1 signaling
generate_cc_pathway_heatmap(object.list, pathway = "IL1", 
                 file_name = "figure4/figure4e_2.png", 
                 figure_path = FIGURE_PATH)

# PTN signaling
generate_cc_pathway_heatmap(object.list, pathway = "PTN", 
                 file_name = "figure4/figure4f_1.png", 
                 figure_path = FIGURE_PATH)

# MK signaling
generate_cc_pathway_heatmap(object.list, pathway = "MK", 
                 file_name = "figure4/figure4f_2.png", 
                 figure_path = FIGURE_PATH)



