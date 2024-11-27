library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(SeuratWrappers)
library(monocle3)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(magick)

combine_perturb <- JoinLayers(combine_perturb)

cds <- as.cell_data_set(combine_perturb)

cds <- cluster_cells(cds)

# Get cell metadata
head(colData(cds))

# Get feature/gene metadata
fData(cds)
rownames(fData(cds))[1:10]

fData(cds)$gene_short_name <- rownames(fData(cds))

# Get counts
head(counts(cds))


# Retrieve clustering information from Surat object

# Assign partitions
recreate.partitions <- c(rep(1,length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
recreate.partitions

cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions

# Assign cluster information
list.cluster <- combine_perturb@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster

# Assign UMAP coordinates
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- combine_perturb@reductions$tsne@cell.embeddings


# No trajectory plot
cluster.before.traj <- plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F, 
                                 group_label_size = 5) + theme(legend.position = "right")
cluster.before.traj



# Learn Trajectory
cds <- learn_graph(cds, use_partition = T)


plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F,
           label_branch_points = T, label_roots = T, label_leaves = F,
           group_label_size = 5)

# Order cells in Pseudotime
cds <- order_cells(cds, reduction_method = "UMAP", root_cells = colnames(cds[, clusters(cds) == 0]))
plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = T,
           label_branch_points = T, label_roots = F, label_leaves = F)


# Cells ordered by Monocle3 Pseudotime
head(pseudotime(cds), 10)

cds$monocle3_pseudotime <- pseudotime(cds)

data.pseudo <- as.data.frame(colData(cds))

ggplot(data.pseudo, aes(monocle3_pseudotime, seurat_clusters, fill = seurat_clusters)) + geom_boxplot()

ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(seurat_clusters, monocle3_pseudotime), fill = seurat_clusters)) + geom_boxplot() + xlab("Pseudotime") + ylab("Clusters") + theme_bw()

saveRDS(cds,file = "cds_combine_perturb_pseudotime.rds")

# Find genes that change as a function of pseudotime
deg <- graph_test(cds, neighbor_graph = "principal_graph")
deg %>% arrange(q_value) %>% filter(status == "OK") %>% head()
write.csv(deg,file = "gene_list_change_as_pseudotime.csv")

FeaturePlot(combine_perturb, features = c("NANOG", "SOX1"),reduction = "tsne")

FeaturePlot(combine_perturb, features = c("ZIC2", "MTA2","AHCTF1"),reduction = "tsne")

# Add pseudotime values into the seuratobject
combine_perturb$pseudotime <- pseudotime(cds)
FeaturePlot(combine_perturb, features = "pseudotime",reduction = "tsne")

dotplot <- DotPlot(Ecto_perturb_GBC,features = c("pseudotime"),group.by = "perturb_target",dot.scale = 8) + coord_flip() + RotatedAxis()
dotplot$data$id <- factor(dotplot$data$id,levels =  dotplot$data[order(dotplot$data$avg.exp.scaled),]$id)

dotplot2 <- DotPlot(Ecto_perturb_GBC,features = c("NANOG"),group.by = "perturb_target",dot.scale = 8) + coord_flip() + RotatedAxis()
dotplot2$data$id <- factor(dotplot2$data$id,levels =  dotplot2$data[order(dotplot2$data$avg.exp.scaled),]$id)

my_genes <- row.names(subset(fData(cds), gene_short_name %in% c("POU5F1"))) 
cds_subset <- cds[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_cells_by = "monocle3_pseudotime" )


my_genes <- row.names(subset(fData(cds), gene_short_name %in% c("NANOG"))) 
cds_subset <- cds[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_cells_by = "monocle3_pseudotime" )


my_genes <- row.names(subset(fData(cds), gene_short_name %in% c("SOX1"))) 
cds_subset <- cds[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_cells_by = "monocle3_pseudotime" )


my_genes <- row.names(subset(fData(cds), gene_short_name %in% c("OTX2"))) 
cds_subset <- cds[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_cells_by = "monocle3_pseudotime" )


my_genes <- row.names(subset(fData(cds), gene_short_name %in% c("ZIC2"))) 
cds_subset <- cds[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_cells_by = "monocle3_pseudotime" )


my_genes <- row.names(subset(fData(cds), gene_short_name %in% c("MTA2"))) 
cds_subset <- cds[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_cells_by = "monocle3_pseudotime" )


my_genes <- row.names(subset(fData(cds), gene_short_name %in% c("AHCTF1"))) 
cds_subset <- cds[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_cells_by = "monocle3_pseudotime" )

TF <- read.table(file = "transcription_factor_gene_list.txt",sep = "\t",header = FALSE)
colnames(TF) <- 'Transcription_factors'
NF <- read.csv(file = "Nuclear_protein_list.csv")

my_genes <- row.names(subset(fData(cds), gene_short_name %in% c("POU5F1","NANOG","OTX2","SOX1",intersect(unique(combine_perturb_GBC$perturb_target),c(TF$Transcription_factors,NF$Gene))))) 
pt.matrix <- exprs(cds)[match(my_genes,rownames(rowData(cds))),order(pseudotime(cds))]
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- my_genes;

#K means with 6 groups
htkm <- Heatmap(
  pt.matrix,
  use_raster = TRUE,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  km = 6,
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)


#Ward.D2 Hierarchical Clustering
hthc <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  use_raster = TRUE,
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)

print(htkm)
print(hthc)


