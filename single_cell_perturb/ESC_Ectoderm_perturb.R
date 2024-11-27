library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(EnhancedVolcano)
library(RColorBrewer)
library(pheatmap)
library(corrplot)
library(gplots)

ESC_perturb.data <- Read10X(data.dir = "ESC_perturb/STARsolo/filtered/")

ESC_perturb <- CreateSeuratObject(counts = ESC_perturb.data, project = "ESC_perturb", min.cells = 3, min.features = 200)

ESC_perturb


Ectoderm_perturb.data <- Read10X(data.dir = "Ectoderm_perturb/STARsolo/filtered/")

Ectoderm_perturb <- CreateSeuratObject(counts = Ectoderm_perturb.data, project = "Ectoderm_perturb", min.cells = 3, min.features = 200)

Ectoderm_perturb


combine_perturb <-merge(ESC_perturb, y = Ectoderm_perturb,add.cell.ids = c("ESC","Ecto"),project = "combine_perturb")

combine_perturb

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
combine_perturb[["percent.mt"]] <- PercentageFeatureSet(combine_perturb, pattern = "^MT-")

# Show QC metrics for the first 5 cells
head(combine_perturb@meta.data, 5)


# Visualize QC metrics as a violin plot
VlnPlot(combine_perturb, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(combine_perturb, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(combine_perturb, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

combine_perturb <- subset(combine_perturb, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)

cell_barcodes <- rownames(combine_perturb@meta.data)

names(cell_barcodes) <- colnames(combine_perturb)

combine_perturb <- AddMetaData(
  object = combine_perturb,
  metadata = cell_barcodes,
  col.name = 'cell.barcodes'
)

combine_perturb <- NormalizeData(combine_perturb, normalization.method = "LogNormalize", scale.factor = 10000)

combine_perturb <- FindVariableFeatures(combine_perturb, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(combine_perturb), 10)

plot3 <- VariableFeaturePlot(combine_perturb)
plot4 <- LabelPoints(plot = plot3, points = top10, repel = TRUE)
plot3 + plot4

#Scaling the data
#'regress out' heterogeneity associated with mitochondrial contamination
combine_perturb <- ScaleData(combine_perturb, vars.to.regress = "percent.mt")

#Perform linear dimensional reduction
combine_perturb <- RunPCA(combine_perturb, features = VariableFeatures(object = combine_perturb))

# Examine and visualize PCA results a few different ways
print(combine_perturb[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(combine_perturb, dims = 1:2, reduction = "pca")

DimPlot(combine_perturb, reduction = "pca")

DimHeatmap(combine_perturb, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(combine_perturb, dims = 1:15, cells = 500, balanced = TRUE)


# Determine the 'dimensionality' of the dataset
combine_perturb <- JackStraw(combine_perturb, num.replicate = 100)
combine_perturb <- ScoreJackStraw(combine_perturb, dims = 1:20)

JackStrawPlot(combine_perturb, dims = 1:20,ymax = 0.8)

ElbowPlot(combine_perturb)

#Cluster the cells
combine_perturb <- FindNeighbors(combine_perturb, dims = 1:20)
combine_perturb <- FindClusters(combine_perturb, resolution = 0.2)
# resolution can adjust from 0.2 to 1.5

head(Idents(combine_perturb), 5)

#Run non-linear dimensional reduction (UMAP/tSNE)
combine_perturb <- RunUMAP(combine_perturb, dims = 1:20)

DimPlot(combine_perturb, reduction = "umap")

combine_perturb <- RunTSNE(combine_perturb, dims = 1:20)

DimPlot(combine_perturb, reduction = "tsne")

UMAPPlot(combine_perturb,label=TRUE)

TSNEPlot(combine_perturb,label=TRUE)


#find each cluster characters
VlnPlot(combine_perturb, features = c("POU5F1"))
FeaturePlot(combine_perturb,features = c("POU5F1"))
FeaturePlot(combine_perturb,features = c("POU5F1"),reduction = "tsne")

VlnPlot(combine_perturb, features = c("SOX2"))
FeaturePlot(combine_perturb,features = c("SOX2"))
FeaturePlot(combine_perturb,features = c("SOX2"),reduction = "tsne")

VlnPlot(combine_perturb, features = c("NANOG"))
FeaturePlot(combine_perturb,features = c("NANOG"))
FeaturePlot(combine_perturb,features = c("NANOG"),reduction = "tsne")

VlnPlot(combine_perturb, features = c("SOX1"))
FeaturePlot(combine_perturb,features = c("SOX1"))
FeaturePlot(combine_perturb,features = c("SOX1"),reduction = "tsne")

DotPlot(combine_perturb,features = c("POU5F1","SOX2","NANOG"))

DotPlot(combine_perturb,features = c("POU5F1","NANOG","SOX1","OTX2")) + RotatedAxis()

VlnPlot(combine_perturb, features = c("EGFP"))
FeaturePlot(combine_perturb,features = c("EGFP"))

VlnPlot(combine_perturb, features = c("CAS9"))
FeaturePlot(combine_perturb,features = c("CAS9"))

VlnPlot(combine_perturb, features = c("MCHERRY"))
FeaturePlot(combine_perturb,features = c("MCHERRY"))

# Mesoderm markers
DotPlot(combine_perturb,features = c("MIXL1","TBXT","TBX6","PDGFRA","MESP1","MEOX1","NKX2-5")) + RotatedAxis()
# Ectoderm markers
DotPlot(combine_perturb,features = c("NCAM1","PAX6","SOX1","OTX2","NES","FOXJ3","FABP7","LHX2")) + RotatedAxis()

# Endoderm markers
DotPlot(combine_perturb,features = c("CDH1","NODAL","SMAD2","SMAD3","FOXA2","CDH2","GATA6")) + RotatedAxis()

# Housekeeping genes
DotPlot(combine_perturb,features = c("GAPDH","ACTB","TBP","UBC")) + RotatedAxis()

table(Idents(combine_perturb))

ridge <- RidgePlot(combine_perturb, features = c("POU5F1"), sort = T)
ridge$data$ident <- factor(ridge$data$ident,levels = c("0","3","1","2")) 

VlnPlot(combine_perturb,features = c("POU5F1","NANOG","SOX1","OTX2")) + RotatedAxis()
dotplot <- DotPlot(combine_perturb,features = c("POU5F1","NANOG","SOX1","OTX2")) + RotatedAxis()
dotplot$data$id <- factor(dotplot$data$id,levels =  c('0','3','1','2'))

combine_perturb <- JoinLayers(combine_perturb)

saveRDS(combine_perturb, file = "combine_perturb.rds")

# find markers for every cluster compared to all remaining cells, report only the positive ones
cluster.markers <- FindAllMarkers(combine_perturb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top15 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
DoHeatmap(combine_perturb, features = top15$gene) + NoLegend()


variable_features <- FindMarkers(combine_perturb,ident.1 = 2,ident.2 = 0)

write.csv(variable_features,file = "variable_features_between_Ectoderm_and_ESC.csv")

sig_ESC_features <- rownames(variable_features[variable_features$p_val_adj < 0.001 & variable_features$avg_log2FC > 2,])

sig_Ecto_features <- rownames(variable_features[variable_features$p_val_adj < 0.001 & variable_features$avg_log2FC < -2,])

sig_variable_features <- rownames(variable_features[variable_features$p_val_adj < 0.001 & abs(variable_features$avg_log2FC) > 2,])

DoHeatmap(subset(combine_perturb,downsample = 100), features = sig_variable_features[1:100]) + NoLegend()

# subselect the perturb cells
perturb_infor <- read.csv(file = "ESC_CBC_GBC_reference.csv",row.names = 1)
perturb_infor <- perturb_infor[perturb_infor$cb_num_mismatches == 0 & perturb_infor$fb_num_mismatches == 0,]
perturb_infor <- perturb_infor[,c(2,8)]
# perturb_infor <- perturb_infor %>% dplyr::group_by(cell_barcode) %>% dplyr::count(target) %>% dplyr::filter(n > 4)
UniqueCBC <- perturb_infor %>%
  dplyr::group_by(cell_barcode) %>%
  dplyr::summarise(counts = n()) %>%
  dplyr::filter(counts ==1) 
perturb_infor <- perturb_infor[perturb_infor$cell_barcode %in% UniqueCBC$cell_barcode,]
perturb_infor[perturb_infor$target == 'NonTarget',]$target <- 'NT_ESC'
perturb_cells <- paste("ESC",perturb_infor$cell_barcode,sep = "_")

perturb_infor2 <- read.csv(file = "Ectoderm_CBC_GBC_reference.csv",row.names = 1)
perturb_infor2 <- perturb_infor2[perturb_infor2$cb_num_mismatches == 0 & perturb_infor2$fb_num_mismatches == 0,]
perturb_infor2 <- perturb_infor2[,c(2,8)]
# perturb_infor2 <- perturb_infor2 %>% dplyr::group_by(cell_barcode) %>% dplyr::count(target) %>% dplyr::filter(n > 4)
UniqueCBC2 <- perturb_infor2 %>%
  dplyr::group_by(cell_barcode) %>%
  dplyr::summarise(counts = n()) %>%
  dplyr::filter(counts ==1) 
perturb_infor2 <- perturb_infor2[perturb_infor2$cell_barcode %in% UniqueCBC2$cell_barcode,]
perturb_infor2[perturb_infor2$target == 'NonTarget',]$target <- 'NT_Ecto'
perturb_cells2 <- paste("Ecto",perturb_infor2$cell_barcode,sep = "_")

combine_perturb_GBC <- subset(combine_perturb, subset = cell.barcodes  %in% c(perturb_cells,perturb_cells2))

perturb_target1 <- perturb_infor$target
names(perturb_target1) <- paste("ESC",perturb_infor$cell_barcode,sep = "_")

perturb_target2 <- perturb_infor2$target
names(perturb_target2) <- paste("Ecto",perturb_infor2$cell_barcode,sep = "_")

perturb_target <- c(perturb_target1,perturb_target2)

combine_perturb_GBC <- AddMetaData(
  object = combine_perturb_GBC,
  metadata = perturb_target,
  col.name = 'perturb_target'
)

saveRDS(combine_perturb_GBC, file = "combine_perturb_GBC.rds")


# get the correlation matrix
TF <- read.table(file = "transcription_factor_gene_list.txt",sep = "\t",header = FALSE)
colnames(TF) <- 'Transcription_factors'
NF <- read.csv(file = "Nuclear_protein_list.csv")


edata <- AggregateExpression(combine_perturb_GBC, group.by = 'perturb_target')$RNA

edata <- as.data.frame(edata)

edata <- edata[rownames(edata) %in% sig_variable_features,]

edata <- edata[,colnames(edata) %in% c(TF$Transcription_factors,NF$Gene,"NT-ESC","NT-Ecto")]

corr<-cor(edata)

corrplot(corr,method = "color",tl.cex = 0.4,tl.col = "black",is.corr = FALSE)

'%notin%' <- Negate('%in%')

# filter cells with at least 30 cells per perturb
ESC_perturb_GBC <- subset(combine_perturb_GBC, subset = orig.ident == 'ESC_perturb')
saveRDS(ESC_perturb_GBC, file = "ESC_perturb_GBC.rds")

# filter cells with at least 30 cells per perturb
Ecto_perturb_GBC <- subset(combine_perturb_GBC, subset = orig.ident == 'Ectoderm_perturb')
saveRDS(Ecto_perturb_GBC, file = "Ectoderm_perturb_GBC.rds")

# get the correlation matrix of esc perturb
edata <- AverageExpression(ESC_perturb_GBC, group.by = 'perturb_target')$RNA

edata <- as.data.frame(edata)

# edata <- edata[rownames(edata) %in% sig_variable_features,]

edata <- edata[,colnames(edata) %in% c(NF$Gene,TF$Transcription_factors,"NT-ESC")]

hc <- hclust(dist(t(edata)),method = "com")

plot(hc)

colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)
my_palette <- colorRampPalette(c("blue",'white','red'))(n=100)

avr_exp_by_perturb_matrix <- t(edata)

center <- avr_exp_by_perturb_matrix["NT-ESC",]

avr_exp_by_perturb_matrix_scale <- scale(avr_exp_by_perturb_matrix, center = center,scale = T)

pheatmap(avr_exp_by_perturb_matrix_scale,
         col=my_palette,
         show_colnames = FALSE,
         fontsize_row = 10,
         cluster_cols = F,
         cluster_rows = T)

pseudotime_by_target <- FetchData(object = ESC_perturb_GBC, vars = c('perturb_target', 'pseudotime'))

pseudotime_by_target <- pseudotime_by_target[pseudotime_by_target$perturb_target %in% c(NF$Gene,TF$Transcription_factors,"NT_ESC"),]

pseudotime_by_target <- pseudotime_by_target[pseudotime_by_target$perturb_target %in% c(TF$Transcription_factors,"NT_ESC"),]

ggplot(pseudotime_by_target, aes(y=reorder(perturb_target,pseudotime), x=pseudotime)) + 
  geom_boxplot(outlier.shape = NA) + 
  coord_flip() + 
  RotatedAxis() +
  ylab("Perturb Target")+
  xlab("Pseudotime")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.text.x=element_text(angle = 45,vjust = 0.5,hjust = 0.5))

# ESC perturb distribution plot
perturb_target_ratio <- table(ESC_perturb_GBC$perturb_target,ESC_perturb_GBC$seurat_clusters)
perturb_target_ratio <- t(prop.table(perturb_target_ratio,margin = 1))

write.csv(perturb_target_ratio,file = "ESC_perturb_GBC_ratio.csv",row.names = TRUE)

ESC_dist <- read.csv(file = "ESC_perturb_GBC_ratio.csv",row.names = 1)

cluster <- rep(rownames(ESC_dist),each = dim(ESC_dist)[2])
perturb_target <- rep(colnames(ESC_dist),times = dim(ESC_dist)[1])
ratio <- c(t(ESC_dist))

combine_data <- data.frame(cluster,perturb_target,ratio)

# Stacked + percent
dist_plot<- ggplot(combine_data, aes(fill=cluster, y=ratio, x=perturb_target)) + 
  geom_bar(position="fill", stat="identity") +
  theme_bw()
dist_plot + RotatedAxis()


# Ecto perturb distribution plot
perturb_target_ratio <- table(Ecto_perturb_GBC$perturb_target,Ecto_perturb_GBC$seurat_clusters)
perturb_target_ratio <- t(prop.table(perturb_target_ratio,margin = 1))

write.csv(perturb_target_ratio,file = "Ecto_perturb_GBC_ratio.csv",row.names = TRUE)

Ecto_dist <- read.csv(file = "Ecto_perturb_GBC_ratio.csv",row.names = 1)

cluster <- rep(rownames(Ecto_dist),each = dim(Ecto_dist)[2])
perturb_target <- rep(colnames(Ecto_dist),times = dim(Ecto_dist)[1])
ratio <- c(t(Ecto_dist))

combine_data <- data.frame(cluster,perturb_target,ratio)

# Stacked + percent
dist_plot<- ggplot(combine_data, aes(fill=cluster, y=ratio, x=perturb_target)) + 
  geom_bar(position="fill", stat="identity") +
  coord_flip() + 
  RotatedAxis() +
  theme_bw()

# ectoderm perturb
edata_Ecto <- AggregateExpression(Ecto_perturb_GBC, group.by = 'perturb_target')$RNA

edata_Ecto <- as.data.frame(edata_Ecto)

# edata_Ecto <- edata_Ecto[rownames(edata_Ecto) %in% sig_Ecto_features,]

edata_Ecto <- edata_Ecto[,colnames(edata_Ecto) %in% c(NF$Gene,TF$Transcription_factors,"NT-Ecto")]

corr_Ecto<-cor(edata_Ecto)

corrplot(corr_Ecto,method = "color",tl.cex = 0.4,tl.col = "black",is.corr = FALSE)


edata_Ecto <- AverageExpression(Ecto_perturb_GBC, group.by = 'perturb_target')$RNA

edata_Ecto <- as.data.frame(edata_Ecto)

# edata_Ecto <- edata_Ecto[rownames(edata_Ecto) %in% sig_Ecto_features,]

edata_Ecto <- edata_Ecto[,colnames(edata_Ecto) %in% c(NF$Gene,TF$Transcription_factors,'NT-Ecto')]

hc <- hclust(dist(t(edata_Ecto)))

plot(hc)

colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)
my_palette <- colorRampPalette(c("blue",'white','red'))(n=100)

avr_exp_by_perturb_matrix <- t(edata_Ecto)

center <- avr_exp_by_perturb_matrix["NT-Ecto",]

avr_exp_by_perturb_matrix_scale <- scale(avr_exp_by_perturb_matrix, center = center,scale = T)

pheatmap(avr_exp_by_perturb_matrix_scale,
         col=my_palette,
         show_colnames = FALSE,
         fontsize_row = 10,
         cluster_cols = F,
         cluster_rows = T)

pseudotime_by_target <- FetchData(object = Ecto_perturb_GBC, vars = c('perturb_target', 'pseudotime'))

pseudotime_by_target <- pseudotime_by_target[pseudotime_by_target$perturb_target %in% c(NF$Gene,TF$Transcription_factors,"NT_Ecto"),]

pseudotime_by_target <- pseudotime_by_target[pseudotime_by_target$perturb_target %in% c(TF$Transcription_factors,"NT_Ecto"),]

ggplot(pseudotime_by_target, aes(y=reorder(perturb_target,pseudotime), x=pseudotime)) + 
  geom_boxplot(outlier.shape = NA) + 
  coord_flip() + 
  RotatedAxis() +
  ylab("Perturb Target")+
  xlab("Pseudotime")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.text.x=element_text(angle = 45,vjust = 0.5,hjust = 0.5))

# find markers that separate cells in the 'NonTarget' group compared with '' group  (metadata variable 'perturb_target')
dotplot <- DotPlot(subset(ESC_perturb_GBC,subset = perturb_target %in% c(NF$Gene,TF$Transcription_factors,'NT_ESC')),features = c("NANOG"),group.by = "perturb_target") + coord_flip() + RotatedAxis()
dotplot$data$id <- factor(dotplot$data$id,levels =  dotplot$data[order(dotplot$data$avg.exp.scaled),]$id)

dotplot <- DotPlot(subset(ESC_perturb_GBC,subset = perturb_target %in% c(TF$Transcription_factors,'NT_ESC')),features = c("NANOG"),group.by = "perturb_target") + coord_flip() + RotatedAxis()
dotplot$data$id <- factor(dotplot$data$id,levels =  dotplot$data[order(dotplot$data$avg.exp.scaled),]$id)

dotplot <- DotPlot(subset(ESC_perturb_GBC,subset = perturb_target %in% c(NF$Gene,TF$Transcription_factors,'NT_ESC')),features = c("POU5F1"),group.by = "perturb_target") + coord_flip() + RotatedAxis()
dotplot$data$id <- factor(dotplot$data$id,levels =  dotplot$data[order(dotplot$data$avg.exp.scaled),]$id)


dotplot <- DotPlot(subset(Ecto_perturb_GBC,subset = perturb_target %in% c(NF$Gene,TF$Transcription_factors,'NT_Ecto')),features = c("NANOG"),group.by = "perturb_target") + coord_flip() + RotatedAxis()
dotplot$data$id <- factor(dotplot$data$id,levels =  dotplot$data[order(dotplot$data$avg.exp.scaled),]$id)

dotplot <- DotPlot(subset(Ecto_perturb_GBC,subset = perturb_target %in% c(TF$Transcription_factors,'NT_Ecto')),features = c("NANOG"),group.by = "perturb_target") + coord_flip() + RotatedAxis()
dotplot$data$id <- factor(dotplot$data$id,levels =  dotplot$data[order(dotplot$data$avg.exp.scaled),]$id)


dotplot <- DotPlot(subset(combine_perturb_GBC,subset = perturb_target %in% c(TF$Transcription_factors,'NT_ESC',"NT_Ecto")),features = c("NANOG"),group.by = "perturb_target") + coord_flip() + RotatedAxis()
dotplot$data$id <- factor(dotplot$data$id,levels =  dotplot$data[order(dotplot$data$avg.exp.scaled),]$id)

dotplot <- DotPlot(subset(combine_perturb_GBC,subset = perturb_target %in% c(TF$Transcription_factors,'NT_ESC',"NT_Ecto")),features = c("pseudotime"),group.by = "perturb_target") + coord_flip() + RotatedAxis()
dotplot$data$id <- factor(dotplot$data$id,levels =  dotplot$data[order(dotplot$data$avg.exp.scaled),]$id)


VlnPlot(subset(Ecto_perturb_GBC,subset = perturb_target %in% c(TF$Transcription_factors,'NT_Ecto')),features = c("NANOG"),group.by = "perturb_target") + RotatedAxis()
