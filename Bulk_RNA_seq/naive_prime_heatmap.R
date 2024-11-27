# load the data
#load the phenotype table
pdata <- read.csv(file = "phenotype_table.csv",row.names = 1)
# load the expression data
edata <- read.csv(file = "expression_table.csv",row.names = 1)

edata <- edata[,1:6]
pdata <- pdata[1:6,]

# table for factor/character variables
pdata$group=as.factor(pdata$group)

table(pdata$group)


library(edgeR)

genelist<-rownames(edata)

y<-DGEList(counts=edata,genes=genelist)

# filtering by using filterByExpr
group<-pdata$group

design<-model.matrix(~0+group)

colnames(design)<-levels(group)

design

keep<-filterByExpr(y,design)
table(keep)

y<-y[keep,,keep.lib.sizes=FALSE]

# Normalization
y <- calcNormFactors(y, method="TMM")

y$samples$group <- pdata$group

# heatmap of selected genes
library(RColorBrewer)
library(gplots)
logCPM<-cpm(y,prior.count = 2,log = TRUE)

table_human<-read.table("table_human_index.csv",header = TRUE,sep = ",",row.names = 1)

id<-match(rownames(logCPM),table_human$ensembl_gene_id)

rownames(logCPM)<-table_human$external_gene_name[id]

selectLab <- c("HMGN1","HMGN4","CXXC1","ZNF8","TRMT2A","UVRAG","ABHD17A","NMNAT3","PHLDA3","POLDIP3","PRR5","SRRM1")

index1<-which(rownames(logCPM) %in% selectLab)

logCPM1<-logCPM[index1,]

my_palette <- colorRampPalette(c("blue",'white','red'))(n=100)

heatmap.2(logCPM1, col=my_palette,
          scale="row", key=T, keysize=1, symkey=T,
          dendrogram = "column",
          density.info="none", trace="none",
          cexCol=1,cexRow = 1)


selectLab2 <- c("CD7","IL6ST","DPPA3","DPPA5","KLF5","KHDC1L","DUSP6","FAT3","THY1","STC1","KLHL4","ZDHHC22","PTPRZ1","CYTL1","SOX11")

index2<-which(rownames(logCPM) %in% selectLab2)

logCPM2<-logCPM[index2,]

my_palette <- colorRampPalette(c("blue",'white','red'))(n=100)

heatmap.2(logCPM2, col=my_palette,
          scale="row", key=T, keysize=1, symkey=T,
          dendrogram = "column",
          density.info="none", trace="none",
          cexCol=1,cexRow = 1)


