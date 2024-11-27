# load the data
#load the phenotype table
pdata <- read.csv(file = "phenotype_table.csv",row.names = 1)
# load the expression data
edata <- read.csv(file = "expression_table.csv",row.names = 1)

edata <- edata[,-1]

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
library(ggplot2)

logCPM<-cpm(y,prior.count = 2,log = TRUE)

table_human<-read.table("table_human_index.csv",header = TRUE,sep = ",",row.names = 1)

id<-match(rownames(logCPM),table_human$ensembl_gene_id)

rownames(logCPM)<-table_human$external_gene_name[id]

selectLab <- c("NCAM1","PAX6","SOX1","OTX2","NES","FOXJ3","FABP7","LHX2","SOX2","POU5F1","NANOG","ZIC2")

index1<-which(rownames(logCPM) %in% selectLab)

logCPM1<-logCPM[index1,]

my_palette <- colorRampPalette(c("blue",'white','red'))(n=100)

heatmap.2(logCPM1, col=my_palette,
          scale="row", key=T, keysize=1, symkey=T,
          dendrogram = "both",
          density.info="none", 
          trace="none",
          cexCol=1,
          adjCol = c(1,0.5),
          srtCol = 45,
          cexRow = 1)
