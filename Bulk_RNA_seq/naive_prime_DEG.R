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



# DEG analysis use EdgeR package
# Put the data into a DGEList object
library(edgeR)

genelist<-rownames(edata)

y<-DGEList(counts=edata,genes=genelist)


# Filtering
countsPerMillion <- cpm(y)
countCheck <- countsPerMillion > 1
keep <- which(rowSums(countCheck) > 1)
y <- y[keep, ]

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


# the DEA result for all the genes
# dea <- lrt$table

y <- estimateDisp(y, design, robust = TRUE)

plotBCV(y)

fit<-glmQLFit(y,design,robust = TRUE)

head(fit$coefficients)

plotQLDisp(fit)

naive_vs_prime<-makeContrasts(naive-prime,levels = design)

res<-glmQLFTest(fit,contrast = naive_vs_prime)

topTags(res)

is.de<-decideTestsDGE(res)

summary(is.de)

plotMD(res,status = is.de)


toptag <- topTags(res, n = nrow(y$genes), p.value = 1)

dea <- toptag$table 

dea <- dea[order(dea$FDR, -abs(dea$logFC), decreasing = FALSE), ]  # sort the table: ascending of FDR then descending of absolute valued of logFC

table_human<-read.table("table_human_index.csv",header = TRUE,sep = ",",row.names = 1)

id <- match(dea$genes,table_human$ensembl_gene_id)

dea$symbol<-table_human$external_gene_name[id]

deg<-dea[dea$PValue<.05,]

write.csv(dea,file = "Naive_prime_DEG_edgeR.csv")

write.csv(deg,file = "Naive_prime_DEG_sig_edgeR.csv")



# Make a basic volcano plot
library(EnhancedVolcano)
EnhancedVolcano(dea,
                lab = dea$symbol,
                selectLab = c("HMGN1","HMGN4","CXXC1","ZNF8","TRMT2A","UVRAG","ABHD17A","NMNAT3","PHLDA3","POLDIP3","PRR5","SRRM1"),
                x='logFC',
                y='PValue',
                title = 'naive vs prime',
                pCutoff = 1e-10,
                FCcutoff = 1,
                pointSize = 1.5,
                labSize = 2.5,
                shape = 19,
                colAlpha = 0.8,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                caption = 'FC cutoff, 1; p-value cutoff, 1e-10')

EnhancedVolcano(dea,
                lab = dea$symbol,
                x='logFC',
                y='PValue',
                title = 'naive vs prime',
                pCutoff = 1e-10,
                FCcutoff = 1,
                pointSize = 1.5,
                labSize = 2.5,
                shape = 19,
                colAlpha = 0.8,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                caption = 'FC cutoff, 1; p-value cutoff, 1e-10')
