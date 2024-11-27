# load the data
library(dplyr)
library(DEP)
library(ggplot2)
library(gplots)

dCas9_nanog_mass<-read.csv("dCas9_NANOG_TBP_AbK_BocK.csv",header = TRUE,row.names = 1)

dim(dCas9_nanog_mass)

colnames(dCas9_nanog_mass)

# normalize the intensity to dCas9
# dCas9_nanog_mass <- dCas9_nanog_mass %>% mutate(across(colnames(dCas9_nanog_mass)[6:17], ~./.[Genes=="cas9"]))

# check whether there any duplicated gene names
dCas9_nanog_mass$Genes %>% duplicated() %>% any()


# make a table of duplicated gene names
dCas9_nanog_mass %>% group_by(Genes) %>%
  summarize(frequency=n()) %>%
  arrange(desc(frequency)) %>%
  filter(frequency>1)

# make unique names using the annotation in the Genes column as primary names and the annotation in "Protein.IDs" as name for those that do not have an gene name.
dCas9_nanog_mass_unique<-make_unique(dCas9_nanog_mass,"Genes","Protein.Ids",delim=";")

dCas9_nanog_mass_unique %>% duplicated() %>% any()


# filtering
Nuclear_protein <- read.csv(file = "Nuclear_factors_list.csv",header = TRUE)
dCas9_nanog_mass_unique <- dCas9_nanog_mass_unique[dCas9_nanog_mass_unique$name %in% Nuclear_protein$Gene,]

# imputed_data <- dCas9_nanog_mass_unique %>% mutate_if(is.numeric, ~replace(., is.na(.), min(., na.rm = TRUE)))
# write.csv(imputed_data,file = "imputed_with_min_value_dCas9_Nanog_TBP.csv")

imputed_data <- dCas9_nanog_mass_unique %>% mutate_if(is.numeric, ~replace(., is.na(.), 0))
write.csv(imputed_data,file = "imputed_with_zero_value_dCas9_Nanog_TBP.csv")

Nanog_specific <- imputed_data %>%
  mutate(Nanog_AbK_mean = (NanogAbk_1+NanogAbk_2+NanogAbk_3)/3,
         Nanog_BocK_mean = (NanogBock_1 + NanogBock_2 + NanogBock_3)/3,
         TBP_AbK_mean = (TBPAbk_1 + TBPAbk_2 + TBPAbk_3)/3,
         TBP_BocK_mean = (TBPBock_1 + TBPBock_2 + TBPBock_3)/3) %>%
  filter((Nanog_AbK_mean > Nanog_BocK_mean) | Nanog_BocK_mean == 0)

TBP_specific <- imputed_data %>%
  mutate(Nanog_AbK_mean = (NanogAbk_1 + NanogAbk_2 + NanogAbk_3)/3,
         Nanog_BocK_mean = (NanogBock_1 + NanogBock_2 + NanogBock_3)/3,
         TBP_AbK_mean = (TBPAbk_1 + TBPAbk_2 + TBPAbk_3)/3,
         TBP_BocK_mean = (TBPBock_1 + TBPBock_2 + TBPBock_3)/3) %>%
  filter((TBP_AbK_mean > TBP_BocK_mean) | TBP_BocK_mean == 0)

# Generate a SummarizedExperiment object using an experimental design
dCas9_nanog_mass_unique1 <- dCas9_nanog_mass_unique[dCas9_nanog_mass_unique$name %in% Nanog_specific$name,]
# Generate a SummarizedExperiment object using an experimental design
data_columns1 <- grep("NanogBock|NanogAbk", colnames(dCas9_nanog_mass_unique1))
experimental_design1 <- data.frame(label=c("NanogAbk_1","NanogAbk_2","NanogAbk_3",
                                          "NanogBock_1","NanogBock_2","NanogBock_3"),
                                  condition=rep(c("Nanog_Abk","Nanog_Bock"),each=3),
                                  replicate=rep(c(1,2,3),times=2))
data_se1<-make_se(dCas9_nanog_mass_unique1,data_columns1,experimental_design1)
data_se1
saveRDS(data_se1,file = "data_se1.rds")

# Generate a SummarizedExperiment object using an experimental design
dCas9_nanog_mass_unique2 <- dCas9_nanog_mass_unique[dCas9_nanog_mass_unique$name %in% TBP_specific$name,]
# Generate a SummarizedExperiment object using an experimental design
data_columns2 <- grep("TBPBock|TBPAbk", colnames(dCas9_nanog_mass_unique1))
experimental_design2 <- data.frame(label=c("TBPAbk_1","TBPAbk_2","TBPAbk_3",
                                           "TBPBock_1","TBPBock_2","TBPBock_3"),
                                   condition=rep(c("TBP_Abk","TBP_Bock"),each=3),
                                   replicate=rep(c(1,2,3),times=2))
data_se2<-make_se(dCas9_nanog_mass_unique2,data_columns2,experimental_design2)
data_se2
saveRDS(data_se2,file = "data_se2.rds")


# Plot a barplot of the protein identification overlap between samples
plot_frequency(data_se1)

plot_frequency(data_se2)

# Filter for proteins that are identified in all replicates of at least one condition
data_filt1 <- filter_missval(data_se1, thr = 0)

data_filt2 <- filter_missval(data_se2, thr = 0)

# Plot a barplot of the number of identified proteins per samples
plot_numbers(data_filt1)

plot_numbers(data_filt2)


# Plot a barplot of the protein identification overlap between samples
plot_coverage(data_filt1)

plot_coverage(data_filt2)


# normalization
# Normalize the data
data_norm1 <- normalize_vsn(data_filt1)

data_norm2 <- normalize_vsn(data_filt2)

# Visualize normalization by boxplots for all samples before and after normalization
plot_normalization(data_filt1, data_norm1)

plot_normalization(data_filt2, data_norm2)


# Impute data for missing values
# Plot a heatmap of proteins with missing values
plot_missval(data_filt1)

plot_missval(data_filt2)

# Plot intensity distributions and cumulative fraction of proteins with and without missing values
plot_detect(data_filt1)

plot_detect(data_filt2)

# Impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
data_imp1 <- impute(data_norm1, fun = "min")

data_imp2 <- impute(data_norm2, fun = "min")

# Plot intensity distributions before and after imputation
plot_imputation(data_norm1, data_imp1)

plot_imputation(data_norm2, data_imp2)



# Differential enrichment analysis
# Differential enrichment analysis  based on linear models and empherical Bayes statistics

# Test every sample versus control
data_diff1 <- test_diff(data_imp1, type = "manual", test = c("Nanog_Abk_vs_Nanog_Bock"))

data_diff2 <- test_diff(data_imp2, type = "manual", test = c("TBP_Abk_vs_TBP_Bock"))

# Denote significant proteins based on user defined cutoffs
dep1 <- add_rejections(data_diff1, alpha = 0.05, lfc = log2(1))

dep2 <- add_rejections(data_diff2, alpha = 0.05, lfc = log2(1))

saveRDS(dep1,file = "dep1.rds")

saveRDS(dep2,file = "dep2.rds")

# Plot the first and second principal components
plot_pca(dep1, x = 1, y = 2, n =100, point_size = 4) + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

plot_pca(dep2, x = 1, y = 2, n = 100, point_size = 4) + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

# Plot the Pearson correlation matrix
plot_cor(dep, significant = TRUE, lower = 0, upper = 1, pal = "Reds")


# Heatmap of all significant proteins
# Plot a heatmap of all significant proteins with the data centered per protein
plot_heatmap(dep1, type = "centered", kmeans = TRUE, 
             k = 6, col_limit = 4, show_row_names = TRUE,
             indicate = c("condition", "replicate"))

library(RColorBrewer)
my_palette <- colorRampPalette(c("blue",'white','red'))(n=100)
plot_data <- as.data.frame(dep1@assays@data@listData)
plot_data <- as.matrix(plot_data)
heatmap.2(plot_data[1:40,], col=my_palette,
          scale="row", key=T, keysize=1, symkey=T,
          dendrogram = "col",
          density.info="none", 
          trace="none",
          cexCol=0.5,
          adjCol = c(1,0.5),
          srtCol = 45,
          adjRow = c(0,1),
          cexRow = 0.5)

plot_heatmap(dep2, type = "centered", kmeans = TRUE, 
             k = 6, col_limit = 4, show_row_names = TRUE,
             indicate = c("condition", "replicate"))

# Plot a barplot for POU5F1 and SOX2
plot_single(dep1, proteins = c("POU5F1","SALL4","SOX2")) + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

plot_single(dep1, proteins = c("POLR2A","POLR2B","POLR2C","POLR2D","POLR2E","POLR2F","POLR2I","POLR2J","POLR2H","POLR2K","POLR2L")) + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

plot_single(dep2, proteins = c("POLR2A","POLR2B","POLR2C","POLR2D","POLR2E","POLR2F","POLR2I","POLR2J","POLR2H","POLR2K","POLR2L")) + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

plot_single(dep2, proteins = c("ZNF787","SRSF11","CARHSP1")) + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

# Generate a results table
data_results1 <- get_results(dep1)

# Number of significant proteins
data_results1 %>% filter(significant) %>% nrow()

colnames(data_results1)

# Generate a results table
data_results2 <- get_results(dep2)

# Number of significant proteins
data_results2 %>% filter(significant) %>% nrow()

colnames(data_results2)

library(ggplot2)
library(EnhancedVolcano)
EnhancedVolcano(data_results1,
                lab = data_results1$name,
                selectLab = data_results1[data_results1$Nanog_Abk_vs_Nanog_Bock_ratio > 3.5 | data_results1$Nanog_Abk_vs_Nanog_Bock_p.val < 1e-6,]$name,
                x = 'Nanog_Abk_vs_Nanog_Bock_ratio',
                y = 'Nanog_Abk_vs_Nanog_Bock_p.val',
                pCutoff = 1e-5,
                FCcutoff = 3,
                xlim = c(-5,5),
                ylim = c(0,10),
                pointSize = 1.5,
                labSize = 5,
                shape = 19,
                colAlpha = 0.9,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                title = NULL,
                caption = NULL,
                subtitle = NULL)

EnhancedVolcano(data_results2,
                lab = data_results2$name,
                x = 'TBP_Abk_vs_TBP_Bock_ratio',
                y = 'TBP_Abk_vs_TBP_Bock_p.val',
                pCutoff = 1e-5,
                FCcutoff = 3,
                xlim = c(-5,5),
                ylim = c(0,10),
                pointSize = 1.5,
                labSize = 5,
                shape = 19,
                colAlpha = 0.9,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                title = NULL,
                caption = NULL,
                subtitle = NULL)

data_results1_top <- data_results1[data_results1$Nanog_Abk_vs_Nanog_Bock_ratio > 2 & data_results1$Nanog_Abk_vs_Nanog_Bock_p.val < 1e-6,c("name","Nanog_Abk_vs_Nanog_Bock_ratio","Nanog_Abk_vs_Nanog_Bock_p.val")]

data_results2_sub <- data_results2[data_results2$name %in% data_results1_top$name,c("name","TBP_Abk_vs_TBP_Bock_ratio","TBP_Abk_vs_TBP_Bock_p.val")]

data_merge <- merge(data_results1_top,data_results2_sub, by = 'name',all.x = TRUE)

dat <- data.frame(
  group = rep(c("NANOG", "TBP"), each=26),
  x = c(data_merge$name,data_merge$name),
  y = c(data_merge$Nanog_Abk_vs_Nanog_Bock_ratio,-data_merge$TBP_Abk_vs_TBP_Bock_ratio))

library(ggplot2)
ggplot(dat, aes(x=x,y=y, fill=group)) + 
  geom_bar(stat="identity", position="identity") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + 
  theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=0.5))

