# load the data
library(dplyr)
library(DEP)
library(ggplot2)
library(gplots)

dCas9_oct4_mass<-read.csv("dCas9_Oct4_sgRNA_Abk_Bock.csv",header = TRUE,row.names = 1)

dim(dCas9_oct4_mass)

colnames(dCas9_oct4_mass)

# normalize the intensity to dCas9
# dCas9_oct4_mass <- dCas9_oct4_mass %>% mutate(across(colnames(dCas9_oct4_mass)[6:17], ~./.[Genes=="cas9"]))

# check whether there any duplicated gene names
dCas9_oct4_mass$Genes %>% duplicated() %>% any()


# make a table of duplicated gene names
dCas9_oct4_mass %>% group_by(Genes) %>%
  summarize(frequency=n()) %>%
  arrange(desc(frequency)) %>%
  filter(frequency>1)

# make unique names using the annotation in the Genes column as primary names and the annotation in "Protein.IDs" as name for those that do not have an gene name.
dCas9_oct4_mass_unique<-make_unique(dCas9_oct4_mass,"Genes","Protein.Ids",delim=";")

dCas9_oct4_mass_unique %>% duplicated() %>% any()

# filtering
Nuclear_protein <- read.csv(file = "Nuclear_factors_list.csv",header = TRUE)
dCas9_oct4_mass_unique <- dCas9_oct4_mass_unique[dCas9_oct4_mass_unique$name %in% Nuclear_protein$Gene,]

# imputed_data <- dCas9_oct4_mass_unique %>% mutate_if(is.numeric, ~replace(., is.na(.), min(., na.rm = TRUE)))
# write.csv(imputed_data,file = "imputed_with_min_value_dCas9_OCT4_Naive_Prime.csv")

imputed_data <- dCas9_oct4_mass_unique %>% mutate_if(is.numeric, ~replace(., is.na(.), 0))
write.csv(imputed_data,file = "imputed_with_zero_value_dCas9_OCT4_relative_to_cas9.csv")

imputed_data_filter <- imputed_data %>%
  mutate(Prime_BocK_mean = (Prime1_1 + Prime1_2 + Prime1_3)/3,
         Prime_AbK_mean = (Prime2_1 + Prime2_2 + Prime2_3)/3,
         Naive_BocK_mean = (Naive1_1 + Naive1_2 + Naive1_3)/3,
         Naive_AbK_mean = (Naive2_1 + Naive2_2 + Naive2_3)/3) %>%
  filter((Prime_BocK_mean < Prime_AbK_mean) & (Naive_BocK_mean < Naive_AbK_mean))

# Generate a SummarizedExperiment object using an experimental design
dCas9_oct4_mass_unique <- dCas9_oct4_mass_unique[dCas9_oct4_mass_unique$name %in% imputed_data_filter$name,]
data_columns <- grep("Prime|Naive", colnames(dCas9_oct4_mass_unique))

experimental_design <- data.frame(label=c("Prime1_1","Prime1_2","Prime1_3",
                                          "Prime2_1","Prime2_2","Prime2_3",
                                          "Naive1_1","Naive1_2","Naive1_3",
                                          "Naive2_1","Naive2_2","Naive2_3"),
                                  condition=rep(c("Prime_Bock","Prime_Abk","Naive_Bock","Naive_Abk"),each=3),
                                  replicate=rep(c(1,2,3),times=4))

data_se<-make_se(dCas9_oct4_mass_unique,data_columns,experimental_design)

data_se

saveRDS(data_se,file = "data_se_dCas9_OCT4_all_new.rds")


# filter on missing values

# Plot a barplot of the protein identification overlap between samples
plot_frequency(data_se)

# Filter for proteins that are identified in all replicates of at least one condition
data_filt <- filter_missval(data_se, thr = 0)

# Less stringent filtering:
# Filter for proteins that are identified in 2 out of 3 replicates of at least one condition
data_filt2 <- filter_missval(data_se, thr = 1)


# Plot a barplot of the number of identified proteins per samples
plot_numbers(data_filt)

# Plot a barplot of the protein identification overlap between samples
plot_coverage(data_filt)


# normalization
# Normalize the data
data_norm <- normalize_vsn(data_filt)

# Visualize normalization by boxplots for all samples before and after normalization
plot_normalization(data_filt, data_norm)

meanSdPlot(data_filt)

# Impute data for missing values
# Plot a heatmap of proteins with missing values
plot_missval(data_filt)

# Plot intensity distributions and cumulative fraction of proteins with and without missing values
plot_detect(data_filt)

# norm_data <- as.data.frame(data_norm@assays@data@listData)
# write.csv(norm_data,file = "norm_data_dCas9_OCT4_naive_prime.csv")
# Impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
data_imp <- impute(data_norm, fun = "min")

# Plot intensity distributions before and after imputation
plot_imputation(data_norm, data_imp)

# Differential enrichment analysis
# Differential enrichment analysis  based on linear models and empherical Bayes statistics

# Test every sample versus control
data_diff <- test_diff(data_imp, type = "manual", test = c("Naive_Abk_vs_Naive_Bock",
                                                           "Naive_Abk_vs_Prime_Abk"))

# Denote significant proteins based on user defined cutoffs
dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(1))

# Plot the first and second principal components
plot_pca(dep, x = 1, y = 2, n = 500, point_size = 4)+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

plot_pca(dep[,c(4,5,6,10,11,12)], x = 1, y = 2, n = 300, point_size = 4)+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

# Plot the Pearson correlation matrix
plot_cor(dep, significant = TRUE, lower = 0, upper = 1, pal = "Reds")


# Heatmap of all significant proteins
# Plot a heatmap of all significant proteins with the data centered per protein
plot_heatmap(dep, type = "centered", kmeans = TRUE, 
             k = 6, col_limit = 4, show_row_names = TRUE,
             indicate = c("condition", "replicate"))

plot_heatmap(dep[,c(4,5,6,10,11,12)], type = "centered", kmeans = TRUE, 
             k = 6, col_limit = 4, show_row_names = TRUE,
             indicate = c("condition", "replicate"))


plot_heatmap(dep[,c(4,5,6,10,11,12)], type = "centered", kmeans = TRUE, 
             k = 6, col_limit = 1, show_row_names = TRUE,
             indicate = c("condition", "replicate"))

# Plot a heatmap of all significant proteins (rows) and the tested contrasts (columns)
plot_heatmap(dep, type = "contrast", kmeans = TRUE, 
             k = 6, col_limit = 10, show_row_names = FALSE)


# Plot a volcano plot for the contrast
plot_volcano(dep, contrast = "Naive_Abk_vs_Naive_Bock", label_size = 2,add_names = TRUE)

plot_volcano(dep, contrast = "Naive_Abk_vs_Prime_Abk", label_size = 2,add_names = TRUE)


# Plot a barplot for POU5F1 and SOX2
plot_single(dep, proteins = c("DNMT3L","DNMT3A","DNMT3B","SMARCA4","ASH2L","POU5F1"))+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

plot_single(dep, proteins = c("ASXL2","ZNF8"))+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

plot_single(dep, proteins = c("BRD3","BRD4"))+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

plot_single(dep, proteins = c("SOX2"))+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
plot_single(dep, proteins = c("BRD3"))+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
plot_single(dep, proteins = c("BRD4"))+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
plot_single(dep, proteins = c("CDK9"))+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

plot_single(dep, proteins = c("ZNF8"))+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
plot_single(dep, proteins = c("DNMT3L"))+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
plot_single(dep, proteins = c("SOX15"))+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
plot_single(dep, proteins = c("SMARCA4"))+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

# Plot a barplot for the protein SOX2 with the data centered
plot_single(dep, proteins = "SOX2", type = "centered")  + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
plot_single(dep, proteins = "POU5F1", type = "centered") + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
plot_single(dep, proteins = "SMARCA4", type = "centered") + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
plot_single(dep, proteins = "DNMT3A", type = "centered") + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
plot_single(dep, proteins = "DNMT3B", type = "centered") + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

plot_single(dep, proteins = "DNMT3L", type = "centered") + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
plot_single(dep, proteins = "SOX15", type = "centered") + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
plot_single(dep, proteins = "ZNF8", type = "centered") + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

plot_single(dep[,c(4,5,6,10,11,12)], proteins = "POU5F1", type = "centered") + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

plot_single(dep, proteins = "SALL4", type = "centered") + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
plot_single(dep[,c(4,5,6,10,11,12)], proteins = "SALL4", type = "centered") + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

plot_single(dep, proteins = "ASH2L", type = "centered") + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
plot_single(dep[,c(4,5,6,10,11,12)], proteins = "ASH2L", type = "centered") + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

plot_single(dep, proteins = "SOX2", type = "centered")  + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
plot_single(dep[,c(4,5,6,10,11,12)], proteins = "SOX2", type = "centered") + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

plot_single(dep, proteins = "ZNF8", type = "centered")  + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
plot_single(dep[,c(4,5,6,10,11,12)], proteins = "ZNF8", type = "centered") + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

plot_single(dep, proteins = "DNMT3L", type = "centered")  + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
plot_single(dep[,c(4,5,6,10,11,12)], proteins = "DNMT3L", type = "centered") + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

plot_single(dep, proteins = "SOX15", type = "centered")  + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
plot_single(dep[,c(4,5,6,10,11,12)], proteins = "SOX15", type = "centered") + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

# The positive transcription elongation factor b (P-TEFb) 
# is composed of cyclins T1 or T2 and cyclin-dependent kinase 9 that regulate the elongation phase of transcription by RNA polymerase II
plot_single(dep, proteins = "CDK9", type = "centered")+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
plot_single(dep[,c(4,5,6,10,11,12)], proteins = "CDK9", type = "centered") + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

plot_single(dep, proteins = "CCNT2", type = "centered")+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
plot_single(dep[,c(4,5,6,10,11,12)], proteins = "CCNT2", type = "centered") + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

# The BET Family Member BRD4 Interacts with OCT4 and Regulates Pluripotency Gene Expression
plot_single(dep, proteins = "BRD4", type = "centered")+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
plot_single(dep[,c(4,5,6,10,11,12)], proteins = "BRD4", type = "centered") + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

plot_single(dep, proteins = "BRD3", type = "centered")+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
plot_single(dep[,c(4,5,6,10,11,12)], proteins = "BRD3", type = "centered") + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

plot_single(dep, proteins = "BRD2", type = "centered")+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
plot_single(dep[,c(4,5,6,10,11,12)], proteins = "BRD2", type = "centered") + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

saveRDS(dep,file = "dep_dCas9_OCT4_Naive_Prime.rds")

# Generate a results table
data_results <- get_results(dep)

# Number of significant proteins
data_results %>% filter(significant) %>% nrow()

colnames(data_results)

library(ggplot2)
library(EnhancedVolcano)
EnhancedVolcano(data_results,
                lab = data_results$name,
                x = 'Naive_Abk_vs_Prime_Abk_ratio',
                y = 'Naive_Abk_vs_Prime_Abk_p.val',
                pCutoff = 0.01,
                FCcutoff = 2,
                xlim = c(-5,5),
                ylim = c(0,10),
                pointSize = 1.5,
                labSize = 5,
                shape = 19,
                colAlpha = 1,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                title = NULL,
                caption = NULL,
                subtitle = NULL)

library(plot3D)
scatter3D(x = data_results$Naive_Abk_vs_Prime_Abk_ratio,
          y = - log10(data_results$Naive_Abk_vs_Prime_Abk_p.val),
          z = data_results$Naive_Abk_centered,
          phi = 0, bty = "f",pch = 20, cex = 1, ticktype = "detailed",col = 'black',
          xlab = "log2 FC Naive vs Prime",
          ylab ="-log10 PValue", 
          zlab = "Naive_Abk_centered")

data_results_select <- data_results[data_results$Naive_Abk_vs_Prime_Abk_ratio > 2 & data_results$Naive_Abk_vs_Prime_Abk_p.val < 0.001,]
points3D(x = data_results_select$Naive_Abk_vs_Prime_Abk_ratio,
         y = - log10(data_results_select$Naive_Abk_vs_Prime_Abk_p.val),
         z = data_results_select$Naive_Abk_centered,
         phi = 0,pch = 20, cex = 1,add = TRUE, col = "red")
text3D(x = data_results_select$Naive_Abk_vs_Prime_Abk_ratio,
       y = - log10(data_results_select$Naive_Abk_vs_Prime_Abk_p.val),
       z = data_results_select$Naive_Abk_centered,
       phi = 0, labels = data_results_select$name,add = TRUE, cex = 1.5, col = "red")


data_results_select2 <- data_results[data_results$Naive_Abk_vs_Prime_Abk_ratio > 2 & data_results$Naive_Abk_vs_Prime_Abk_p.val < 0.001 & data_results$Abk_minus_Bock > 2,]
points3D(x = data_results_select2$Naive_Abk_vs_Prime_Abk_ratio,
         y = - log10(data_results_select2$Naive_Abk_vs_Prime_Abk_p.val),
         z = data_results_select2$Abk_minus_Bock,
         phi = 0,pch = 20, cex = 1,add = TRUE,col = "red")
text3D(x = data_results_select2$Naive_Abk_vs_Prime_Abk_ratio,
       y = - log10(data_results_select2$Naive_Abk_vs_Prime_Abk_p.val),
       z = data_results_select2$Abk_minus_Bock,
       phi = 0, labels = data_results_select2$name,add = TRUE, cex = 1.5,col = "red")



data_results_sub <- data_results[,c(1,2,4,6,11)]
norm_data <- as.data.frame(data_imp@assays@data@listData)
norm_data$name <- rownames(norm_data)
norm_data$Abk_minus_Bock <- (norm_data$Naive_Abk_1 + norm_data$Naive_Abk_2 + norm_data$Naive_Abk_3)/3 - (norm_data$Naive_Bock_1 + norm_data$Naive_Bock_2 + norm_data$Naive_Bock_3)/3
data_results <- merge(data_results_sub,norm_data,all.x = TRUE)
library(plot3D)
scatter3D(x = data_results_new$Naive_Abk_vs_Prime_Abk_ratio,
          y = - log10(data_results_new$Naive_Abk_vs_Prime_Abk_p.val),
          z = data_results_new$Abk_minus_Bock,
          phi = 0, bty = "f",pch = 20, cex = 1, ticktype = "detailed",col = 'black',
          xlab = "log2 FC Naive vs Prime",
          ylab ="-log10 PValue", 
          zlab = "AbK minus BocK")

data_results_select <- data_results_new[data_results_new$Naive_Abk_vs_Prime_Abk_ratio > 2 & data_results_new$Naive_Abk_vs_Prime_Abk_p.val < 0.001,]
points3D(x = data_results_select$Naive_Abk_vs_Prime_Abk_ratio,
         y = - log10(data_results_select$Naive_Abk_vs_Prime_Abk_p.val),
         z = data_results_select$Abk_minus_Bock,
         phi = 0,pch = 20, cex = 1,add = TRUE, col = "darkgray")
text3D(x = data_results_select$Naive_Abk_vs_Prime_Abk_ratio,
       y = - log10(data_results_select$Naive_Abk_vs_Prime_Abk_p.val),
       z = data_results_select$Abk_minus_Bock,
       phi = 0, labels = data_results_select$name,add = TRUE, cex = 1.5, col = "darkgray")


data_results_select2 <- data_results_new[data_results_new$Naive_Abk_vs_Prime_Abk_ratio > 2 & data_results_new$Naive_Abk_vs_Prime_Abk_p.val < 0.001 & data_results_new$Abk_minus_Bock > 2,]
points3D(x = data_results_select2$Naive_Abk_vs_Prime_Abk_ratio,
         y = - log10(data_results_select2$Naive_Abk_vs_Prime_Abk_p.val),
         z = data_results_select2$Abk_minus_Bock,
         phi = 0,pch = 20, cex = 1,add = TRUE,col = "red")
text3D(x = data_results_select2$Naive_Abk_vs_Prime_Abk_ratio,
       y = - log10(data_results_select2$Naive_Abk_vs_Prime_Abk_p.val),
       z = data_results_select2$Abk_minus_Bock,
       phi = 0, labels = data_results_select2$name,add = TRUE, cex = 1.5,col = "red")

