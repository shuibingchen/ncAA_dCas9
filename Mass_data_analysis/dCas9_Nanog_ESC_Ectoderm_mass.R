# load the data
library(dplyr)
# install for DEP
# install.packages("rvest")
# install.packages("R.utils")
# install.packages("ncdf4")
# install.packages("mvtnorm")
library(DEP)
library(ggplot2)
library(gplots)

mass_data<-read.csv("Mass_dCas9_Nanog.csv",header = TRUE)

# normalize the intensity to dCas9
# mass_data <- mass_data %>% mutate(across(colnames(mass_data)[7:18], ~./.[Genes=="cas9"]))

# data preparation

# check whether there any duplicated gene names
mass_data$Genes %>% duplicated() %>% any()

# make a table of duplicated gene names
mass_data %>% group_by(Genes) %>%
  summarize(frequency=n()) %>%
  arrange(desc(frequency)) %>%
  filter(frequency>1)

# make unique names using the annotation in the Genes column as primary names and the annotation in "Protein.IDs" as name for those that do not have an gene name.
mass_data_unique<-make_unique(mass_data,"Genes","Protein.Ids",delim=";")

mass_data_unique %>% duplicated() %>% any()

# filtering
Nuclear_protein <- read.csv(file = "Nuclear_factors_list.csv",header = TRUE)
mass_data_unique <- mass_data_unique[mass_data_unique$name %in% Nuclear_protein$Gene,]


# imputed_data <- mass_data_unique %>% mutate_if(is.numeric, ~replace(., is.na(.), min(., na.rm = TRUE)))
# write.csv(imputed_data,file = "imputed_with_min_value_dCas9_ESC_Ectoderm.csv")

imputed_data <- mass_data_unique %>% mutate_if(is.numeric, ~replace(., is.na(.), 0))
write.csv(imputed_data,file = "imputed_with_zero_value_dCas9_ESC_Ectoderm_relative_to_cas9.csv")

imputed_data_filter <- imputed_data %>%
  mutate(ESC_BocK_mean = (ESCNanogBock1 + ESCNanogBock2 + ESCNanogBock3)/3,
         ESC_AbK_mean = (ESCNanogAbk1 + ESCNanogAbk2 + ESCNanogAbk3)/3,
         Ecto_BocK_mean = (EctoNanogBock1 + EctoNanogBock2 + EctoNanogBock3)/3,
         Ecto_AbK_mean = (EctoNanogAbk1 + EctoNanogAbk2 + EctoNanogAbk3)/3) %>%
  filter((ESC_BocK_mean <= ESC_AbK_mean) & (Ecto_BocK_mean <= Ecto_AbK_mean))

# Generate a SummarizedExperiment object using an experimental design
data_columns <- grep("Abk|Bock", colnames(mass_data_unique))

experimental_design <- data.frame(label=c("EctoNanogAbk1","EctoNanogAbk2","EctoNanogAbk3",
                                          "EctoNanogBock1","EctoNanogBock2","EctoNanogBock3",
                                          "ESCNanogAbk1","ESCNanogAbk2","ESCNanogAbk3",
                                          "ESCNanogBock1","ESCNanogBock2","ESCNanogBock3"),
                                  condition=rep(c("EctoAbk","EctoBock","ESCAbk","ESCBock"),each=3),
                                  replicate=rep(c("1","2","3"),times=4))

data_se<-make_se(mass_data_unique,data_columns,experimental_design)

data_se

saveRDS(data_se,file = "data_se_dCas9_Nanog.rds")


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

meanSdPlot(data_filt)
# Visualize normalization by boxplots for all samples before and after normalization
plot_normalization(data_filt, data_norm)


# Impute data for missing values
# Plot a heatmap of proteins with missing values
plot_missval(data_filt)

plot_missval(data_filt[,1:6])

plot_missval(data_filt[,7:12])

# Plot intensity distributions and cumulative fraction of proteins with and without missing values
plot_detect(data_filt)

# norm_data <- as.data.frame(data_norm@assays@data@listData)
# write.csv(norm_data,file = "norm_data_dCas9_ESC_Ectoderm.csv")

# Impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
#  data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)

data_imp <- impute(data_norm, fun = "min")

# Impute missing data using random draws from a manually defined left-shifted Gaussian distribution (for MNAR)
data_imp_man <- impute(data_norm, fun = "man", shift = 1.8, scale = 0.3)

# Impute missing data using the k-nearest neighbour approach (for MAR)
data_imp_knn <- impute(data_norm, fun = "knn", rowmax = 0.9)

# Plot intensity distributions before and after imputation
plot_imputation(data_norm, data_imp)

plot_imputation(data_norm, data_imp_man)

plot_imputation(data_norm, data_imp_knn)


# Differential enrichment analysis
# Differential enrichment analysis  based on linear models and empherical Bayes statistics

# Test manually defined comparisons
data_diff_manual <- test_diff(data_imp, type = "manual", 
                              test = c("EctoAbk_vs_EctoBock",
                                       "ESCAbk_vs_ESCBock",
                                       "ESCAbk_vs_EctoAbk"))


# Denote significant proteins based on user defined cutoffs
dep <- add_rejections(data_diff_manual, alpha = 0.05, lfc = log2(1))

saveRDS(dep,file = "dep_dCas9_Nanog.rds")
# Plot the first and second principal components
plot_pca(dep, x = 1, y = 2, n = 500, point_size = 4)+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

plot_pca(dep[,c(1,2,3,7,8,9)], x = 1, y = 2, n = 500, point_size = 4)+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

# Plot the Pearson correlation matrix
plot_cor(dep, significant = TRUE, lower = 0, upper = 1, pal = "Reds")


# Heatmap of all significant proteins
# Plot a heatmap of all significant proteins with the data centered per protein
plot_heatmap(dep, type = "centered", kmeans = TRUE, 
             k = 6, col_limit = 4, show_row_names = TRUE,
             indicate = c("condition", "replicate"))

plot_heatmap(dep[,c(1,2,3,7,8,9)], type = "centered", kmeans = TRUE, 
             k = 6, col_limit = 4, show_row_names = TRUE,
             indicate = c("condition", "replicate"))

# Plot a heatmap of all significant proteins (rows) and the tested contrasts (columns)
plot_heatmap(dep, type = "contrast", kmeans = TRUE, 
             k = 6, col_limit = 10, show_row_names = FALSE)


# Plot a volcano plot for the contrast "Ubi6 vs Ctrl""
plot_volcano(dep, contrast = "EctoAbk_vs_EctoBock", label_size = 2,add_names = TRUE)

plot_volcano(dep, contrast = "ESCAbk_vs_ESCBock", label_size = 2,add_names = TRUE)


# Plot a barplot for POU5F1 and SOX2
plot_single(dep, proteins = c("POU5F1","SOX2"))


# Plot a barplot for the protein SOX2 with the data centered
plot_single(dep, proteins = "BOLA2B", type = "centered")

plot_single(dep, proteins = "HDAC2", type = "centered")

# Oct4 binds to Nanog promoter region to matain pluripotency
library(ggplot2)
plot_single(dep, proteins = "POU5F1", type = "centered") + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

# Sox2 is thought to be required to drive ectoderm differentiation in pluripotent cells (Thomson et al., 2011)
plot_single(dep, proteins = "SOX2", type = "centered") + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

# Generate a results table
data_results <- get_results(dep)

# Number of significant proteins
data_results %>% filter(significant) %>% nrow()

colnames(data_results)

library(EnhancedVolcano)
EnhancedVolcano(data_results,
                lab = data_results$name,
                x = 'ESCAbk_vs_EctoAbk_ratio',
                y = 'ESCAbk_vs_EctoAbk_p.val',
                pCutoff = 1e-2,
                FCcutoff = 1.5,
                xlim = c(-8,8),
                ylim = c(0,12),
                pointSize = 1,
                labSize = 5,
                colAlpha = 0.9,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                title = NULL,
                subtitle = NULL,
                caption = NULL)

write.csv(data_results,file = "dCas9_Nanog_Abk_vs_BocK_DEA.csv")


# single cell RNA seq gene expression data
data_results <- read.csv(file = "dCas9_Nanog_Abk_vs_BocK_DEA.csv",row.names = 1,header = T)


ESC_specific <- data_results[data_results$ESCAbk_vs_EctoAbk_ratio >= 1 & data_results$ESCAbk_vs_EctoAbk_p.val < 1e-2,]$name

Ecto_specific <- data_results[data_results$ESCAbk_vs_EctoAbk_ratio <= -1 & data_results$ESCAbk_vs_EctoAbk_p.val < 1e-2,]$name
