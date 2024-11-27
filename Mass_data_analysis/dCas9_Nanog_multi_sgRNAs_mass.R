# load the data
library(dplyr)
library(DEP)
library(ggplot2)
library(gplots)

dCas9_nanog_mass<-read.csv("Nanog_dCas9_Abk_Bock_multi_sgRNAs.csv",header = TRUE)

dim(dCas9_nanog_mass)

colnames(dCas9_nanog_mass)

# data preparation

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

imputed_data <- dCas9_nanog_mass_unique %>% mutate_if(is.numeric, ~replace(., is.na(.), 0))
write.csv(imputed_data,file = "imputed_with_zero_value_dCas9_Nanog_multi_sgRNA.csv")

imputed_data_filter <- imputed_data %>%
  mutate(BocK_mean = (Bock_1 + Bock_2 + Bock_3)/3,
         AbK_mean = (Abk_1 + Abk_2 + Abk_3)/3) %>%
  filter((BocK_mean <= AbK_mean) & AbK_mean!=0 )

# Generate a SummarizedExperiment object using an experimental design
data_columns <- grep("Abk|Bock", colnames(dCas9_nanog_mass_unique))

experimental_design <- data.frame(label=c("Bock_1","Bock_2","Bock_3","Abk_1","Abk_2","Abk_3"),
                                  condition=c("Bock","Bock","Bock","Abk","Abk","Abk"),
                                  replicate=c(1,2,3,1,2,3))

data_se<-make_se(dCas9_nanog_mass_unique,data_columns,experimental_design)

data_se

saveRDS(data_se,file = "data_se_dCas9_Nanog_Abk_Bock.rds")

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

plot_numbers(data_filt[,4:6])

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

# Plot intensity distributions and cumulative fraction of proteins with and without missing values
plot_detect(data_filt)


# Impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
data_imp <- impute(data_norm, fun = "min")
# data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)

# Plot intensity distributions before and after imputation
plot_imputation(data_norm, data_imp)


# Differential enrichment analysis
# Differential enrichment analysis  based on linear models and empherical Bayes statistics

# Test every sample versus control
data_diff <- test_diff(data_imp, type = "control", control = "Bock")


# Denote significant proteins based on user defined cutoffs
dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(1))

# Plot the first and second principal components
plot_pca(dep, x = 1, y = 2, n = 500, point_size = 4)+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

# Plot the Pearson correlation matrix
plot_cor(dep, significant = TRUE, lower = 0, upper = 1, pal = "Reds")


# Heatmap of all significant proteins
# Plot a heatmap of all significant proteins with the data centered per protein
plot_heatmap(dep, type = "centered", kmeans = FALSE, 
             k = 6, col_limit = 4, show_row_names = TRUE,
             indicate = c("condition", "replicate"))

# Plot a heatmap of all significant proteins (rows) and the tested contrasts (columns)
plot_heatmap(dep, type = "contrast", kmeans = FALSE, 
             k = 6, col_limit = 10, show_row_names = FALSE)


# Plot a volcano plot for the contrast "Ubi6 vs Ctrl""
plot_volcano(dep, contrast = "Abk_vs_Bock", label_size = 2,add_names = TRUE)


# Plot a barplot for POU5F1 and SOX2
plot_single(dep, proteins = c("POU5F1", "SOX2"))


# Plot a barplot for the protein SOX2 with the data centered
plot_single(dep, proteins = "SOX2", type = "centered")+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

plot_single(dep, proteins = "POU5F1", type = "centered")+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

plot_single(dep, proteins = "ASF1A", type = "centered")+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

plot_single(dep, proteins = c("POU5F1","SOX2"))+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

plot_single(dep, proteins = c("POU5F1","SOX2","ASF1A"))+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

# RNA polemerase components
plot_single(dep, proteins = c("POLR2A","POLR2B","POLR2C","POLR2D","POLR2E","POLR2F","POLR2G","POLR2H","POLR2I","POLR2J","POLR2K","POLR2L"))+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

plot_single(dep, proteins = c("POLR2C","POLR2E","POLR2H","POLR2K"))+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())


# Generate a results table
saveRDS(dep,file = "dep_dCas9_Nanog_Abk_Bock.rds")

data_results <- get_results(dep)

# Number of significant proteins
data_results %>% filter(significant) %>% nrow()

colnames(data_results)

# Generate a wide data.frame
df_wide <- get_df_wide(dep)
# Generate a long data.frame
df_long <- get_df_long(dep)

library(EnhancedVolcano)
EnhancedVolcano(data_results,
                lab = data_results$name,
                x = 'Abk_vs_Bock_ratio',
                y = 'Abk_vs_Bock_p.val',
                pCutoff = 0.001,
                FCcutoff = 1.5,
                xlim = c(-6,6),
                ylim = c(0,8),
                pointSize = 1.5,
                labSize = 5,
                colAlpha = 0.9,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                title = 'dCas9_Nanog_Abk_vs_Bock',
                subtitle = NULL,
                caption = NULL)

