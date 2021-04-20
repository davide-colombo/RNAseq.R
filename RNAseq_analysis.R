# @Author: Davide Colombo
# @Date: 19th April, 2021

# @Description: a project for RNA sequencing analysis

library(dplyr)
library(ggplot2)

# Read the data file
rna_data <- read.table("data/dataset.txt", header = T, sep = "\t")

# ==============================================================================

#                                   DOMANDA 2
# filtrare i geni non espressi (considera un gene espresso se #reads > 5 
# in tutti e 8 i soggetti)

# ==============================================================================

# Find the indices of the low expressed gene
bool_exp <- rna_data[, -1] > 5
exp_index <- which(apply(bool_exp, 1, all) == T)

# Divide the data set into expressed and not expressed gene
exp_gene <- rna_data[exp_index, ]

# ==============================================================================

#                          SEQUENCING DEPTH INSPECTION

# ==============================================================================

# Get the count of reads for each library
Nk <- colSums(exp_gene[, -1]) / 10^6
knitr::kable(Nk, format = "markdown")

# Try to understand if there are any significative difference in the number of reads
# obtained in each samples. This is important because we have to take in account
# that the higher the number of reads in a library, the more could be the bias in
# the effective level of expression of gene.

# In practice, if the distribution of FC values computed on normalized data is not
# centered, this means that the normalization represents a bias for our measures. 

# Min and max normalized number of reads for each group
min_max_Nk_norm <- (Nk - min(Nk)) / (max(Nk) - min(Nk))
ggplot(mapping = aes(x = min_max_Nk_norm)) + 
      geom_histogram(color = "blue", alpha = 0.9, fill = "lightblue", 
                     position = "identity", bins = 30) + 
      labs(title = "Distribution of Nk values (Min-Max normalized)",
           x = "Normalized Nk values")

# Mean and standard deviation normalized number of reads for each group
mean_sd_Nk_norm <- (Nk - mean(Nk)) / sd(Nk)
ggplot(mapping = aes(x = mean_sd_Nk_norm)) + 
      geom_histogram(color = "blue", alpha = 0.9, fill = "lightblue", 
                     position = "identity", bins = 30) + 
      labs(title = "Distribution of Nk values (Mean-sd normalized)",
           x = "Normalized Nk values")

# ==============================================================================

#                           LIBRARY COMPOSITION

# ==============================================================================

# It could be that two libraries have the same number of reads but there are
# genes with high level of expression that increase the bias due to the fact that
# a different number of reads per gene with the same total number of reads increases
# the differences in expression, so the expression between samples appears larger.

# In order to use the Counts per million normalization, we might incurr in biases
# because of the library compositions. So before computing the CPM normalization, 
# I decide to remove the gene with 0 reads in at least one sample.

# In this way, the focus stays on the housekeeping genes.

# de.bool <- rna_data[, -1] == 0
# de.index <- which(apply(de.bool, 1, any) == T)
# hk_genes <- rna_data[-de.index, ]
# 
# # Now, compute the normalization based on Counts per million
# 
# hk_Nk <- colSums(hk_genes[, -1]) / 10^6
# knitr::kable(hk_Nk, format = "markdown")
# 
# hk_cpm <- as.data.frame(t(apply(hk_genes[, -1], 1, '/', hk_Nk))) %>%
#             mutate(Gene_id = hk_genes$Gene_id) %>%
#             relocate(Gene_id, .before = everything())

# The distribution of reads does not change very much, this suggest that
# the samples probably do not belongs to different tissues.

# ==============================================================================

#                          COUNTS PER MILLION NORMALIZATION
#                             (WITHIN GROUPS VARIABILITY)

# ==============================================================================

source("stat_function/within_Group_Summary.R")
source("plot_function/library_histogram.R")

# Counts per million normalization
cpm_rna <- as.data.frame(t(apply(exp_gene[, -1], 1, '/', Nk))) %>%
            mutate(Gene_id = exp_gene$Gene_id) %>%
            relocate(Gene_id, .before = everything())

# Extract the name of the libraries (samples)
library_names <- colnames(cpm_rna[, -1])

# Compute the summary
cpm_summary_list <- lapply(library_names, within_Group_Summary, data = cpm_rna[, -1], unit = "Counts per million")

# Transform in data frame
cpm_summary <- NULL
for(i in 1:length(cpm_summary_list)){
      cpm_summary <- rbind(cpm_summary, cpm_summary_list[[i]])
}

# Summary for each library
knitr::kable(cpm_summary, format = "markdown")

# Plot the distribution of the gene expression for each library
histogram_list <- lapply(library_names, library_histogram, data = cpm_rna, unit = "CPM, unormalized", c(0, 1500))
ggpubr::ggarrange(plotlist = histogram_list, nrow = 2, ncol = 4)

# ==============================================================================

#                          CPM BASED GENE RANK EXPRESSION

# ==============================================================================

source("stat_function/select_N_most_expressed.R")
source("stat_function/gene_rank_expression.R")

# Select the 30 most expressed gene in each library
cpm_most_exp <- lapply(library_names, select_N_most_expressed, data = cpm_rna, N = 30)

# Assign a weight based on the rank of the gene in each library
cpm_ranks <- NULL
for(i in 1:length(cpm_most_exp)) {
      cpm_ranks <- gene_rank_expression(data = cpm_most_exp[[i]], ranks = cpm_ranks)
}

# Order the counts in descending order
cpm_ranks <- cpm_ranks[order(-cpm_ranks[, 2]), ]
knitr::kable(cpm_ranks, format = "markdown")

source("plot_function/library_barplot.R")

# Barplot of the 30 most expressed gene in each library
cpm_barplot_list <- lapply(library_names, library_barplot, 
                           data = cpm_rna[which(cpm_rna$Gene_id %in% cpm_ranks$Gene_id), ], 
                           unit = "Counts per million")
ggpubr::ggarrange(plotlist = cpm_barplot_list, nrow = 2, ncol = 4)

# ==============================================================================

#                         FOLD CHANGE BETWEEN LIBRARIES

# ==============================================================================

source("stat_function/compute_fold_change.R")

# Define the labels of the group of variables to compare
fc_comparison <- data.frame(KO = colnames(cpm_rna)[grep("KO", colnames(cpm_rna))], 
                            WT = colnames(cpm_rna)[grep("WT", colnames(cpm_rna))])

# Compute the fold change between groups (only for gene expressed in each samples)
fc_list <- purrr::pmap(fc_comparison, compute_fold_change, data = cpm_rna)

# Creating a suitable data frame for the 'library_histogram' function
cpm_fc <- NULL
for(i in 1:length(fc_list)) {
      cpm_fc <- cbind(cpm_fc, fc_list[[i]])
}
colnames(cpm_fc) <- apply(fc_comparison, 1, function(x) paste0("log2(", x, "/Nk) ", collapse = " - "))
cpm_fc <- as.data.frame(cpm_fc)

# Taking the name of the variables to compare
fc_names <- colnames(cpm_fc)

# ==============================================================================

# Got a problem with this: we have compute the CPM normalization on all the original
# data. If there are some difference in the expression level of genes due to the
# library composition, then the distribution of values will be biased.

# # Find the NaN values and indices
# na.bool <- is.na(cpm_fc)
# na.index <- which(apply(na.bool, 1, any))
# 
# # Remove the NaN values
# filtered_cpm_fc <- cpm_fc[-na.index, ]
# 
# # Find the -inf/inf values and indices
# inf.bool <- t(apply(filtered_cpm_fc, 1, is.infinite))
# inf.index <- which(apply(inf.bool, 1, any))
# 
# # Remove the -inf/inf values
# filtered_cpm_fc <- filtered_cpm_fc[-inf.index, ]

# ==============================================================================

# Compute the summary statistics between groups
fc_group_list <- lapply(fc_names, within_Group_Summary, data = cpm_fc, unit = "Counts per million")

fc_group_summary <- NULL
for(i in 1:length(fc_group_list)) {
      fc_group_summary <- rbind(fc_group_summary, fc_group_list[[i]])
}

# The summary of the distribution of Fold change Counts per million values between groups
knitr::kable(fc_group_summary, format = "markdown")

# Plot the histogram of the fold change between groups
fc_hist_list <- lapply(fc_names, library_histogram, data = cpm_fc, "CPM, normalized", c(-4, 4))
ggpubr::ggarrange(plotlist = fc_hist_list, nrow = 2, ncol = 2)

# ==============================================================================

#                                   DOMANDA 3
# trasformare i valori di espressione dei geni espressi in scala log2

# ==============================================================================

log2_exp_gene <- log2(exp_gene[, -1]) %>%
                  mutate(Gene_id = exp_gene$Gene_id) %>%
                  relocate(Gene_id, .before = everything())

# ==============================================================================

#                                   DOMANDA 4
# disegnare boxplot per ogni campione prima e dopo la “log-trasformazione” per
# testare l’effetto della stessa

# ==============================================================================

source("plot_function/library_boxplot.R")

# Select the name of the target variables for which we have to compute the box plot
box_target <- colnames(exp_gene[, -1])

# Boxplot list on the original scale of the data
boxlist <- lapply(box_target, library_boxplot, data = exp_gene, title = "Original scale")
ggpubr::ggarrange(plotlist = boxlist, nrow = 2, ncol = 4)

# Boxplot list after log transform
log2_boxlist <- lapply(box_target, library_boxplot, data = log2_exp_gene, title = "log2 scale")
ggpubr::ggarrange(plotlist = log2_boxlist, ncol = 4, nrow = 2)










