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
bool_expression <- rna_data[, -1] <= rep(5, ncol(rna_data)-1)
not_exp_index <- which(apply(bool_expression, 1, all) == T)

# Divide the data set into expressed and not expressed gene
not_exp_gene <- rna_data[not_exp_index, ]
exp_gene <- rna_data[-not_exp_index, ]

# ==============================================================================

#                          COUNTS PER MILLION NORMALIZATION
#                             (WITHIN GROUPS VARIABILITY)

# ==============================================================================

source("stat_function/within_Group_Summary.R")
source("plot_function/library_histogram.R")

# Get the count of reads for each library
Nk <- colSums(rna_data[, -1])

# Counts per million normalization
cpm_library <- as.data.frame(t(apply(rna_data[, -1], 1, '/', Nk)))
cpm_rna <- cpm_library %>%
            mutate(Gene_id = rna_data$Gene_id) %>%
            relocate(Gene_id, .before = everything())

# Extract the name of the libraries (samples)
library_names <- colnames(cpm_rna[, -1])

# Compute the summary
cpm_summary_list <- lapply(library_names, within_Group_Summary, data = cpm_rna[, -1], unit = "cpm")

# Transform in data frame
cpm_summary <- NULL
for(i in 1:length(cpm_summary_list)){
      cpm_summary <- rbind(cpm_summary, cpm_summary_list[[i]])
}

# Summary for each library
knitr::kable(cpm_summary, format = "markdown")

# Plot the distribution of the gene expression for each library
histogram_list <- lapply(library_names, library_histogram, data = cpm_rna)
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

# Select the name of the target variables for which we have to compute the box plot
box_target <- colnames(exp_gene[, -1])

# Boxplot list on the original scale of the data
boxlist <- lapply(box_target, compute_BOXPLOT, data = exp_gene, title = "Original scale")
ggpubr::ggarrange(plotlist = boxlist, nrow = 2, ncol = 4)

# Boxplot list after log transform
log2_boxlist <- lapply(box_target, compute_BOXPLOT, data = log2_exp_gene, title = "log2 scale")
ggpubr::ggarrange(plotlist = log2_boxlist, ncol = 4, nrow = 2)










