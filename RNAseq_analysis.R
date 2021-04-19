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










