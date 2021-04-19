# @Author: Davide Colombo
# @Date: 19th April, 2021

# @Description: a project for RNA sequencing analysis

library(dplyr)
library(ggplot2)

# Read the data file
rna_data <- read.table("data/dataset.txt", header = T, sep = "\t")

# Find the indices of the low expressed gene
bool_expression <- rna_data[, -1] <= rep(5, ncol(rna_data)-1)
not_exp_index <- which(apply(bool_expression, 1, all) == T)

# Divide the data set into expressed and not expressed gene
not_exp_gene <- rna_data[not_exp_index, ]
exp_gene <- rna_data[-not_exp_index, ]










