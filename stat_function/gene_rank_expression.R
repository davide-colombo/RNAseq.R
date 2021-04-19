# @Author: Davide Colombo
# @Date: 19th April, 2021

# @Description: a function that takes the output of the 'select_N_most_expressed'
#               function and assign to each gene a score based on the rank position
#               of the gene in the i-th library

# @Arguments: 
# data:           the data frame to extract the data
# ranks:         the data frame that keeps the rank for each gene

# @Return:
# ranks

gene_rank_expression <- function(data, ranks) {
      weight <- seq(nrow(data), 1, by = -1) / nrow(data)
      
      if(is.null(ranks)) {
            ranks <- data %>%
                        select(Gene_id) %>%
                        mutate(rank_exp = weight)
      } else {
            for(gene_idx in 1:nrow(data)) {
                  gene_pos <- which(ranks[, "Gene_id"] == data[gene_idx, "Gene_id"])
                  if(length(gene_pos) != 0) {
                        ranks[gene_pos, "rank_exp"] <- ranks[gene_pos, "rank_exp"] + weight[gene_idx]
                  }
            }
      }
      
      return (ranks)
}
