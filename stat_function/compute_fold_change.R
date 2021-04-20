# @Author: Davide Colombo
# @Date: 19th April, 2021

# @Description: a function that computes the fold change between two groups

# @Arguments: 
# data:           the data frame to extract the data
# ...:            the groups to compare

# @Return:
# fc:             fold change vector between groups

compute_fold_change <- function(data, ...) {
      groups <- list(...)
      g1 <- groups[[1]]
      g2 <- groups[[2]]
      fc <- log2(data[, g1]) - log2(data[, g2])
      return (fc)
}
