
compute_fold_change <- function(data, ...) {
      groups <- list(...)
      g1 <- groups[[1]]
      g2 <- groups[[2]]
      fc <- log2(data[, g1]) - log2(data[, g2])
      return (fc)
}
