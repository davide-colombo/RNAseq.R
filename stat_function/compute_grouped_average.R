
compute_grouped_average <- function(data, ...) {
      group <- list(...)
      g1 <- group[[1]]
      g2 <- group[[2]]
      fc_data <- cbind(log2(data[, g1]), log2(data[, g2]))
      return (rowMeans(fc_data))
}
