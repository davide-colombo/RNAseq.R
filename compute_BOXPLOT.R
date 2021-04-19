
compute_BOXPLOT <- function(data, target, title) {
      fig <- ggplot2::ggplot(data = data, mapping = aes(x = data[, target])) +
            geom_boxplot() + 
            labs(title = as.character(title), 
                 x = as.name(target))
      return (fig)
}
