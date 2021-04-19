# @Author: Davide Colombo
# @Date: 19th April, 2021

# @Description: a function that returns the histogram plot of the target variable
#               present in data object

# @Arguments:
# data:           the data frame to extract the data
# target:         the name of the target variable

# @Return:
# df:             the histogram plot take with ggplot2 library

library_boxplot <- function(data, target, title) {
      fig <- ggplot2::ggplot(data = data, mapping = aes(x = data[, target])) +
            geom_boxplot() + 
            labs(title = as.character(title), 
                 x = as.name(target))
      return (fig)
}
