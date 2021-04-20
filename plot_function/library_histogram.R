# @Author: Davide Colombo
# @Date: 19th April, 2021

# @Description: a function that returns the histogram plot of the target variable
#               present in data object

# @Arguments:
# data:           the data frame to extract the data
# target:         the name of the target variable

# @Return:
# df:             the histogram plot take with ggplot2 library

library_histogram <- function(data, target, unit) {
      ggplot2::ggplot(data = data, mapping = aes(x = data[, target])) +
            geom_histogram(aes(y = ..density..), color = "green", alpha = 0.6, 
                           lwd = 0.2, bins = 100, position = "identity") +
            labs(title = paste0("Histogram of ", as.character(target), " library, (", unit, ")"), 
                 x = as.character(target))
}
