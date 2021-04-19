# @Author: Davide Colombo
# @Date: 19th April, 2021

# @Description: a function that returns the histogram plot of the target variable
#               present in data object

# @Arguments:
# data:           the data frame to extract the data
# target:         the name of the target variable

# @Return:
# df:             the histogram plot take with ggplot2 library

library_histogram <- function(data, target) {
      ggplot2::ggplot(data = data, mapping = aes(x = data[, target])) +
            geom_histogram(color = "green", bins = 100) +
            labs(title = paste0("CPM histogram of ", as.character(target), " library"), 
                 x = as.character(target))
}
