# @Author: Davide Colombo
# @Date: 19th April, 2021

# @Description: a function that return a boxplot object of the 30 most important
#               gene expressed in each library (based on counts per million)

# @Arguments:
# data:           the data frame to extract the data
# target:         the name of the target variable
# unit:           the unit variable in which data are captured

# @Return:
#                 the boxplot taken with ggplot2 library

library_barplot <- function(data, target, unit) {
      df <- data %>%
                  select("Gene_id", as.name(target))
      ggplot2::ggplot(data = df, mapping = aes(x = df[, 1], y = df[, 2])) +
            geom_col(color = "yellow", fill = "red") +
            coord_flip() + 
            labs(title = paste0("CPM distribution of the m. exp. 30 gene in ", as.character(target)), 
                 x = paste0(as.name(target)),
                 y = paste0(as.character(unit))) +
            theme(axis.text.x = element_text(angle = 60, hjust = 1))
}
