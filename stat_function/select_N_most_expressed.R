# @Author: Davide Colombo
# @Date: 19th April, 2021

# @Description: a function that selects the top N values in the target row of data 
#               object

# @Arguments: 
# data:           the data frame to extract the data
# target:         the name of the target variable
# N:              the number of rows to select

# @Return:
# df:             a data frame that contains the Gene_id and target columns
#                 of the selected N rows

select_N_most_expressed <- function(data, target, N) {
      df <- data %>%
            select("Gene_id", as.name(target)) %>%
            arrange(desc(data[, target]))
      return (df[1:N, ])
}
