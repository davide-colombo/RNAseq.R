# @Author: Davide Colombo
# @Date: 19th April, 2021

# @Description: a function that computes some basic statistics on the target variable
#               of the data object

# @Arguments: 
# data:           the data frame to extract the data
# target:         the name of the target variable
# unit:           the data measure unit

# @Return:
# df:             a data frame that contains a basic statistic summary of the target
#                 variable

within_Group_Summary <- function(data, target, unit) {
      df <- data.frame(library = target,
                       min = min(data[, target], na.rm = T), 
                       max = max(data[, target], na.rm = T))
      quantiles <- quantile(data[, target], probs = seq(0, 1, .25), na.rm = T)
      quantile_df <- data.frame(q.25 = quantiles[2],
                                q.50 = quantiles[3],
                                mean = mean(data[, target], na.rm = T),
                                q.75 = quantiles[4])
      rownames(quantile_df) <- NULL
      df <- cbind(df, quantile_df)
      df <- df %>%
            relocate(max, .after = everything()) %>%
            mutate(unit = rep(unit, nrow(df))) %>%
            mutate(sd = sd(data[, target], na.rm = T)) %>%
            mutate(N = nrow(data))
      return (df)
}
