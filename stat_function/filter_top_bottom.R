
filter_top_bottom <- function(data, tresh) {
      
      index_tresh <- round(tresh*nrow(data))
      return (data[index_tresh:(nrow(data)-index_tresh), ])
}