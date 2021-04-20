
source("stat_function/compute_fold_change.R")
source("stat_function/compute_grouped_average.R")

edgeR_norm <- function(data, exptresh) {
      
      # First step: filtering data and removing genes with 0 expressions
      exp_bool <- data[, -1] > exptresh
      exp_index <- which(apply(exp_bool, 1, all) == T)
      exp_data <- data[exp_index, ]
      
      # Second step: divide for the number of reads per million
      nreads <- colSums(exp_data[, -1]) / 10^6
      exp_data <- as.data.frame(t(apply(exp_data[, -1], 1, '/', nreads))) %>%
                  mutate(Gene_id = exp_data[, 1]) %>%
                  relocate(Gene_id, .before = everything())
      
      # Third step: selecting the reference sample as the sample which 75th quantile
      # is the most similar to the average of the 75th quantiles between all samples
      q75 <- apply(exp_data[, -1], 2, quantile, probs = 0.75)
      avg75 <- mean(q75)
      ref_sample <- which.min(abs(q75 - avg75))
      
      # Fourth step: compute the log fold change between reference sample and all other samples
      ref_vs_other <- list(ref = colnames(exp_data)[ref_sample+1],
                           other = colnames(exp_data)[-c(1, (ref_sample+1))]) %>% 
                           purrr::cross_df()
      fc_list <- purrr::pmap(ref_vs_other, compute_fold_change, data = exp_data)
      
      for(i in 1:length(fc_list)) {
            current <- fc_list[[i]]
            current <- as.data.frame(current) %>%
                        mutate(Gene_id = exp_data[, 1]) %>%
                        relocate(Gene_id, .before = everything())
            
            ordered <- current[order(current[, 2]), ]
            colnames(ordered)[2] <- paste0("log2(", ref_vs_other[i, ], "/Nk) ", collapse = " - ")
            fc_list[[i]] <- ordered
      }
      
      # # Fifth step: remove the -inf values in the log2 fold change computation
      # inf_bool <- t(apply(fc_df, 1, is.infinite))
      # inf_index <- which(apply(inf_bool, 1, any))
      # fc_df <- fc_df[-inf_index, ]
      
      # Sixth step: compute the grouped average of the fold change
      avg_fc_list <- purrr::pmap(ref_vs_other, compute_grouped_average, data = exp_data)
      
      for(i in 1:length(avg_fc_list)) {
            current <- avg_fc_list[[i]]
            current <- as.data.frame(current) %>%
                        mutate(Gene_id = exp_data[, 1]) %>%
                        relocate(Gene_id, .before = everything())
            ordered <- current[order(current[, 2]), ]
            colnames(ordered) <- paste0("Avg(log2(", ref_vs_other[i, ], "/Nk) ", collapse = " - ")
            avg_fc_list[[i]] <- ordered
      }
      
      # Seventh step: filter biased genes and highly/lowly expressed genes
      # fc_tresh <- round(0.3 * nrow(fc_df))
      # fc_filt <- fc_df[(fc_tresh+1):(nrow(fc_df)-fc_tresh), ]
      # 
      # avg_fc_tresh <- round(0.05*nrow(avg_fc_df))
      # avg_fc_filt <- avg_fc_df[(avg_fc_tresh+1):(nrow(avg_fc_df)-avg_fc_tresh), ]
      
      
      res_list <- vector(mode = "list", length = 5)
      res_list[[1]] <- fc_list
      res_list[[2]] <- avg_fc_list
      # res_list[[3]] <- fc_filt
      # res_list[[4]] <- avg_fc_filt
      res_list[[5]] <- exp_data
      
      return (res_list)
}







