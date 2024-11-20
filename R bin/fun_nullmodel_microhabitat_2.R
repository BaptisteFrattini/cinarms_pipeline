#' Null model
#'
#' @param metadata_data_mean data
#' @return the path to the subseted raw data file
#' @export

fun_null_model_microhab <- function(meta_data){
  
  # meta_data = targets::tar_read(metadata_data)
  library(ggplot2)
  library(tidyr)
  library(tibble)
  library(forcats)
  library(reshape2)
  library(ggplot2)
  library(EcoSimR)
  library(effsize)
  library(dplyr)
  library(canaper)
  # install.packages("effsize")
  meta = read.csv(meta_data[grepl("metadata", meta_data)], header = TRUE)
  data = read.csv(meta_data[!grepl("metadata", meta_data)], header = TRUE)
  data <- data[ , colSums(data) != 0]
  meta_red <- data.frame(prefixe = substr(meta$arms_name,1,5), 
                         set = meta$arms_name, 
                         Orientation = meta$Orientation,
                         Open_Close = meta$o_c)
  
  rownames(data) <- paste0(meta$arms_name,meta$Orientation,meta$o_c,meta$plate_number)
  
  unique_prefixes <- meta_red$set
  
  # Compute series of observed index ####
  ## Create empty tables ####
  jaccard_results <- as.data.frame(matrix(NA, nrow = 96*15, ncol = 1))
  jaccard_results$names <-  rep(sort(unique(meta$arms_name)), each = 96)
  
  turnover_results <- as.data.frame(matrix(NA, nrow = 96*15, ncol = 1))
  turnover_results$names <-  rep(sort(unique(meta$arms_name)), each = 96)
  
  nestedness_results <- as.data.frame(matrix(NA, nrow = 96*15, ncol = 1))
  nestedness_results$names <-  rep(sort(unique(meta$arms_name)), each = 96)
  
  bray_results <- as.data.frame(matrix(NA, nrow = 96*15, ncol = 1))
  bray_results$names <-  rep(sort(unique(meta$arms_name)), each = 96)
  
  unique_prefixes <- unique(meta$arms_name)

  
  for (prefix in unique_prefixes) {
    #prefix = "CINA1A"
    subset_data <- subset(data, meta$arms_name == prefix)
    
    mat.bray <- vegan::vegdist(subset_data, "bray")
    
    B.pair.pa <- betapart::beta.pair(vegan::decostand(subset_data, "pa"), index.family = "jaccard")
    
    mat.turn <- B.pair.pa$beta.jtu
    mat.nest <- B.pair.pa$beta.jne
    mat.jacc <- B.pair.pa$beta.jac
    
    #Jaccard
    
    df.jacc <- melt(as.matrix(mat.jacc), varnames = c("row", "col"))
    df.jacc$row <- as.character(df.jacc$row)
    df.jacc$col <- as.character(df.jacc$col)
    df.jacc <- subset(df.jacc, row != col)
    df.jacc <- df.jacc %>%
      mutate(
        row = as.character(row),  # Convert 'row' to character
        col = as.character(col),  # Convert 'col' to character
        combined = paste(pmin(row, col), pmax(row, col), sep = "_")
      )
    df.jacc <- df.jacc[!duplicated(df.jacc$combined), ]
    df.jacc$row <- substr(df.jacc$row, 7, 8)
    df.jacc$col <- substr(df.jacc$col, 7, 8)
    df.jacc <- subset(df.jacc, row != col)
    
    jaccard_results$V1[jaccard_results$names == prefix] <- df.jacc$value
   
    
    #Turnover
    
    df.turn <- melt(as.matrix(mat.turn), varnames = c("row", "col"))
    df.turn$row <- as.character(df.turn$row)
    df.turn$col <- as.character(df.turn$col)
    df.turn <- subset(df.turn, row != col)
    df.turn <- df.turn %>%
      mutate(
        row = as.character(row),  # Convert 'row' to character
        col = as.character(col),  # Convert 'col' to character
        combined = paste(pmin(row, col), pmax(row, col), sep = "_")
      )
    df.turn <- df.turn[!duplicated(df.turn$combined), ]
    df.turn$row <- substr(df.turn$row, 7, 8)
    df.turn$col <- substr(df.turn$col, 7, 8)
    df.turn <- subset(df.turn, row != col)
    
    turnover_results$V1[turnover_results$names == prefix] <- df.turn$value
    
    # Nestedness
    df.nest <- melt(as.matrix(mat.nest), varnames = c("row", "col"))
    df.nest$row <- as.character(df.nest$row)
    df.nest$col <- as.character(df.nest$col)
    df.nest <- subset(df.nest, row != col)
    df.nest <- df.nest %>%
      mutate(
        row = as.character(row),  # Convert 'row' to character
        col = as.character(col),  # Convert 'col' to character
        combined = paste(pmin(row, col), pmax(row, col), sep = "_")
      )
    df.nest <- df.nest[!duplicated(df.nest$combined), ]
    df.nest$row <- substr(df.nest$row, 7, 8)
    df.nest$col <- substr(df.nest$col, 7, 8)
    df.nest <- subset(df.nest, row != col)
    
    nestedness_results$V1[nestedness_results$names == prefix] <- df.nest$value
    
    # Bray
    df.bray <- melt(as.matrix(mat.bray), varnames = c("row", "col"))
    df.bray$row <- as.character(df.bray$row)
    df.bray$col <- as.character(df.bray$col)
    df.bray <- subset(df.bray, row != col)
    df.bray <- df.bray %>%
      mutate(
        row = as.character(row),  # Convert 'row' to character
        col = as.character(col),  # Convert 'col' to character
        combined = paste(pmin(row, col), pmax(row, col), sep = "_")
      )
    df.bray <- df.bray[!duplicated(df.bray$combined), ]
    df.bray$row <- substr(df.bray$row, 7, 8)
    df.bray$col <- substr(df.bray$col, 7, 8)
    df.bray <- subset(df.bray, row != col)
    
    bray_results$V1[bray_results$names == prefix] <- df.bray$value
    
  }
  
  
  results <- data.frame(arms = jaccard_results$names,
                        jacc = jaccard_results$V1, 
                        turn = turnover_results$V1, 
                        nest = nestedness_results$V1,
                        bray = bray_results$V1)
  
  
  # Compute series of null index ####
  ## Create empty tables ####
  null_jaccard_results <- as.data.frame(matrix(NA, nrow = 96*15*1000, ncol = 1))
  null_jaccard_results$names <-  rep(sort(unique(meta$arms_name)), each = 96*1000)
  null_jaccard_results$iterations <-  rep(c(1:1000), each = 96*15)
  
  null_turnover_results <- as.data.frame(matrix(NA, nrow = 96*15*1000, ncol = 1))
  null_turnover_results$names <-  rep(sort(unique(meta$arms_name)), each = 96*1000)
  null_turnover_results$iterations <-  rep(c(1:1000), each = 96*15)
  
  null_nestedness_results <- as.data.frame(matrix(NA, nrow = 96*15*1000, ncol = 1))
  null_nestedness_results$names <-  rep(sort(unique(meta$arms_name)), each = 96*1000)
  null_nestedness_results$iterations <-  rep(c(1:1000), each = 96*15)
  
  null_bray_results <- as.data.frame(matrix(NA, nrow = 96*15*1000, ncol = 1))
  null_bray_results$names <-  rep(sort(unique(meta$arms_name)), each = 96*1000)
  null_bray_results$iterations <-  rep(c(1:1000), each = 96*15)
  
  unique_prefixes <- unique(meta$arms_name)
  
  for(i in 1:1000) {
    #i = 1
    for (prefix in unique_prefixes) {
      #prefix = "CINA1A"
      #presence absence
      data_pa <- vegan::decostand(data, "pa")
      null_model_pa <- sim9(data_pa, nReps = 1000, metric = "c_score", algo = "sim2")
      null_model_pa <- null_model_pa$Randomized.Data
      #abondance
      data_round <- round(data*10000, digits = 0)
      null_model_ab <- canaper::cpr_rand_comm(data_round, "quasiswap_count", 10000)
      null_model_ab <- null_model_ab/10000
      null_model_ab <- as.data.frame(vegan::decostand(null_model_ab, "total"))*100
      
      #subset
      null_model_pa <- subset(null_model_pa, meta$arms_name == prefix)
      null_model_ab <- subset(null_model_ab, meta$arms_name == prefix)
      #Presence absence
      #Dissimilarity matrix computing
      B.pair.pa <- betapart::beta.pair(null_model_pa, index.family = "jaccard")
      #Decomposing beta-div
      mat.turn <- B.pair.pa$beta.jtu
      mat.nest <- B.pair.pa$beta.jne
      mat.jacc <- B.pair.pa$beta.jac
      #Jaccard
      df.jacc <- melt(as.matrix(mat.jacc), varnames = c("row", "col"))
      df.jacc$row <- as.character(df.jacc$row)
      df.jacc$col <- as.character(df.jacc$col)
      df.jacc <- subset(df.jacc, row != col)
      df.jacc <- df.jacc %>%
        mutate(
          row = as.character(row),  # Convert 'row' to character
          col = as.character(col),  # Convert 'col' to character
          combined = paste(pmin(row, col), pmax(row, col), sep = "_")
        )
      df.jacc <- df.jacc[!duplicated(df.jacc$combined), ]
      df.jacc$row <- substr(df.jacc$row, 7, 8)
      df.jacc$col <- substr(df.jacc$col, 7, 8)
      df.jacc <- subset(df.jacc, row != col)
      
      indices <- null_jaccard_results$names == prefix & null_jaccard_results$iterations == i
      null_jaccard_results$V1[indices] <- df.jacc$value
      
      #Turnover
      df.turn <- melt(as.matrix(mat.turn), varnames = c("row", "col"))
      df.turn$row <- as.character(df.turn$row)
      df.turn$col <- as.character(df.turn$col)
      df.turn <- subset(df.turn, row != col)
      df.turn <- df.turn %>%
        mutate(
          row = as.character(row),  # Convert 'row' to character
          col = as.character(col),  # Convert 'col' to character
          combined = paste(pmin(row, col), pmax(row, col), sep = "_")
        )
      df.turn <- df.turn[!duplicated(df.turn$combined), ]
      df.turn$row <- substr(df.turn$row, 7, 8)
      df.turn$col <- substr(df.turn$col, 7, 8)
      df.turn <- subset(df.turn, row != col)
      
      indices <- null_turnover_results$names == prefix & null_turnover_results$iterations == i
      null_turnover_results$V1[indices] <- df.turn$value
      
      #Nestedness
      df.nest <- melt(as.matrix(mat.nest), varnames = c("row", "col"))
      df.nest$row <- as.character(df.nest$row)
      df.nest$col <- as.character(df.nest$col)
      df.nest <- subset(df.nest, row != col)
      df.nest <- df.nest %>%
        mutate(
          row = as.character(row),  # Convert 'row' to character
          col = as.character(col),  # Convert 'col' to character
          combined = paste(pmin(row, col), pmax(row, col), sep = "_")
        )
      df.nest <- df.nest[!duplicated(df.nest$combined), ]
      df.nest$row <- substr(df.nest$row, 7, 8)
      df.nest$col <- substr(df.nest$col, 7, 8)
      df.nest <- subset(df.nest, row != col)
      
      indices <- null_nestedness_results$names == prefix & null_nestedness_results$iterations == i
      null_nestedness_results$V1[indices] <- df.nest$value
      
      
      #Bray Curtis
      mat.bray <- vegan::vegdist(null_model_ab, "bray")
      
      df.bray <- melt(as.matrix(mat.bray), varnames = c("row", "col"))
      df.bray$row <- as.character(df.bray$row)
      df.bray$col <- as.character(df.bray$col)
      df.bray <- subset(df.bray, row != col)
      df.bray <- df.bray %>%
        mutate(
          row = as.character(row),  # Convert 'row' to character
          col = as.character(col),  # Convert 'col' to character
          combined = paste(pmin(row, col), pmax(row, col), sep = "_")
        )
      df.bray <- df.bray[!duplicated(df.bray$combined), ]
      df.bray$row <- substr(df.bray$row, 7, 8)
      df.bray$col <- substr(df.bray$col, 7, 8)
      df.bray <- subset(df.bray, row != col)
      
      indices <- null_bray_results$names == prefix & null_bray_results$iterations == i
      null_bray_results$V1[indices] <- df.bray$value
      
    }
    
  }
  
  null_results <- data.frame(arms = null_jaccard_results$names,
                             jacc = null_jaccard_results$V1, 
                             turn = null_turnover_results$V1, 
                             nest = null_nestedness_results$V1,
                             bray = null_bray_results$V1)
  
  write.table(null_results, file = "outputs/null_model/null_results_microhab_1000_iter.csv", sep = ";", dec = ",")
  
  ## Compute quantile table ####
  ## Create empty table ####
  
  null_quantiles_tab_jacc <- as.data.frame(matrix(NA, nrow = 15, ncol = 41))
  rownames(null_quantiles_tab_jacc) <-  sort(unique(meta_mean$arms_name))
  
  null_quantiles_tab_turn <- as.data.frame(matrix(NA, nrow = 15, ncol = 41))
  rownames(null_quantiles_tab_turn) <-  sort(unique(meta_mean$arms_name))
  
  null_quantiles_tab_nest <- as.data.frame(matrix(NA, nrow = 15, ncol = 41))
  rownames(null_quantiles_tab_nest) <-  sort(unique(meta_mean$arms_name))
  
  null_quantiles_tab_bray <- as.data.frame(matrix(NA, nrow = 15, ncol = 41))
  rownames(null_quantiles_tab_bray) <-  sort(unique(meta_mean$arms_name))
  
  ## Make the loop ####
  for (prefix in unique_prefixes) {
    #prefix = "RUNA2"
    subset_null <- subset(null_results, null_results$arms == prefix)
    
    
    #jacc # plutot faire des histo
    
    step = 0.025*(max(subset_null$jacc)-min(subset_null$jacc))
    brek = seq(min(subset_null$jacc),max(subset_null$jacc),  step)
    
    a = ggplot(subset_null, aes(x = jacc)) +
      geom_histogram(breaks = brek, fill = "coral", color = "black", alpha = 0.7) +
      labs(title = paste0("Frequency of Jaccard null values for ", prefix), x = "Beta-div index", y = "Frequency") +
      theme_minimal()
    
    quantile_jacc <- quantile(subset_null$jacc, probs = seq(0, 1, 0.025))
    
    null_quantiles_tab_jacc[prefix,] <- quantile_jacc
    
    #nest
    step = 0.025*(max(subset_null$nest)-min(subset_null$nest))
    brek = seq(min(subset_null$nest),max(subset_null$nest),  step)
    
    b = ggplot(subset_null, aes(x = nest)) +
      geom_histogram(breaks = brek, fill = "coral", color = "black", alpha = 0.7) +
      labs(title = paste0("Frequency of Nestedness null values for ", prefix), x = "Beta-div index", y = "Frequency") +
      theme_minimal()
    
    quantile_nest <- quantile(subset_null$nest, probs = seq(0, 1, 0.025))
    
    null_quantiles_tab_nest[prefix,] <- quantile_nest
    
    #turn
    step = 0.025*(max(subset_null$turn)-min(subset_null$turn))
    brek = seq(min(subset_null$turn),max(subset_null$turn),  step)
    
    c = ggplot(subset_null, aes(x = turn)) +
      geom_histogram(breaks = brek, fill = "coral", color = "black", alpha = 0.7) +
      labs(title = paste0("Frequency of Turnover null values for ", prefix), x = "Beta-div index", y = "Frequency") +
      theme_minimal()
    
    quantile_turn <- quantile(subset_null$turn, probs = seq(0, 1, 0.025))
    
    null_quantiles_tab_turn[prefix,] <- quantile_turn
    
    #bray
    step = 0.025*(max(subset_null$bray)-min(subset_null$bray))
    brek = seq(min(subset_null$bray),max(subset_null$bray),  step)
    
    d = ggplot(subset_null, aes(x = bray)) +
      geom_histogram(breaks = brek, fill = "coral", color = "black", alpha = 0.7) +
      labs(title = paste0("Frequency of Bray-Curtis null values for ", prefix), x = "Beta-div index", y = "Frequency") +
      theme_minimal()
    
    quantile_bray <- quantile(subset_null$bray, probs = seq(0, 1, 0.025))
    
    null_quantiles_tab_bray[prefix,] <- quantile_bray
    
  }
  colnames(null_quantiles_tab_jacc) <- names(quantile_jacc)
  colnames(null_quantiles_tab_turn) <- names(quantile_jacc)
  colnames(null_quantiles_tab_nest) <- names(quantile_jacc)
  colnames(null_quantiles_tab_bray) <- names(quantile_jacc)
  
  null_quantiles <- list(jacc = null_quantiles_tab_jacc, 
                         turn = null_quantiles_tab_turn, 
                         nest = null_quantiles_tab_nest, 
                         bray = null_quantiles_tab_bray)
  ##Table of quantiles ####
  # Function to find quantile for each observed value
  find_quantile <- function(observed_value, quantiles) {
    cols <- as.numeric(gsub("%", "", colnames(quantiles)))
    for (i in seq_along(cols)) {
      if (observed_value <= quantiles[, i]) {
        return(colnames(quantiles)[i])
      }
    }
    return("100%")
  }
  
  # Apply the function to compute quantiles for each index
  observed <- results
  result_quantile <- data.frame(arms = result$arms,
                                jacc_quantile = rep(NA, 96*15),
                                turn_quantile = rep(NA, 96*15), 
                                nest_quantile = rep(NA, 96*15),
                                bray_quantile = rep(NA, 96*15))
  result_quantile$jacc_quantile <- sapply(1:nrow(observed), function(i) find_quantile(observed$jacc[i], null_quantiles$jacc[observed$arms[i], , drop=FALSE]))
  result_quantile$turn_quantile <- sapply(1:nrow(observed), function(i) find_quantile(observed$turn[i], null_quantiles$turn[observed$arms[i], , drop=FALSE]))
  result_quantile$nest_quantile <- sapply(1:nrow(observed), function(i) find_quantile(observed$nest[i], null_quantiles$nest[observed$arms[i], , drop=FALSE]))
  result_quantile$bray_quantile <- sapply(1:nrow(observed), function(i) find_quantile(observed$bray[i], null_quantiles$bray[observed$arms[i], , drop=FALSE]))
  
  # Show result
  print(result_quantile)
  
  quantile_levels <- c("0%", "2.5%", "5%", "7.5%", "10%", "12.5%", "15%", "17.5%", "20%", "22.5%", 
                       "25%", "27.5%", "30%", "32.5%", "35%", "37.5%", "40%", "42.5%", "45%", "47.5%",
                       "50%", "52.5%", "55%", "57.5%", "60%", "62.5%", "65%", "67.5%", "70%", "72.5%", 
                       "75%", "77.5%", "80%", "82.5%", "85%", "87.5%", "90%", "92.5%", "95%", "97.5%", 
                       "100%")
  
  # Ensure the "nest" column is a factor with the quantile levels
  result_quantile$nest <- factor(result_quantile$nest, levels = quantile_levels, ordered = TRUE)
  
  # Convert the factor to numeric to use in the plot
  result_quantile$nest_numeric <- as.numeric(result_quantile$nest)
  
  # Create the data frame for plotting nest quantiles
  plot_data <- data.frame(
    Row = 1:nrow(result_quantile),  # X-axis: Rows (Comparisons)
    Quantile = result_quantile$nest_numeric  # Y-axis: Numeric quantile values for 'nest'
  )
  
  plot_data$arms <- result_quantile$arms
  
  write.table(result_quantile, file = "outputs/null_model/Quantiles_1000_iter.csv", sep = ";", dec = ",")
  
  
  ## Create quantile plot (nest) ####
  intra = c("between ARMS of \n the CINA1 batch", "between ARMS of \n the CINA3 batch","between ARMS of \n the CINA2 batch","between ARMS of \n the CINA4 batch","between ARMS of \n the RUNA2 batch")
  
  ## SES ####
  
  null_SES_tab_jacc <- as.data.frame(matrix(NA, nrow = 96*15, ncol = 1))
  null_SES_tab_jacc$names <-  rep(sort(unique(meta_mean$arms_name)), each = 96)
  
  null_SES_tab_turn <- as.data.frame(matrix(NA, nrow = 96*15, ncol = 1))
  null_SES_tab_turn$names <-  rep(sort(unique(meta_mean$arms_name)), each = 96)
  
  null_SES_tab_nest <- as.data.frame(matrix(NA, nrow = 96*15, ncol = 1))
  null_SES_tab_nest$names <-  rep(sort(unique(meta_mean$arms_name)), each = 96)
  
  null_SES_tab_bray <- as.data.frame(matrix(NA, nrow = 96*15, ncol = 1))
  null_SES_tab_bray$names <-  rep(sort(unique(meta_mean$arms_name)), each = 96)
  
  
  for (prefix in unique_prefixes) {
    #prefix = "CINA1"
    subset_null <- subset(null_results, null_results$arms == prefix)
    subset_results <- subset(results, results$arms == prefix)
    
    SES_jacc <- (subset_results$jacc - mean(subset_null$jacc))/sd(subset_null$jacc)
    
    null_SES_tab_jacc$V1[null_SES_tab_jacc$names == prefix] <- SES_jacc
    
    SES_nest <- (subset_results$nest - mean(subset_null$nest))/sd(subset_null$nest)
    
    null_SES_tab_nest$V1[null_SES_tab_nest$names == prefix] <- SES_nest
    
    SES_turn <- (subset_results$turn - mean(subset_null$turn))/sd(subset_null$turn)
    
    null_SES_tab_turn$V1[null_SES_tab_turn$names == prefix] <- SES_turn
    
    SES_bray <- (subset_results$bray - mean(subset_null$bray))/sd(subset_null$bray)
    
    null_SES_tab_bray$V1[null_SES_tab_bray$names == prefix] <- SES_bray
    
  }
  
  SES_results <- data.frame(arms = null_SES_tab_jacc$names,
                            jacc = null_SES_tab_jacc$V1, 
                            turn = null_SES_tab_turn$V1, 
                            nest = null_SES_tab_nest$V1,
                            bray = null_SES_tab_bray$V1)
  
  write.table(SES_results, file = "outputs/null_model/SES_microhab_1000_iter.csv", sep = ";", dec = ",")
  
  ## Plot SES results ####
  
}
