#' Null model
#'
#' @param metadata_data_mean data
#' @return the path to the subseted raw data file
#' @export

fun_null_model <- function(metadata_data_mean){
  
  # metadata_data_mean = targets::tar_read(mean_metadata_data)
  library(ggpubr)
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
  
  meta_mean <- read.csv(metadata_data_mean[grepl("metadata", metadata_data_mean)], header = TRUE)
  
  df_mean <- read.csv(metadata_data_mean[!grepl("metadata", metadata_data_mean)], header = TRUE)
  rownames(df_mean) <- meta_mean$arms
  

  df_mean_pa <- vegan::decostand(df_mean, "pa")
  rownames(df_mean_pa) <- meta_mean$arms
  
  # Compute beta diversity indices ####
    ## Create empty tables ####
  
  jaccard_results <- as.data.frame(matrix(NA, nrow = 15, ncol = 1))
  jaccard_results$names <-  rep(sort(unique(meta_mean$arms_name)), each = 3)
  
  turnover_results <- as.data.frame(matrix(NA, nrow = 15, ncol = 1))
  turnover_results$names <-  rep(sort(unique(meta_mean$arms_name)), each = 3)
  
  nestedness_results <- as.data.frame(matrix(NA, nrow = 15, ncol = 1))
  nestedness_results$names <-  rep(sort(unique(meta_mean$arms_name)), each = 3)
  
  bray_results <- as.data.frame(matrix(NA, nrow = 15, ncol = 1))
  bray_results$names <-  rep(sort(unique(meta_mean$arms_name)), each = 3)
  
  unique_prefixes <- unique(meta_mean$arms_name)
  
    ## Make the loop ####
  
  for (prefix in unique_prefixes) {
    #prefix = "CINA1"
    subset_data_pa <- subset(df_mean_pa, meta_mean$arms_name == prefix)
    # subset_data_pa <- subset_data_pa[, colSums(subset_data_pa) > 0]
    subset_data_ab <- subset(df_mean, meta_mean$arms_name == prefix)
    # subset_data_ab <- subset_data_ab[, colSums(subset_data_ab) > 0]
    
    obs.B.pair.pa <- betapart::beta.pair(subset_data_pa, index.family = "jaccard")
    obs.mat.turn <- obs.B.pair.pa$beta.jtu
    obs.mat.nest <- obs.B.pair.pa$beta.jne
    obs.mat.jacc <- obs.B.pair.pa$beta.jac
    obs.mat.bray <- vegan::vegdist(subset_data_ab, "bray")
    
    #jacc
    obs.jacc <- as.vector(obs.mat.jacc)
    jaccard_results$V1[jaccard_results$names == prefix] <- obs.jacc
    #turn
    obs.turn <- as.vector(obs.mat.turn)
    turnover_results$V1[turnover_results$names == prefix] <- obs.turn
    #nest
    obs.nest <- as.vector(obs.mat.nest)
    nestedness_results$V1[nestedness_results$names == prefix] <- obs.nest
    #bray
    obs.bray <- as.vector(obs.mat.bray)
    bray_results$V1[bray_results$names == prefix] <- obs.bray
    
    
  }
  
  results <- data.frame(arms = jaccard_results$names,
                        jacc = jaccard_results$V1, 
                        turn = turnover_results$V1, 
                        nest = nestedness_results$V1,
                        bray = bray_results$V1)
  
  # Compute series of null data ####
    ## Create empty tables ####
  null_jaccard_results <- as.data.frame(matrix(NA, nrow = 15*100, ncol = 1))
  null_jaccard_results$names <-  rep(sort(unique(meta_mean$arms_name)), each = 300)
  null_jaccard_results$iterations <-  rep(c(1:100), each = 3)
  
  null_turnover_results <- as.data.frame(matrix(NA, nrow = 15*100, ncol = 1))
  null_turnover_results$names <-  rep(sort(unique(meta_mean$arms_name)), each = 300)
  null_turnover_results$iterations <-  rep(c(1:100), each = 3)
  
  null_nestedness_results <- as.data.frame(matrix(NA, nrow = 15*100, ncol = 1))
  null_nestedness_results$names <-  rep(sort(unique(meta_mean$arms_name)), each = 300)
  null_nestedness_results$iterations <-  rep(c(1:100), each = 3)
  
  null_bray_results <- as.data.frame(matrix(NA, nrow = 15*100, ncol = 1))
  null_bray_results$names <-  rep(sort(unique(meta_mean$arms_name)), each = 300)
  null_bray_results$iterations <-  rep(c(1:100), each = 3)
  
  unique_prefixes <- unique(meta_mean$arms_name)
  
    ## Make the loop ####
  
  for (prefix in unique_prefixes) {
    #prefix = "CINA1"
    
    for(i in 1:100) {
    #i = 1
    subset_data_pa <- subset(df_mean_pa, meta_mean$arms_name == prefix)
    # subset_data_pa <- subset_data_pa[, colSums(subset_data_pa) > 0]
    subset_data_ab <- subset(df_mean, meta_mean$arms_name == prefix)
    # subset_data_ab <- subset_data_ab[, colSums(subset_data_ab) > 0]
      #Presence absence

    null_model_pa <- sim9(subset_data_pa, nReps = 1000, metric = "c_score", algo = "sim2")
    null_model_pa <- null_model_pa$Randomized.Data
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
    indices <- null_nestedness_results$names == prefix & null_nestedness_results$iterations == i
    null_nestedness_results$V1[indices] <- df.nest$value
    
      #abondance
    subset_data_ab <- round(subset_data_ab*10000, digits = 0)
    null_model_ab <- canaper::cpr_rand_comm(subset_data_ab, "quasiswap_count", 10000)
    null_model_ab <- null_model_ab/10000
    null_model_ab <- as.data.frame(vegan::decostand(null_model_ab, "total"))*100
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
    indices <- null_bray_results$names == prefix & null_bray_results$iterations == i
    null_bray_results$V1[indices] <- df.bray$value
    
    }
  
  }

  null_results <- data.frame(arms = null_jaccard_results$names,
                             jacc = null_jaccard_results$V1, 
                             turn = null_turnover_results$V1, 
                             nest = null_nestedness_results$V1,
                             bray = null_bray_results$V1)
  
  # Compute quantile table ####
    ## Create empty table ####
    
  null_quantiles_tab_jacc <- as.data.frame(matrix(NA, nrow = 5, ncol = 41))
  rownames(null_quantiles_tab_jacc) <-  sort(unique(meta_mean$arms_name))
    
  null_quantiles_tab_turn <- as.data.frame(matrix(NA, nrow = 5, ncol = 41))
  rownames(null_quantiles_tab_turn) <-  sort(unique(meta_mean$arms_name))
  
  null_quantiles_tab_nest <- as.data.frame(matrix(NA, nrow = 5, ncol = 41))
  rownames(null_quantiles_tab_nest) <-  sort(unique(meta_mean$arms_name))
  
  null_quantiles_tab_bray <- as.data.frame(matrix(NA, nrow = 5, ncol = 41))
  rownames(null_quantiles_tab_bray) <-  sort(unique(meta_mean$arms_name))
    
    ## Make the loop ####

  
  for (prefix in unique_prefixes) {
    #prefix = "CINA1"
    subset_null <- subset(null_results, null_results$arms == prefix)
    
    
    #jacc # plutot faire des histo
    
    step = 0.025*(max(subset_null$jacc)-min(subset_null$jacc))
    brek = seq(min(subset_null$jacc),max(subset_null$jacc),  step)
    
    print(ggplot(subset_null, aes(x = jacc)) +
          geom_histogram(breaks = brek, fill = "coral", color = "black", alpha = 0.7) +
          labs(title = paste0("Frequency of Jaccard null values for ", prefix), x = "Beta-div index", y = "Frequency") +
          theme_minimal())
    
    quantile_jacc <- quantile(subset_null$jacc, probs = seq(0, 1, 0.025))
    
    null_quantiles_tab_jacc[prefix,] <- quantile_jacc
    
    #nest
    step = 0.025*(max(subset_null$nest)-min(subset_null$nest))
    brek = seq(min(subset_null$nest),max(subset_null$nest),  step)
    
    print(ggplot(subset_null, aes(x = nest)) +
          geom_histogram(breaks = brek, fill = "coral", color = "black", alpha = 0.7) +
          labs(title = paste0("Frequency of Nestedness null values for ", prefix), x = "Beta-div index", y = "Frequency") +
          theme_minimal())
    
    quantile_nest <- quantile(subset_null$nest, probs = seq(0, 1, 0.025))
    
    null_quantiles_tab_nest[prefix,] <- quantile_nest
    
    #turn
    step = 0.025*(max(subset_null$turn)-min(subset_null$turn))
    brek = seq(min(subset_null$turn),max(subset_null$turn),  step)
    
    print(ggplot(subset_null, aes(x = turn)) +
      geom_histogram(breaks = brek, fill = "coral", color = "black", alpha = 0.7) +
      labs(title = paste0("Frequency of Turnover null values for ", prefix), x = "Beta-div index", y = "Frequency") +
      theme_minimal())
    
    quantile_turn <- quantile(subset_null$turn, probs = seq(0, 1, 0.025))
    
    null_quantiles_tab_turn[prefix,] <- quantile_turn
    
    #bray
    step = 0.025*(max(subset_null$bray)-min(subset_null$bray))
    brek = seq(min(subset_null$bray),max(subset_null$bray),  step)
    
    print(ggplot(subset_null, aes(x = bray)) +
      geom_histogram(breaks = brek, fill = "coral", color = "black", alpha = 0.7) +
      labs(title = paste0("Frequency of Bray-Curtis null values for ", prefix), x = "Beta-div index", y = "Frequency") +
      theme_minimal())
    
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
  
  # Quantile checking function ####
  check_quantiles <- function(observed, quantiles) {
    if (observed <= quantiles["2.5%"]) {
      return("YES [0-2.5%]")
    } else if (observed >= quantiles["97.5%"]) {
      return("YES [97.5-100%]")
    } else {
      return("NO [2.5-97.5%]")
    }
  }
  
  
  # Create a new data frame to store the results
  tab_res <- data.frame(arms = results$arms, jacc = NA, turn = NA, nest = NA, bray = NA)
  
  # Create a table of results ####
  # Loop through each observed value and check against the corresponding quantile data frame
  for (i in 1:nrow(results)) {
    arms <- results$arms[i]
    tab_res$jacc[i] <- check_quantiles(results$jacc[i], null_quantiles$jacc[arms, ])
    tab_res$turn[i] <- check_quantiles(results$turn[i], null_quantiles$turn[arms, ])
    tab_res$nest[i] <- check_quantiles(results$nest[i], null_quantiles$turn[arms, ])
    tab_res$bray[i] <- check_quantiles(results$bray[i], null_quantiles$bray[arms, ])
  }
  
  # Print the results
  print(tab_res)
  
  
return(NULL)

}