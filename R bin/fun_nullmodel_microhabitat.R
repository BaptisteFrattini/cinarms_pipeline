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
  
  unique_prefixes <- meta_red$set
  
  # Compute series of observed index ####
  ## Create empty tables ####
  jaccard_results <- as.data.frame(matrix(NA, nrow = 96*15, ncol = 1))
  jaccard_results$names <-  rep(sort(unique(meta$arms_name)), each = 96)
  
  turnover_results <- as.data.frame(matrix(NA, nrow = 96*15, ncol = 1))
  turnover_results$names <-  rep(sort(unique(meta$arms_name)), each = 96)
  
  nestedness_results <- as.data.frame(matrix(NA, nrow = 96*15, ncol = 1))
  nestedness_results$names <-  rep(sort(unique(meta_mean$arms_name)), each = 96)
  
  bray_results <- as.data.frame(matrix(NA, nrow = 96*15, ncol = 1))
  bray_results$names <-  rep(sort(unique(meta_mean$arms_name)), each = 96)
  
  unique_prefixes <- unique(meta_mean$arms_name)

  
  for (prefix in unique_prefixes) {
    #prefix = "CINA1A"
    subset_data <- subset(data, meta$arms_name == prefix)
    subset_meta_red <- subset(meta_red, meta_red$set == prefix)
    
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
    df.jacc$row <- substr(df.jacc$row, 6, 7)
    df.jacc$col <- substr(df.jacc$col, 6, 7)
    df.jacc <- subset(df.jacc, row != col)
    
    jaccard_results[[prefix]] <- df.jacc$value
   
    
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
    df.turn$row <- substr(df.turn$row, 6, 7)
    df.turn$col <- substr(df.turn$col, 6, 7)
    df.turn <- subset(df.turn, row != col)
    
    turnover_results[[prefix]] <- df.turn$value
    
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
    df.nest$row <- substr(df.nest$row, 6, 7)
    df.nest$col <- substr(df.nest$col, 6, 7)
    df.nest <- subset(df.nest, row != col)
    
    nestedness_results[[prefix]] <- df.nest$value
    
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
    df.bray$row <- substr(df.bray$row, 6, 7)
    df.bray$col <- substr(df.bray$col, 6, 7)
    df.bray <- subset(df.bray, row != col)
    
    bray_results[[prefix]] <- df.bray$value
    
  }
  
  
  results <- list(bray = bray_results,
                  jacc = jaccard_results,
                  turn = turnover_results,
                  nest = nestedness_results)
  
  
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
      df.jacc$row <- substr(df.jacc$row, 6, 7)
      df.jacc$col <- substr(df.jacc$col, 6, 7)
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
      df.turn$row <- substr(df.turn$row, 6, 7)
      df.turn$col <- substr(df.turn$col, 6, 7)
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
      df.nest$row <- substr(df.nest$row, 6, 7)
      df.nest$col <- substr(df.nest$col, 6, 7)
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
      df.bray$row <- substr(df.bray$row, 6, 7)
      df.bray$col <- substr(df.bray$col, 6, 7)
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
  
    
}
