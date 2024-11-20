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
  jaccard_results <- as.data.frame(matrix(NA, nrow = 96*15, ncol = 2))
  jaccard_results$names <-  rep(sort(unique(meta$arms_name)), each = 96)
  colnames(jaccard_results) <- c("same", "diff", "names")
  
  turnover_results <- as.data.frame(matrix(NA, nrow = 96*15, ncol = 2))
  turnover_results$names <-  rep(sort(unique(meta$arms_name)), each = 96)
  colnames(turnover_results) <- c("same", "diff", "names")
  
  nestedness_results <- as.data.frame(matrix(NA, nrow = 96*15, ncol = 2))
  nestedness_results$names <-  rep(sort(unique(meta$arms_name)), each = 96)
  colnames(turnover_results) <- c("same", "diff", "names")
  
  bray_results <- as.data.frame(matrix(NA, nrow = 96*15, ncol = 2))
  bray_results$names <-  rep(sort(unique(meta$arms_name)), each = 96)
  colnames(turnover_results) <- c("same", "diff", "names")
  
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
    
    jaccard_results$diff[jaccard_results$names == prefix] <- df.jacc$value
    
    
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
    
    turnover_results$diff[turnover_results$names == prefix] <- df.turn$value
    
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
    
    nestedness_results$diff[nestedness_results$names == prefix] <- df.nest$value
    
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
    
    bray_results$diff[bray_results$names == prefix] <- df.bray$value
    
  }
  
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
    df.jacc$row <- substr(df.jacc$row, 7, 8)
    df.jacc$col <- substr(df.jacc$col, 7, 8)
    df.jacc <- subset(df.jacc, row == col)
    
    jaccard_results$same[jaccard_results$names == prefix] <- df.jacc$value
    
    
    #Turnover
    
    df.turn <- melt(as.matrix(mat.turn), varnames = c("row", "col"))
    df.turn$row <- as.character(df.turn$row)
    df.turn$col <- as.character(df.turn$col)
    df.turn <- subset(df.turn, row != col)
    df.turn$row <- substr(df.turn$row, 7, 8)
    df.turn$col <- substr(df.turn$col, 7, 8)
    df.turn <- subset(df.turn, row == col)
    
    turnover_results$same[turnover_results$names == prefix] <- df.turn$value
    
    # Nestedness
    df.nest <- melt(as.matrix(mat.nest), varnames = c("row", "col"))
    df.nest$row <- as.character(df.nest$row)
    df.nest$col <- as.character(df.nest$col)
    df.nest <- subset(df.nest, row != col)
    df.nest$row <- substr(df.nest$row, 7, 8)
    df.nest$col <- substr(df.nest$col, 7, 8)
    df.nest <- subset(df.nest, row == col)
    
    nestedness_results$same[nestedness_results$names == prefix] <- df.nest$value
    
    # Bray
    df.bray <- melt(as.matrix(mat.bray), varnames = c("row", "col"))
    df.bray$row <- as.character(df.bray$row)
    df.bray$col <- as.character(df.bray$col)
    df.bray <- subset(df.bray, row != col)
    df.bray$row <- substr(df.bray$row, 7, 8)
    df.bray$col <- substr(df.bray$col, 7, 8)
    df.bray <- subset(df.bray, row == col)
    
    bray_results$same[bray_results$names == prefix] <- df.bray$value
    
  }
  
  
  results <- data.frame(arms = jaccard_results$names,
                        jacc_diff = jaccard_results$diff, 
                        turn_diff = turnover_results$diff, 
                        nest_diff = nestedness_results$diff,
                        bray_diff = bray_results$diff,
                        jacc_same = jaccard_results$same, 
                        turn_same = turnover_results$same, 
                        nest_same = nestedness_results$same,
                        bray_same = bray_results$same)
  
  data_long <- results %>%
    pivot_longer(cols = -arms, names_to = "metric_type", values_to = "value") %>%
    separate(metric_type, into = c("metric", "type"), sep = "_")
  
  # Create a new column combining arms and type
  data_long$arms_type <- paste(data_long$arms, data_long$type, sep = "_")
  
  
  custom_colors <- c(
    "CINA1A_diff" = "#CC66CC", "CINA1A_same" = "#CC66CC",
    "CINA1B_diff" = "#CC66CC", "CINA1B_same" = "#CC66CC",
    "CINA1C_diff" = "#CC66CC", "CINA1C_same" = "#CC66CC",
    "CINA3A_diff" = "#CC66CC", "CINA3A_same" = "#CC66CC",
    "CINA3B_diff" = "#CC66CC", "CINA3B_same" = "#CC66CC",
    "CINA3C_diff" = "#CC66CC", "CINA3C_same" = "#CC66CC",
    "CINA2A_diff" = "#1B9E77", "CINA2A_same" = "#1B9E77",
    "CINA2B_diff" = "#1B9E77", "CINA2B_same" = "#1B9E77",
    "CINA2C_diff" = "#1B9E77", "CINA2C_same" = "#1B9E77",
    "CINA4A_diff" = "#1B9E77", "CINA4A_same" = "#1B9E77",
    "CINA4B_diff" = "#1B9E77", "CINA4B_same" = "#1B9E77",
    "CINA4C_diff" = "#1B9E77", "CINA4C_same" = "#1B9E77",
    "RUNA2A_diff" = "#FF7F00", "RUNA2A_same" = "#FF7F00",
    "RUNA2B_diff" = "#FF7F00", "RUNA2B_same" = "#FF7F00",
    "RUNA2C_diff" = "#FF7F00", "RUNA2C_same" = "#FF7F00"
  )
  
  desired_order <- c("CINA1A", "CINA1B", "CINA1C", "CINA3A", "CINA3B", "CINA3C", 
                     "CINA2A", "CINA2B", "CINA2C", "CINA4A", "CINA4B", "CINA4C", 
                     "RUNA2A", "RUNA2B", "RUNA2C")
  
  data_long$arms <- factor(data_long$arms, levels = desired_order)
  data_long$metric <- factor(data_long$metric, levels = c("bray", "jacc", "turn", "nest"))
  
  kgds <- ggplot(data_long, aes(x = arms, y = value, fill = arms_type)) +
    geom_boxplot(position = position_dodge(width = 0.8)) +
    stat_summary(fun = mean, geom = "point", shape = 18, color = "red", size = 3, 
                 position = position_dodge(width = 0.8)) +  # Red diamond for mean
    scale_fill_manual(values = custom_colors) +  # Apply custom colors
    facet_wrap(~metric, scales = "free_y", ncol = 4) +  # 4 columns
    labs(title = "Boxplots of Metrics by Arms with Custom Colors and Means",
         x = "Arms", y = "Value") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")  # Remove the legend
  
  path_to_boxplot_UDOC <- paste0("outputs/beta/boxplot_betadiv_microhab_25_10_24.pdf")
  ggsave(filename =  path_to_boxplot_UDOC, plot = kgds, width = 18.5, height = 7.5)
  
  return(NULL) 
}
 
 