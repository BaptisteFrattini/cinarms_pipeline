#' Mean all species cover of each plates of each arms 
#'
#' @param raw_data the path to the raw data file
#' @param arms_id the ID of the arms to subset for
#'
#' @return the path to the subseted derived data file
#' @export 

beta_microhab <- function(meta_data){
  
  # meta_data = targets::tar_read(metadata_data)
  library(ggpubr)
  library(forcats)
  library(reshape2) 
  library(ggplot2)
  library(dplyr)
  library(betapart)
  meta = read.csv(meta_data[grepl("metadata", meta_data)], header = TRUE)
  meta <- meta[,-4]
  data = read.csv(meta_data[!grepl("metadata", meta_data)], header = TRUE)
  data <- data[ , colSums(data) != 0]
  
  meta_red <- data.frame(prefixe = substr(meta$arms_name,1,5), arms = meta$arms_name, Orientation = meta$Orientation)
  
  df_mean_microhab <- data %>% 
    group_by(meta$Orientation,meta$arms_name) %>% 
    summarise_all(mean, na.rm = TRUE)
  
  
  
  # Initialize empty lists to store results
  
  bray_results <- list()
  jaccard_results <- list()
  turnover_results <- list()
  nestedness_results <- list()
  
  # Iterate through unique "prefixe" values
  unique_prefixes <- unique(meta_red$prefixe)
  
  for (prefix in unique_prefixes) {
    # Subset the data for the current prefix
    subset_data <-  data[meta_red$prefixe == prefix, ]
    subset_meta <-  meta_red[meta_red$prefixe == prefix, ]
    row.names(subset_data) <- paste0(rownames(subset_data), subset_meta$Orientation)
    mat.bray <- vegan::vegdist(subset_data, "bray")

    B.pair.pa <- betapart::beta.pair(vegan::decostand(subset_data, "pa"), index.family = "jaccard")
    
    mat.turn <- B.pair.pa$beta.jtu
    mat.nest <- B.pair.pa$beta.jne
    mat.jacc <- 1-B.pair.pa$beta.jac
    
    #Turnover
    df.turn <- melt(as.matrix(mat.turn), varnames = c("row", "col"))
    df.turn$row <- as.character(df.turn$row)
    df.turn$col <- as.character(df.turn$col)
    df.turn <- subset(df.turn, row != col)
    df.turn <- df.turn %>%
      mutate(
        row = as.character(row),  # Convert 'row' to character
        col = as.character(col),  # Convert 'col' to character
        combined = paste(pmin(row, col), pmax(row, col), sep = "")
      )
    df.turn <- df.turn[!duplicated(df.turn$combined), ]
    
    df.turn$row <- substr(df.turn$row, nchar(df.turn$row), nchar(df.turn$row))
    df.turn$col <- substr(df.turn$col, nchar(df.turn$col),nchar(df.turn$col))
    df.turn$same_value <- ifelse(df.turn$row == df.turn$col, "Yes", "No")
    
    #Nestedness
    df.nest <- melt(as.matrix(mat.nest), varnames = c("row", "col"))
    df.nest$row <- as.character(df.nest$row)
    df.nest$col <- as.character(df.nest$col)
    df.nest <- subset(df.nest, row != col)
    df.nest <- df.nest %>%
      mutate(
        row = as.character(row),  # Convert 'row' to character
        col = as.character(col),  # Convert 'col' to character
        combined = paste(pmin(row, col), pmax(row, col), sep = "")
      )
    df.nest <- df.nest[!duplicated(df.nest$combined), ]
    
    df.nest$row <- substr(df.nest$row, nchar(df.nest$row), nchar(df.nest$row))
    df.nest$col <- substr(df.nest$col, nchar(df.nest$col),nchar(df.nest$col))
    df.nest$same_value <- ifelse(df.nest$row == df.nest$col, "Yes", "No")
    
    #Jaccard
    df.jacc <- melt(as.matrix(mat.jacc), varnames = c("row", "col"))
    df.jacc$row <- as.character(df.jacc$row)
    df.jacc$col <- as.character(df.jacc$col)
    df.jacc <- subset(df.jacc, row != col)
    df.jacc <- df.jacc %>%
      mutate(
        row = as.character(row),  # Convert 'row' to character
        col = as.character(col),  # Convert 'col' to character
        combined = paste(pmin(row, col), pmax(row, col), sep = "")
      )
    df.jacc <- df.jacc[!duplicated(df.jacc$combined), ]
    
    df.jacc$row <- substr(df.jacc$row, nchar(df.jacc$row), nchar(df.jacc$row))
    df.jacc$col <- substr(df.jacc$col, nchar(df.jacc$col),nchar(df.jacc$col))
    df.jacc$same_value <- ifelse(df.jacc$row == df.jacc$col, "Yes", "No")
    
    # Bray Curtis
    df.bray <- melt(as.matrix(mat.bray), varnames = c("row", "col"))
    df.bray$row <- as.character(df.bray$row)
    df.bray$col <- as.character(df.bray$col)
    df.bray <- subset(df.bray, row != col)
    df.bray <- df.bray %>%
      mutate(
        row = as.character(row),  # Convert 'row' to character
        col = as.character(col),  # Convert 'col' to character
        combined = paste(pmin(row, col), pmax(row, col), sep = "")
      )
    df.bray <- df.bray[!duplicated(df.bray$combined), ]
    
    df.bray$row <- substr(df.bray$row, nchar(df.bray$row), nchar(df.bray$row))
    df.bray$col <- substr(df.bray$col, nchar(df.bray$col),nchar(df.bray$col))
    df.bray$same_value <- ifelse(df.bray$row == df.bray$col, "Yes", "No")
    
    # Store the results in the respective lists
    bray_results[[prefix]] <- df.bray
    jaccard_results[[prefix]] <- df.jacc
    turnover_results[[prefix]] <- df.turn
    nestedness_results[[prefix]] <-  df.nest
  }
  
  jaccard_results
  turnover_results
  nestedness_results
  bray_results
  
  df_jaccard <- data.frame(Same_value = jaccard_results$CINA1$same_value, 
                           CINA1 = jaccard_results$CINA1$value, 
                           CINA3 = jaccard_results$CINA3$value, 
                           CINA2 = jaccard_results$CINA2$value, 
                           CINA4 = jaccard_results$CINA4$value, 
                           RUNA2 = jaccard_results$RUNA2$value)
 
  df_jaccard_No <- subset(df_jaccard, Same_value == "No")
  df_jaccard_No <- melt(as.matrix(df_jaccard_No[,-1]))
  colnames(df_jaccard_No) <- c("num", "set", "value")
  
  df_turnover <- data.frame(Same_value = turnover_results$CINA1$same_value, 
                           CINA1 = turnover_results$CINA1$value, 
                           CINA3 = turnover_results$CINA3$value, 
                           CINA2 = turnover_results$CINA2$value, 
                           CINA4 = turnover_results$CINA4$value, 
                           RUNA2 = turnover_results$RUNA2$value)
  
  df_turnover_No <- subset(df_turnover, Same_value == "No")
  df_turnover_No <- melt(as.matrix(df_turnover_No[,-1]))
  colnames(df_turnover_No) <- c("num", "set", "value")
  
  df_nestedness <- data.frame(Same_value = nestedness_results$CINA1$same_value, 
                           CINA1 = nestedness_results$CINA1$value, 
                           CINA3 = nestedness_results$CINA3$value, 
                           CINA2 = nestedness_results$CINA2$value, 
                           CINA4 = nestedness_results$CINA4$value, 
                           RUNA2 = nestedness_results$RUNA2$value)
  
  df_nestedness_No <- subset(df_nestedness, Same_value == "No")
  df_nestedness_No <- melt(as.matrix(df_nestedness_No[,-1]))
  colnames(df_nestedness_No) <- c("num", "set", "value")
  
  
  df_bray <- data.frame(Same_value = bray_results$CINA1$same_value, 
                              CINA1 = bray_results$CINA1$value, 
                              CINA3 = bray_results$CINA3$value, 
                              CINA2 = bray_results$CINA2$value, 
                              CINA4 = bray_results$CINA4$value, 
                              RUNA2 = bray_results$RUNA2$value)
  
  df_bray_No <- subset(df_bray, Same_value == "No")
  df_bray_No <- melt(as.matrix(df_bray_No[,-1]))
  colnames(df_bray_No) <- c("num", "set", "value")
  
  #Jaccard
  intra = c("between upward \n and downward faces \n of the CINA1 set", 
            "between upward \n and downward faces \n of the CINA3 set", 
            "between upward \n and downward faces \n of the CINA2 set", 
            "between upward \n and downward faces \n of the CINA4 set", 
            "between upward \n and downward faces \n of the RUNA2 set")
  
  kk <- ggplot(df_jaccard_No, aes(x = fct_relevel(set, "CINA1", "CINA3", "CINA2", "CINA4", "RUNA2"), y = value)) +
    geom_boxplot(fill =  c("#CC66CC","#CC66CC","#1B9E77","#1B9E77","#FF7F00") ) +
    labs(title = "",
         x = "Comparisons",
         y = "Jaccard component") +
    theme(legend.position = "none") +
    scale_x_discrete(labels=intra) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.title.x = element_blank(), axis.title.y = element_text(size=12))
  
  gg <- ggplot(df_bray_No, aes(x = fct_relevel(set, "CINA1", "CINA3", "CINA2", "CINA4", "RUNA2"), y = value)) +
    geom_boxplot(fill =  c("#CC66CC","#CC66CC","#1B9E77","#1B9E77","#FF7F00") ) +
    labs(title = "",
         x = "Comparisons",
         y = "bray component") +
    theme(legend.position = "none") +
    scale_x_discrete(labels=intra) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.title.x = element_blank(), axis.title.y = element_text(size=12))
  
  dd <- ggplot(df_turnover_No, aes(x = fct_relevel(set, "CINA1", "CINA3", "CINA2", "CINA4", "RUNA2"), y = value)) +
    geom_boxplot(fill =  c("#CC66CC","#CC66CC","#1B9E77","#1B9E77","#FF7F00") ) +
    labs(title = "",
         x = "Comparisons",
         y = "Turnover component") +
    theme(legend.position = "none") +
    scale_x_discrete(labels=intra) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.title.x = element_blank(), axis.title.y = element_text(size=12))
  
  ss <- ggplot(df_nestedness_No, aes(x = fct_relevel(set, "CINA1", "CINA3", "CINA2", "CINA4", "RUNA2"), y = value)) +
    geom_boxplot(fill =  c("#CC66CC","#CC66CC","#1B9E77","#1B9E77","#FF7F00") ) +
    labs(title = "",
         x = "Comparisons",
         y = "Nestedness component") +
    theme(legend.position = "none") +
    scale_x_discrete(labels=intra) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.title.x = element_blank(), axis.title.y = element_text(size=12))
    
  fin <- cowplot::plot_grid(kk, gg, dd, ss,
                            ncol = 2,
                            nrow = 2)
  
  path_to_boxplot <- paste0("outputs/beta/boxplot_betadiv_microhab.pdf")
  ggsave(filename =  path_to_boxplot, plot = fin, width = 12, height = 14)
  
  
  
  return ( path_to_boxplot )
  
}