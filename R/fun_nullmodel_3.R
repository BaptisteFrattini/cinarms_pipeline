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
  
  # import data ####
  meta_mean <- read.csv(metadata_data_mean[grepl("metadata", metadata_data_mean)], header = TRUE)
  
  df_mean <- read.csv(metadata_data_mean[!grepl("metadata", metadata_data_mean)], header = TRUE)
  rownames(df_mean) <- meta_mean$arms
  
  
  df_mean_pa <- vegan::decostand(df_mean, "pa")
  rownames(df_mean_pa) <- meta_mean$arms
  
  dis.bray <- vegan::vegdist(df_mean, "bray")
  dis.jacc <- vegan::vegdist(df_mean_pa, "jaccard")
  
  
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
  # Create empty tables ####
  null_jaccard_results <- as.data.frame(matrix(NA, nrow = 3*5*1000, ncol = 1))
  null_jaccard_results$names <-  rep(sort(unique(meta_mean$arms_name)), each = 3*1000)
  null_jaccard_results$iterations <-  rep(c(1:1000), each = 3*5)

  null_turnover_results <- as.data.frame(matrix(NA, nrow = 3*5*1000, ncol = 1))
  null_turnover_results$names <-  rep(sort(unique(meta_mean$arms_name)), each = 3*1000)
  null_turnover_results$iterations <-  rep(c(1:1000), each = 3*5)

  null_nestedness_results <- as.data.frame(matrix(NA, nrow = 3*5*1000, ncol = 1))
  null_nestedness_results$names <-  rep(sort(unique(meta_mean$arms_name)), each = 3*1000)
  null_nestedness_results$iterations <-  rep(c(1:1000), each = 3*5)

  null_bray_results <- as.data.frame(matrix(NA, nrow = 3*5*1000, ncol = 1))
  null_bray_results$names <-  rep(sort(unique(meta_mean$arms_name)), each = 3*1000)
  null_bray_results$iterations <-  rep(c(1:1000), each = 3*5)

  unique_prefixes <- unique(meta_mean$arms_name)

    ## Make the loop ####

  # for(i in 1:1000) {
  #   #i = 1
  #     for (prefix in unique_prefixes) {
  #       #prefix = "CINA1"
  #       #presence absence
  #       null_model_pa <- sim9(df_mean_pa, nReps = 1, metric = "c_score", algo = "sim2")
  #       null_model_pa <- null_model_pa$Randomized.Data
  # 
  #       #abondance
  #       subset_data_ab <- round(df_mean*10000, digits = 0)
  #       null_model_ab <- canaper::cpr_rand_comm(subset_data_ab, "quasiswap_count", 100)
  #       null_model_ab <- null_model_ab/10000
  #       null_model_ab <- as.data.frame(vegan::decostand(null_model_ab, "total"))*100
  #       null_model_pa <- subset(null_model_pa, meta_mean$arms_name == prefix)
  #       # subset_data_pa <- subset_data_pa[, colSums(subset_data_pa) > 0]
  #       null_model_ab <- subset(null_model_ab, meta_mean$arms_name == prefix)
  #       # subset_data_ab <- subset_data_ab[, colSums(subset_data_ab) > 0]
  #       #Presence absence
  #       #Dissimilarity matrix computing
  #       B.pair.pa <- betapart::beta.pair(null_model_pa, index.family = "jaccard")
  #       #Decomposing beta-div
  #       mat.turn <- B.pair.pa$beta.jtu
  #       mat.nest <- B.pair.pa$beta.jne
  #       mat.jacc <- B.pair.pa$beta.jac
  #       #Jaccard
  #       df.jacc <- melt(as.matrix(mat.jacc), varnames = c("row", "col"))
  #       df.jacc$row <- as.character(df.jacc$row)
  #       df.jacc$col <- as.character(df.jacc$col)
  #       df.jacc <- subset(df.jacc, row != col)
  #       df.jacc <- df.jacc %>%
  #         mutate(
  #           row = as.character(row),  # Convert 'row' to character
  #           col = as.character(col),  # Convert 'col' to character
  #           combined = paste(pmin(row, col), pmax(row, col), sep = "_")
  #         )
  #       df.jacc <- df.jacc[!duplicated(df.jacc$combined), ]
  #       indices <- null_jaccard_results$names == prefix & null_jaccard_results$iterations == i
  #       null_jaccard_results$V1[indices] <- df.jacc$value
  # 
  #       #Turnover
  #       df.turn <- melt(as.matrix(mat.turn), varnames = c("row", "col"))
  #       df.turn$row <- as.character(df.turn$row)
  #       df.turn$col <- as.character(df.turn$col)
  #       df.turn <- subset(df.turn, row != col)
  #       df.turn <- df.turn %>%
  #         mutate(
  #           row = as.character(row),  # Convert 'row' to character
  #           col = as.character(col),  # Convert 'col' to character
  #           combined = paste(pmin(row, col), pmax(row, col), sep = "_")
  #         )
  #       df.turn <- df.turn[!duplicated(df.turn$combined), ]
  #       indices <- null_turnover_results$names == prefix & null_turnover_results$iterations == i
  #       null_turnover_results$V1[indices] <- df.turn$value
  # 
  #       #Nestedness
  #       df.nest <- melt(as.matrix(mat.nest), varnames = c("row", "col"))
  #       df.nest$row <- as.character(df.nest$row)
  #       df.nest$col <- as.character(df.nest$col)
  #       df.nest <- subset(df.nest, row != col)
  #       df.nest <- df.nest %>%
  #         mutate(
  #           row = as.character(row),  # Convert 'row' to character
  #           col = as.character(col),  # Convert 'col' to character
  #           combined = paste(pmin(row, col), pmax(row, col), sep = "_")
  #         )
  #       df.nest <- df.nest[!duplicated(df.nest$combined), ]
  #       indices <- null_nestedness_results$names == prefix & null_nestedness_results$iterations == i
  #       null_nestedness_results$V1[indices] <- df.nest$value
  # 
  # 
  #       #Bray Curtis
  #       mat.bray <- vegan::vegdist(null_model_ab, "bray")
  # 
  #       df.bray <- melt(as.matrix(mat.bray), varnames = c("row", "col"))
  #       df.bray$row <- as.character(df.bray$row)
  #       df.bray$col <- as.character(df.bray$col)
  #       df.bray <- subset(df.bray, row != col)
  #       df.bray <- df.bray %>%
  #         mutate(
  #           row = as.character(row),  # Convert 'row' to character
  #           col = as.character(col),  # Convert 'col' to character
  #           combined = paste(pmin(row, col), pmax(row, col), sep = "_")
  #         )
  #       df.bray <- df.bray[!duplicated(df.bray$combined), ]
  #       indices <- null_bray_results$names == prefix & null_bray_results$iterations == i
  #       null_bray_results$V1[indices] <- df.bray$value
  # 
  #   }
  # 
  # }

  null_results <- data.frame(arms = null_jaccard_results$names,
                             jacc = null_jaccard_results$V1,
                             turn = null_turnover_results$V1,
                             nest = null_nestedness_results$V1,
                             bray = null_bray_results$V1)

  # write.table(null_results, file = "outputs/null_model/null_results_1000_iter.csv", sep = ";", dec = ",")
  null_results <- read.table(file = "outputs/null_model/null_results_1000_iter.csv", sep = ";", dec = ",")
  
    ## Compute quantile table ####
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
  result_quantile <- data.frame(arms = results$arms,
                                jacc_quantile = rep(NA, 15),
                                turn_quantile = rep(NA, 15), 
                                nest_quantile = rep(NA, 15),
                                bray_quantile = rep(NA, 15))
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
  result_quantile$jacc <- factor(result_quantile$jacc, levels = quantile_levels, ordered = TRUE)
  result_quantile$turn <- factor(result_quantile$turn, levels = quantile_levels, ordered = TRUE)
  result_quantile$nest <- factor(result_quantile$nest, levels = quantile_levels, ordered = TRUE)
  result_quantile$bray <- factor(result_quantile$bray, levels = quantile_levels, ordered = TRUE)
  # Convert the factor to numeric to use in the plot
  result_quantile$jacc_numeric <- as.numeric(result_quantile$jacc)
  result_quantile$turn_numeric <- as.numeric(result_quantile$turn)
  result_quantile$nest_numeric <- as.numeric(result_quantile$nest)
  result_quantile$bray_numeric <- as.numeric(result_quantile$bray)
  
  # Create the data frame for plotting nest quantiles
  plot_data <- data.frame(
    Row = 1:nrow(result_quantile),  # X-axis: Rows (Comparisons)
    Quantile_jacc = result_quantile$jacc_numeric,
    Quantile_turn = result_quantile$turn_numeric,
    Quantile_nest = result_quantile$nest_numeric,
    Quantile_bray = result_quantile$bray_numeric# Y-axis: Numeric quantile values for 'nest'
  )
  
  plot_data$arms <- result_quantile$arms
  
  write.table(result_quantile, file = "outputs/null_model/Quantiles_1000_iter.csv", sep = ";", dec = ",")
  
  
    ## Create quantile plot (nest) ####
  intra = c("between ARMS of \n the CINA1 batch", "between ARMS of \n the CINA3 batch","between ARMS of \n the CINA2 batch","between ARMS of \n the CINA4 batch","between ARMS of \n the RUNA2 batch")
  
  dd1 <- ggplot(plot_data, 
                aes(x = fct_relevel(arms, "CINA1", "CINA3", "CINA2", "CINA4", "RUNA2"), 
                    y = Quantile_jacc)) +
    geom_jitter(aes(color = arms), width = 0.2, size = 3.5, show.legend = FALSE) +
    labs(title = "Jaccard dissimilarity",
         x = "Comparisons",
         y = "", size = 16) +
    scale_x_discrete(labels = intra) +
    scale_color_manual(values = c("CINA1" = "#CC66CC", "CINA3" = "#CC66CC", 
                                  "CINA2" = "#1B9E77", "CINA4" = "#1B9E77", 
                                  "RUNA2" = "#FF7F00")) +
    scale_y_continuous(breaks = 1:41, labels = quantile_levels) +  # 41 quantiles sur l'axe Y avec labels
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA), # Changer l'arrière-plan
      panel.grid.major = element_line(color = "grey"),  # Quadrillage principal en blanc
      panel.grid.minor = element_line(color = "white", linetype = "dashed"),  # Quadrillage mineur
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14), 
      axis.title.x = element_blank(), 
      axis.text.y = element_text(size = 8.5),
      plot.title = element_text(size = 18, face = "bold")
    ) +
    geom_hline(yintercept = 2, colour = "red") +
    geom_hline(yintercept = 40, colour = "red") +
    geom_hline(yintercept = 20, colour = "darkgrey")
  

  
  
  dd2 <- ggplot(plot_data, 
                aes(x = fct_relevel(arms, "CINA1", "CINA3", "CINA2", "CINA4", "RUNA2"), 
                    y = Quantile_turn)) +
    geom_jitter(aes(color = arms), width = 0.2, size = 3.5, show.legend = FALSE) +
    labs(title = "Turnover component",
         x = "Comparisons",
         y = "", size = 16) +
    scale_x_discrete(labels = intra) +
    scale_color_manual(values = c("CINA1" = "#CC66CC", "CINA3" = "#CC66CC", 
                                  "CINA2" = "#1B9E77", "CINA4" = "#1B9E77", 
                                  "RUNA2" = "#FF7F00")) +
    scale_y_continuous(breaks = 1:41, labels = quantile_levels) +  # 41 quantiles sur l'axe Y avec labels
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA), # Changer l'arrière-plan
      panel.grid.major = element_line(color = "grey"),  # Quadrillage principal en blanc
      panel.grid.minor = element_line(color = "white", linetype = "dashed"),  # Quadrillage mineur
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14), 
      axis.title.x = element_blank(), 
      axis.text.y = element_text(size = 8.5),
      plot.title = element_text(size = 18, face = "bold")
    ) +
    geom_hline(yintercept = 2, colour = "red") +
    geom_hline(yintercept = 40, colour = "red") +
    geom_hline(yintercept = 20, colour = "darkgrey")
  
  
  dd3 <- ggplot(plot_data, 
                aes(x = fct_relevel(arms, "CINA1", "CINA3", "CINA2", "CINA4", "RUNA2"), 
                    y = Quantile_nest)) +
    geom_jitter(aes(color = arms), width = 0.2, size = 3.5, show.legend = FALSE) +
    labs(title = "Nestedness component",
         x = "Comparisons",
         y = "", size = 16) +
    scale_x_discrete(labels = intra) +
    scale_color_manual(values = c("CINA1" = "#CC66CC", "CINA3" = "#CC66CC", 
                                  "CINA2" = "#1B9E77", "CINA4" = "#1B9E77", 
                                  "RUNA2" = "#FF7F00")) +
    scale_y_continuous(breaks = 1:41, labels = quantile_levels) +  # 41 quantiles sur l'axe Y avec labels
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA), # Changer l'arrière-plan
      panel.grid.major = element_line(color = "grey"),  # Quadrillage principal en blanc
      panel.grid.minor = element_line(color = "white", linetype = "dashed"),  # Quadrillage mineur
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14), 
      axis.title.x = element_blank(), 
      axis.text.y = element_text(size = 8.5),
      plot.title = element_text(size = 18, face = "bold")
    ) +
    geom_hline(yintercept = 2, colour = "red") +
    geom_hline(yintercept = 40, colour = "red") +
    geom_hline(yintercept = 20, colour = "darkgrey")
  
  
  dd4 <- ggplot(plot_data, 
                aes(x = fct_relevel(arms, "CINA1", "CINA3", "CINA2", "CINA4", "RUNA2"), 
                    y = Quantile_bray)) +
    geom_jitter(aes(color = arms), width = 0.2, size = 3.5, show.legend = FALSE) +
    labs(title = "Bray-Curtis dissimilarity",
         x = "Comparisons",
         y = "Percentiles", size = 16) +
    scale_x_discrete(labels = intra) +
    scale_color_manual(values = c("CINA1" = "#CC66CC", "CINA3" = "#CC66CC", 
                                  "CINA2" = "#1B9E77", "CINA4" = "#1B9E77", 
                                  "RUNA2" = "#FF7F00")) +
    scale_y_continuous(breaks = 1:41, labels = quantile_levels) +  # 41 quantiles sur l'axe Y avec labels
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA), # Changer l'arrière-plan
      panel.grid.major = element_line(color = "grey"),  # Quadrillage principal en blanc
      panel.grid.minor = element_line(color = "white", linetype = "dashed"),  # Quadrillage mineur
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14), 
      axis.title.x = element_blank(), 
      axis.text.y = element_text(size = 8.5),
      plot.title = element_text(size = 18, face = "bold")
    ) +
    geom_hline(yintercept = 2, colour = "red") +
    geom_hline(yintercept = 40, colour = "red") +
    geom_hline(yintercept = 20, colour = "darkgrey")
  
  
  fin <- cowplot::plot_grid(dd4, dd1, dd2, dd3, 
                            ncol = 4,
                            nrow = 1)
  
  path_to_boxplot <- paste0("outputs/null_model/boxplot_null_model_16_05_2025.pdf")
  ggsave(filename =  path_to_boxplot, plot = fin, width = 21.5, height = 7.5)
  
  

  
  
  # ## SES ####
  # 
  # null_SES_tab_jacc <- as.data.frame(matrix(NA, nrow = 15, ncol = 1))
  # null_SES_tab_jacc$names <-  rep(sort(unique(meta_mean$arms_name)), each = 3)
  # 
  # null_SES_tab_turn <- as.data.frame(matrix(NA, nrow = 15, ncol = 1))
  # null_SES_tab_turn$names <-  rep(sort(unique(meta_mean$arms_name)), each = 3)
  # 
  # null_SES_tab_nest <- as.data.frame(matrix(NA, nrow = 15, ncol = 1))
  # null_SES_tab_nest$names <-  rep(sort(unique(meta_mean$arms_name)), each = 3)
  # 
  # null_SES_tab_bray <- as.data.frame(matrix(NA, nrow = 15, ncol = 1))
  # null_SES_tab_bray$names <-  rep(sort(unique(meta_mean$arms_name)), each = 3)
  # 
  # 
  # for (prefix in unique_prefixes) {
  #   #prefix = "CINA1"
  #   subset_null <- subset(null_results, null_results$arms == prefix)
  #   subset_results <- subset(results, results$arms == prefix)
  # 
  #   SES_jacc <- (subset_results$jacc - mean(subset_null$jacc))/sd(subset_null$jacc)
  # 
  #   null_SES_tab_jacc$V1[null_SES_tab_jacc$names == prefix] <- SES_jacc
  #   
  #   SES_nest <- (subset_results$nest - mean(subset_null$nest))/sd(subset_null$nest)
  #   
  #   null_SES_tab_nest$V1[null_SES_tab_nest$names == prefix] <- SES_nest
  #   
  #   SES_turn <- (subset_results$turn - mean(subset_null$turn))/sd(subset_null$turn)
  #   
  #   null_SES_tab_turn$V1[null_SES_tab_turn$names == prefix] <- SES_turn
  #   
  #   SES_bray <- (subset_results$bray - mean(subset_null$bray))/sd(subset_null$bray)
  #   
  #   null_SES_tab_bray$V1[null_SES_tab_bray$names == prefix] <- SES_bray
  # 
  # }
  # 
  # SES_results <- data.frame(arms = null_SES_tab_jacc$names,
  #                           jacc = null_SES_tab_jacc$V1, 
  #                           turn = null_SES_tab_turn$V1, 
  #                           nest = null_SES_tab_nest$V1,
  #                           bray = null_SES_tab_bray$V1)
  # 
  # write.table(SES_results, file = "outputs/null_model/SES_1000_iter.csv", sep = ";", dec = ",")
  # 
  #   ## Plot SES results ####
  # intra = c("between ARMS of \n the CINA1 batch", "between ARMS of \n the CINA3 batch","between ARMS of \n the CINA2 batch","between ARMS of \n the CINA4 batch","between ARMS of \n the RUNA2 batch")
  # 
  # kk <- ggplot(SES_results, 
  #              aes(x = fct_relevel(arms, "CINA1", "CINA3", "CINA2", "CINA4", "RUNA2"), 
  #                  y = jacc)) +
  #   geom_jitter(aes(color = arms), width = 0.2, size = 2, show.legend = FALSE) +
  #   labs(title = "Jaccard dissimilarity",
  #        x = "Comparisons",
  #        y = "SES") +
  #   theme(legend.position = "none") +
  #   scale_x_discrete(labels = intra) +
  #   scale_color_manual(values = c("CINA1" = "#CC66CC", "CINA3" = "#CC66CC", 
  #                                 "CINA2" = "#1B9E77", "CINA4" = "#1B9E77", 
  #                                 "RUNA2" = "#FF7F00")) +
  #   theme_classic() +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), 
  #         axis.title.x = element_blank(), 
  #         axis.text.y = element_text(size = 12)) +
  #   geom_hline(yintercept = -1.96, colour = "red") +
  #   geom_hline(yintercept = 1.96, colour = "red") +
  #   geom_hline(yintercept = 0, colour = "darkgrey") +
  #   ylim(min = -4.5, max = 4.5) 
  # 
  # gg <- ggplot(SES_results, 
  #              aes(x = fct_relevel(arms, "CINA1", "CINA3", "CINA2", "CINA4", "RUNA2"), 
  #                  y = turn)) +
  #   geom_jitter(aes(color = arms), width = 0.2, size = 2, show.legend = FALSE) +
  #   labs(title = "Turnover component",
  #        x = "Comparisons",
  #        y = "SES") +
  #   theme(legend.position = "none") +
  #   scale_x_discrete(labels = intra) +
  #   scale_color_manual(values = c("CINA1" = "#CC66CC", "CINA3" = "#CC66CC", 
  #                                 "CINA2" = "#1B9E77", "CINA4" = "#1B9E77", 
  #                                 "RUNA2" = "#FF7F00")) +
  #   theme_classic() +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), 
  #         axis.title.x = element_blank(), 
  #         axis.text.y = element_text(size = 12)) +
  #   geom_hline(yintercept = -1.96, colour = "red") +
  #   geom_hline(yintercept = 1.96, colour = "red") +
  #   geom_hline(yintercept = 0, colour = "darkgrey") +
  #   ylim(min = -4.5, max = 4.5) 
  # 
  # # dd <- ggplot(SES_results, 
  # #              aes(x = fct_relevel(arms, "CINA1", "CINA3", "CINA2", "CINA4", "RUNA2"), 
  # #                  y = nest)) +
  # #   geom_jitter(aes(color = arms), width = 0.2, size = 2, show.legend = FALSE) +
  # #   labs(title = "Nestedness component",
  # #        x = "Comparisons",
  # #        y = "") +
  # #   theme(legend.position = "none") +
  # #   scale_x_discrete(labels = intra) +
  # #   scale_color_manual(values = c("CINA1" = "#CC66CC", "CINA3" = "#CC66CC", 
  # #                                 "CINA2" = "#1B9E77", "CINA4" = "#1B9E77", 
  # #                                 "RUNA2" = "#FF7F00")) +
  # #   theme_classic() +
  # #   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), 
  # #         axis.title.x = element_blank(), 
  # #         axis.title.y = element_blank(),
  # #         axis.text.y = element_text(size = 12)) +
  # #   geom_hline(yintercept = -1.96, colour = "red") +
  # #   geom_hline(yintercept = 1.96, colour = "red") +
  # #   geom_hline(yintercept = 0, colour = "darkgrey") +
  # #   ylim(min = -4.5, max = 4.5)
  # 
  # ss <- ggplot(SES_results, 
  #              aes(x = fct_relevel(arms, "CINA1", "CINA3", "CINA2", "CINA4", "RUNA2"), 
  #                  y = bray)) +
  #   geom_jitter(aes(color = arms), width = 0.2, size = 2, show.legend = FALSE) +
  #   labs(title = "Bray-Curtis dissimilarity",
  #        x = "Comparisons",
  #        y = "SES") +
  #   theme(legend.position = "none") +
  #   scale_x_discrete(labels = intra) +
  #   scale_color_manual(values = c("CINA1" = "#CC66CC", "CINA3" = "#CC66CC", 
  #                                 "CINA2" = "#1B9E77", "CINA4" = "#1B9E77", 
  #                                 "RUNA2" = "#FF7F00")) +
  #   theme_classic() +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), 
  #         axis.title.x = element_blank(), 
  #         axis.text.y = element_text(size = 12)) +
  #   geom_hline(yintercept = -1.96, colour = "red") +
  #   geom_hline(yintercept = 1.96, colour = "red") +
  #   geom_hline(yintercept = 0, colour = "darkgrey") +
  #   ylim(min = -4.5, max = 4.5)  
  # 
  # 
  # fin <- cowplot::plot_grid(ss, kk, gg, dd,
  #                           ncol = 4,
  #                           nrow = 1)
  # 
  # path_to_boxplot <- paste0("outputs/null_model/boxplot_null_model_full_22_10_24.pdf")
  # ggsave(filename =  path_to_boxplot, plot = fin, width = 18.5, height = 7.5)
  # 
  # return ####
  
  return(NULL)
  
}
