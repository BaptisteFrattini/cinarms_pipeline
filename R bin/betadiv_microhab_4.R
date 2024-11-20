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
  library(tidyr)
  library(EcoSimR)
  library(tibble)
  meta = read.csv(meta_data[grepl("metadata", meta_data)], header = TRUE)
  meta <- meta[,-4]
  data = read.csv(meta_data[!grepl("metadata", meta_data)], header = TRUE)
  data <- data[ , colSums(data) != 0]
  
  meta_red <- data.frame(prefixe = substr(meta$arms_name,1,5), 
                         set = meta$arms_name, 
                         Orientation = meta$Orientation,
                         Open_Close = meta$o_c)
  data

  #### Comparaison UDOC sans pooler ####
  
  unique_prefixes <- meta_red$set

  turnover_results <- as.data.frame(matrix(NA, nrow = 6*15, ncol = 1))
  turnover_results$names <- rep(sort(unique(meta$arms_name)), each = 6)
  
  nestedness_results <- as.data.frame(matrix(NA, nrow = 6*15, ncol = 1))
  nestedness_results$names <-  rep(sort(unique(meta$arms_name)), each = 6)
  
  jaccard_results <- as.data.frame(matrix(NA, nrow = 6*15, ncol = 1))
  jaccard_results$names <-  rep(sort(unique(meta$arms_name)), each = 6)
  
  bray_results <- as.data.frame(matrix(NA, nrow = 6*15, ncol = 1))
  bray_results$names <-  rep(sort(unique(meta$arms_name)), each = 6)
  
  
  unique_prefixes <- unique(meta_red$set)


  for (prefix in unique_prefixes) {
    #prefix = "CINA1A"

    rownames(data) <- paste0(meta$arms_name, "_", meta$Orientation, meta$o_c, meta$plate_number)
    subset_data <- subset(data, meta$arms_name == prefix)
    
    subset_meta_red <- subset(meta_red, meta_red$prefixe == prefix)
   
    mat.bray <- vegan::vegdist(subset_data, "bray")

    B.pair.pa <- betapart::beta.pair(vegan::decostand(subset_data, "pa"), index.family = "jaccard")

    mat.turn <- B.pair.pa$beta.jtu
    mat.nest <- B.pair.pa$beta.jne
    mat.jacc <- B.pair.pa$beta.jac

    #Turnover

    df.turn <- melt(as.matrix(mat.turn), varnames = c("row", "col"))
 
    df.turn$row <- substr(df.turn$row, 1, nchar(as.character(df.turn$row)) - 1)
    df.turn$col <- substr(df.turn$col, 1, nchar(as.character(df.turn$col)) - 1)

    
    df.turn <- subset(df.turn, row != col)
    df.turn <- df.turn %>%
      mutate(
        row = as.character(row),  # Convert 'row' to character
        col = as.character(col),  # Convert 'col' to character
        combined = paste(pmin(row, col), pmax(row, col), sep = "_")
      )
    df.turn <- df.turn[!duplicated(df.turn$combined), ]
    df.turn <- df.turn[,-4]

    turnover_results$V1[turnover_results$names == prefix] <- df.turn$value

    # Nestedness
    df.nest <- melt(as.matrix(mat.nest), varnames = c("row", "col"))
    
    df.nest$row <- substr(df.nest$row, 1, nchar(as.character(df.nest$row)) - 1)
    df.nest$col <- substr(df.nest$col, 1, nchar(as.character(df.nest$col)) - 1)
    
    df.nest <- subset(df.nest, row != col)
    df.nest <- df.nest %>%
      mutate(
        row = as.character(row),  # Convert 'row' to character
        col = as.character(col),  # Convert 'col' to character
        combined = paste(pmin(row, col), pmax(row, col), sep = "_")
      )
    df.nest <- df.nest[!duplicated(df.nest$combined), ]
    df.nest <- df.nest[,-4]

    nestedness_results$V1[nestedness_results$names == prefix] <- df.nest$value

    # Jaccard
    df.jacc <- melt(as.matrix(mat.jacc), varnames = c("row", "col"))
    
    df.jacc$row <- substr(df.jacc$row, 1, nchar(as.character(df.jacc$row)) - 1)
    df.jacc$col <- substr(df.jacc$col, 1, nchar(as.character(df.jacc$col)) - 1)
    
    df.jacc <- subset(df.jacc, row != col)
    df.jacc <- df.jacc %>%
      mutate(
        row = as.character(row),  # Convert 'row' to character
        col = as.character(col),  # Convert 'col' to character
        combined = paste(pmin(row, col), pmax(row, col), sep = "_")
      )
    df.jacc <- df.jacc[!duplicated(df.jacc$combined), ]
    df.jacc <- df.jacc[,-4]

    jaccard_results$V1[jaccard_results$names == prefix] <- df.jacc$value

    # Bray
    df.bray <- melt(as.matrix(mat.bray), varnames = c("row", "col"))
    
    df.bray$row <- substr(df.bray$row, 1, nchar(as.character(df.bray$row)) - 1)
    df.bray$col <- substr(df.bray$col, 1, nchar(as.character(df.bray$col)) - 1)
    
    df.bray <- subset(df.bray, row != col)
    df.bray <- df.bray %>%
      mutate(
        row = as.character(row),  # Convert 'row' to character
        col = as.character(col),  # Convert 'col' to character
        combined = paste(pmin(row, col), pmax(row, col), sep = "_")
      )
    df.bray <- df.bray[!duplicated(df.bray$combined), ]
    df.bray <- df.bray[,-4]

    bray_results$V1[bray_results$names == prefix] <- df.bray$value



  }


  results <- data.frame(bray = bray_results$V1,
                        jacc = jaccard_results$V1,
                        turn = turnover_results$V1,
                        nest = nestedness_results$V1)

  rownames(results) <- rownames(bray_results)
  
  #### compute null data ####
  unique_prefixes <- meta_red$set
  
  turnover_results <- as.data.frame(matrix(NA, nrow = 6*15, ncol = 1))
  turnover_results$names <- rep(sort(unique(meta$arms_name)), each = 6)
  
  nestedness_results <- as.data.frame(matrix(NA, nrow = 6*15, ncol = 1))
  nestedness_results$names <-  rep(sort(unique(meta$arms_name)), each = 6)
  
  jaccard_results <- as.data.frame(matrix(NA, nrow = 6*15, ncol = 1))
  jaccard_results$names <-  rep(sort(unique(meta$arms_name)), each = 6)
  
  bray_results <- as.data.frame(matrix(NA, nrow = 6*15, ncol = 1))
  bray_results$names <-  rep(sort(unique(meta$arms_name)), each = 6)
  
  
  unique_prefixes <- unique(meta_red$set)
  
  
  num_iterations <- 1000  
  data_t_pa <- vegan::decostand(t(data), "pa")

  null_model <- EcoSimR::sim9(data_t_pa, nReps = num_iterations, metric = "c_score", algo = "sim2")
  
  summary(null_model)
  null_model_data <- t(null_model$Randomized.Data)
  
  plot(null_model,type="burn_in")
  plot(null_model,type="hist")
  plot(null_model,type="cooc")
  
  data_rounded <- round(data*10000, digits = 0)
  data_null_bray <- canaper::cpr_rand_comm(data_rounded, "quasiswap_count", 10000)
  data_null_bray <- data_null_bray/10000
  
  
  data_null_bray <- as.data.frame(vegan::decostand(data_null_bray, "total"))*100
  
  for (prefix in unique_prefixes) {
    #prefix = "CINA1A"
    
    rownames(null_model_data) <- paste0(meta$arms_name, "_", meta$Orientation, meta$o_c, meta$plate_number)
    subset_data <- subset(null_model_data, meta$arms_name == prefix)
    
    subset_meta_red <- subset(meta_red, meta_red$prefixe == prefix)
    
    mat.bray <- vegan::vegdist(data_null_bray, "bray")
    
    B.pair.pa <- betapart::beta.pair(vegan::decostand(null_model_data, "pa"), index.family = "jaccard")
    
    mat.turn <- B.pair.pa$beta.jtu
    mat.nest <- B.pair.pa$beta.jne
    mat.jacc <- B.pair.pa$beta.jac
    
    #Turnover
    
    df.turn <- melt(as.matrix(mat.turn), varnames = c("row", "col"))
    
    df.turn$row <- substr(df.turn$row, 1, nchar(as.character(df.turn$row)) - 1)
    df.turn$col <- substr(df.turn$col, 1, nchar(as.character(df.turn$col)) - 1)
    
    
    df.turn <- subset(df.turn, row != col)
    df.turn <- df.turn %>%
      mutate(
        row = as.character(row),  # Convert 'row' to character
        col = as.character(col),  # Convert 'col' to character
        combined = paste(pmin(row, col), pmax(row, col), sep = "_")
      )
    df.turn <- df.turn[!duplicated(df.turn$combined), ]
    df.turn <- df.turn[,-4]
    
    turnover_results$V1[turnover_results$names == prefix] <- df.turn$value
    
    # Nestedness
    df.nest <- melt(as.matrix(mat.nest), varnames = c("row", "col"))
    
    df.nest$row <- substr(df.nest$row, 1, nchar(as.character(df.nest$row)) - 1)
    df.nest$col <- substr(df.nest$col, 1, nchar(as.character(df.nest$col)) - 1)
    
    df.nest <- subset(df.nest, row != col)
    df.nest <- df.nest %>%
      mutate(
        row = as.character(row),  # Convert 'row' to character
        col = as.character(col),  # Convert 'col' to character
        combined = paste(pmin(row, col), pmax(row, col), sep = "_")
      )
    df.nest <- df.nest[!duplicated(df.nest$combined), ]
    df.nest <- df.nest[,-4]
    
    nestedness_results$V1[nestedness_results$names == prefix] <- df.nest$value
    
    # Jaccard
    df.jacc <- melt(as.matrix(mat.jacc), varnames = c("row", "col"))
    
    df.jacc$row <- substr(df.jacc$row, 1, nchar(as.character(df.jacc$row)) - 1)
    df.jacc$col <- substr(df.jacc$col, 1, nchar(as.character(df.jacc$col)) - 1)
    
    df.jacc <- subset(df.jacc, row != col)
    df.jacc <- df.jacc %>%
      mutate(
        row = as.character(row),  # Convert 'row' to character
        col = as.character(col),  # Convert 'col' to character
        combined = paste(pmin(row, col), pmax(row, col), sep = "_")
      )
    df.jacc <- df.jacc[!duplicated(df.jacc$combined), ]
    df.jacc <- df.jacc[,-4]
    
    jaccard_results$V1[jaccard_results$names == prefix] <- df.jacc$value
    
    # Bray
    df.bray <- melt(as.matrix(data_null_bray), varnames = c("row", "col"))
    
    df.bray$row <- substr(df.bray$row, 1, nchar(as.character(df.bray$row)) - 1)
    df.bray$col <- substr(df.bray$col, 1, nchar(as.character(df.bray$col)) - 1)
    
    df.bray <- subset(df.bray, row != col)
    df.bray <- df.bray %>%
      mutate(
        row = as.character(row),  # Convert 'row' to character
        col = as.character(col),  # Convert 'col' to character
        combined = paste(pmin(row, col), pmax(row, col), sep = "_")
      )
    df.bray <- df.bray[!duplicated(df.bray$combined), ]
    df.bray <- df.bray[,-4]
    
    bray_results$V1[bray_results$names == prefix] <- df.bray$value
    
    
    
  }
  
  null_results <- data.frame(jacc = jaccard_results$V1, 
                             turn = turnover_results$V1, 
                             nest = nestedness_results$V1,
                             bray = bray_results$V1)
  
  rownames(null_results) <- rownames(jaccard_results)
  
  
  ## Compute the null deviation data ##
  
  ## JACCARD
  
  null.dev.jacc <- (results$jacc - mean(null_results$jacc))/sd(null_results$jacc)
  
  ## NESTEDNESS
  
  null.dev.nest <- (results$nest - mean(null_results$nest))/sd(null_results$nest)
  
  ## TURNOVER
  
  null.dev.turn <- (results$turn - mean(null_results$turn))/sd(null_results$turn)
  
  ## BRAY
  
  null.dev.bray <- (results$bray - mean(null_results$bray))/sd(null_results$bray)
  
  SES_data <- data.frame(ses_jacc = null.dev.jacc, ses_turn = null.dev.turn, ses_nest = null.dev.nest, ses_bray = null.dev.bray)

  #Jaccard
  intra = c("between different \n microhabitats of \n the CINA1 batch", 
            "between different \n microhabitats of \n the CINA3 batch", 
            "between different \n microhabitats of \n the CINA2 batch", 
            "between different \n microhabitats of \n the CINA4 batch", 
            "between different \n microhabitats of \n the RUNA2 batch")
  
  set_name <- sort(rep(unique(substr(meta$arms_name, 1, 5)),3*6))
  
  SES_data$comp <- set_name
  SES_data$imm_time <- rep(c("6 months","6 months","6 months", 
                         "1 year", "1 year", "1 year", 
                         "6 months","6 months","6 months",
                         "1 year", "1 year", "1 year", 
                         "2 years","2 years","2 years"), each = 6)
  
  #Jaccard
  
  means <- aggregate(ses_jacc ~  set_name, SES_data, mean)
  
  ll <- ggplot(SES_data, aes(x = fct_relevel(set_name, "CINA1", "CINA3", "CINA2", "CINA4", "RUNA2"), y = ses_jacc)) +
    geom_boxplot(color = c("#CC66CC", "#CC66CC", "#1B9E77", "#1B9E77", "#FF7F00"), fill = NA) +
    geom_jitter(aes(color = set_name), width = 0.2, size = 1.5, alpha = 0.6, show.legend = FALSE) +  # Points colorés selon les catégories
    scale_color_manual(values = c("CINA1" = "#CC66CC", "CINA3" = "#CC66CC", "CINA2" = "#1B9E77", "CINA4" = "#1B9E77", "RUNA2" = "#FF7F00")) +
    labs(title = "Jaccard dissimilarity",
         x = "Comparisons",
         y = "") +
    theme(legend.position = "none") +
    scale_x_discrete(labels = intra) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), 
          axis.title.x = element_blank(), 
          axis.title.y = element_blank(),
          axis.text.y = element_text(size = 12)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11), 
          axis.title.x = element_blank(), 
          axis.title.y = element_text(size = 12)) +
    geom_hline(yintercept = -1.96, colour = "red") +
    geom_hline(yintercept = 1.96, colour = "red") +
    geom_hline(yintercept = 0, colour = "darkgrey") +
    ylim(min = -4.5, max = 4.5) +
    stat_summary(fun = mean, colour = "darkred", geom = "point", 
                 shape = 18, size = 3, show.legend = FALSE)

  # Turnover
  means <- aggregate(ses_turn ~  set_name, SES_data, mean)
  
  mm <- ggplot(SES_data, aes(x = fct_relevel(set_name, "CINA1", "CINA3", "CINA2", "CINA4", "RUNA2"), y = ses_turn)) +
    geom_boxplot(color = c("#CC66CC", "#CC66CC", "#1B9E77", "#1B9E77", "#FF7F00"), fill = NA) +
    geom_jitter(aes(color = set_name), width = 0.2, size = 1.5, alpha = 0.6, show.legend = FALSE) +  # Points colorés selon les catégories
    scale_color_manual(values = c("CINA1" = "#CC66CC", "CINA3" = "#CC66CC", "CINA2" = "#1B9E77", "CINA4" = "#1B9E77", "RUNA2" = "#FF7F00")) +
    labs(title = "Turnover",
         x = "Comparisons",
         y = "") +
    theme(legend.position = "none") +
    scale_x_discrete(labels = intra) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), 
          axis.title.x = element_blank(), 
          axis.title.y = element_blank(),
          axis.text.y = element_text(size = 12)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11), 
          axis.title.x = element_blank(), 
          axis.title.y = element_text(size = 12)) +
    geom_hline(yintercept = -1.96, colour = "red") +
    geom_hline(yintercept = 1.96, colour = "red") +
    geom_hline(yintercept = 0, colour = "darkgrey") +
    ylim(min = -4.5, max = 4.5) +
    stat_summary(fun = mean, colour = "darkred", geom = "point", 
                 shape = 18, size = 3, show.legend = FALSE)
  # nest
  
  means <- aggregate(ses_nest ~  set_name, SES_data, mean)
  
  nn <- ggplot(SES_data, aes(x = fct_relevel(set_name, "CINA1", "CINA3", "CINA2", "CINA4", "RUNA2"), y = ses_nest)) +
    geom_boxplot(color = c("#CC66CC", "#CC66CC", "#1B9E77", "#1B9E77", "#FF7F00"), fill = NA) +
    geom_jitter(aes(color = set_name), width = 0.2, size = 1.5, alpha = 0.6, show.legend = FALSE) +  # Points colorés selon les catégories
    scale_color_manual(values = c("CINA1" = "#CC66CC", "CINA3" = "#CC66CC", "CINA2" = "#1B9E77", "CINA4" = "#1B9E77", "RUNA2" = "#FF7F00")) +
    labs(title = "Nestedness",
         x = "Comparisons",
         y = "") +
    theme(legend.position = "none") +
    scale_x_discrete(labels = intra) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), 
          axis.title.x = element_blank(), 
          axis.title.y = element_blank(),
          axis.text.y = element_text(size = 12)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11), 
          axis.title.x = element_blank(), 
          axis.title.y = element_text(size = 12)) +
    geom_hline(yintercept = -1.96, colour = "red") +
    geom_hline(yintercept = 1.96, colour = "red") +
    geom_hline(yintercept = 0, colour = "darkgrey") +
    ylim(min = -4.5, max = 4.5) +
    stat_summary(fun = mean, colour = "darkred", geom = "point", 
                 shape = 18, size = 3, show.legend = FALSE)
  
  means <- aggregate(ses_bray ~  set_name, SES_data, mean)
  
  oo <- ggplot(SES_data, aes(x = fct_relevel(set_name, "CINA1", "CINA3", "CINA2", "CINA4", "RUNA2"), y = ses_bray)) +
    geom_boxplot(color = c("#CC66CC", "#CC66CC", "#1B9E77", "#1B9E77", "#FF7F00"), fill = NA) +
    geom_jitter(aes(color = set_name), width = 0.2, size = 1.5, alpha = 0.6, show.legend = FALSE) +  # Points colorés selon les catégories
    scale_color_manual(values = c("CINA1" = "#CC66CC", "CINA3" = "#CC66CC", "CINA2" = "#1B9E77", "CINA4" = "#1B9E77", "RUNA2" = "#FF7F00")) +
    labs(title = "Bray-Curtis dissimilarity",
         x = "Comparisons",
         y = "") +
    theme(legend.position = "none") +
    scale_x_discrete(labels = intra) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), 
          axis.title.x = element_blank(), 
          axis.title.y = element_blank(),
          axis.text.y = element_text(size = 12)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11), 
          axis.title.x = element_blank(), 
          axis.title.y = element_text(size = 12)) +
    geom_hline(yintercept = -1.96, colour = "red") +
    geom_hline(yintercept = 1.96, colour = "red") +
    geom_hline(yintercept = 0, colour = "darkgrey") +
    ylim(min = -4.5, max = 4.5) +
    stat_summary(fun = mean, colour = "darkred", geom = "point", 
                 shape = 18, size = 3, show.legend = FALSE)
  
  fin <- cowplot::plot_grid(oo, ll, mm, nn,
                            ncol = 4,
                            nrow = 1)
  
  path_to_boxplot_SES_UDOC <- paste0("outputs/null_model/boxplot_betadiv_microhab_UDOC_SES_16_10_2024.pdf")
  ggsave(filename =  path_to_boxplot_SES_UDOC, plot = fin, width = 18.5, height = 7.5)

  
  return (path_to_boxplot_SES_UDOC)
  
}
