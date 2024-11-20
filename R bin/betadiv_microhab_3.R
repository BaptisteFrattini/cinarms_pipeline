#' Mean all species cover of each plates of each arms 
#'
#' @param raw_data the path to the raw data file
#' @param arms_id the ID of the arms to subset for
#'
#' @return the path to the subseted derived data file
#' @export 

beta_microhab_3 <- function(meta_data){
  
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
  
  # unique_prefixes <- meta_red$prefixe
  # 
  # bray_results <- data.frame(NULL, NULL)
  # turnover_results <-data.frame(NULL, NULL)
  # nestedness_results <- data.frame(NULL, NULL)
  # jaccard_results <- data.frame(NULL, NULL)
  # 
  # 
  # for (prefix in unique_prefixes) {
  #   #prefix = "CINA1"
  #   subset_data <- subset(data, meta$arms == prefix)
  #   subset_meta_red <- subset(meta_red, meta_red$prefixe == prefix)
  #   
  #   mat.bray <- vegan::vegdist(subset_data, "bray")
  #   
  #   B.pair.pa <- betapart::beta.pair(vegan::decostand(subset_data, "pa"), index.family = "jaccard") 
  #   
  #   mat.turn <- B.pair.pa$beta.jtu
  #   mat.nest <- B.pair.pa$beta.jne
  #   mat.jacc <- B.pair.pa$beta.jac
  #   
  #   #Turnover
  #   
  #   df.turn <- melt(as.matrix(mat.turn), varnames = c("row", "col"))
  #   df.turn$row <- as.character(df.turn$row)
  #   df.turn$col <- as.character(df.turn$col)
  #   df.turn <- subset(df.turn, row != col)
  #   df.turn <- df.turn %>%
  #     mutate(
  #       row = as.character(row),  # Convert 'row' to character
  #       col = as.character(col),  # Convert 'col' to character
  #       combined = paste(pmin(row, col), pmax(row, col), sep = "_")
  #     )
  #   df.turn <- df.turn[!duplicated(df.turn$combined), ]
  #   df.turn <- df.turn[,-4]
  #   
  #   turnover_results[prefix,1] <- mean(df.turn$value)
  #   
  #   # Nestedness
  #   df.nest <- melt(as.matrix(mat.nest), varnames = c("row", "col"))
  #   df.nest$row <- as.character(df.nest$row)
  #   df.nest$col <- as.character(df.nest$col)
  #   df.nest <- subset(df.nest, row != col)
  #   df.nest <- df.nest %>%
  #     mutate(
  #       row = as.character(row),  # Convert 'row' to character
  #       col = as.character(col),  # Convert 'col' to character
  #       combined = paste(pmin(row, col), pmax(row, col), sep = "_")
  #     )
  #   df.nest <- df.nest[!duplicated(df.nest$combined), ]
  #   df.nest <- df.nest[,-4]
  #   
  #   nestedness_results[prefix,1] <- mean(df.nest$value)
  #   
  #   # Jaccard
  #   df.jacc <- melt(as.matrix(mat.jacc), varnames = c("row", "col"))
  #   df.jacc$row <- as.character(df.jacc$row)
  #   df.jacc$col <- as.character(df.jacc$col)
  #   df.jacc <- subset(df.jacc, row != col)
  #   df.jacc <- df.jacc %>%
  #     mutate(
  #       row = as.character(row),  # Convert 'row' to character
  #       col = as.character(col),  # Convert 'col' to character
  #       combined = paste(pmin(row, col), pmax(row, col), sep = "_")
  #     )
  #   df.jacc <- df.jacc[!duplicated(df.jacc$combined), ]
  #   df.jacc <- df.jacc[,-4]
  #   
  #   jaccard_results[prefix,1] <- mean(df.jacc$value)
  #   
  #   # Bray
  #   df.bray <- melt(as.matrix(mat.bray), varnames = c("row", "col"))
  #   df.bray$row <- as.character(df.bray$row)
  #   df.bray$col <- as.character(df.bray$col)
  #   df.bray <- subset(df.bray, row != col)
  #   df.bray <- df.bray %>%
  #     mutate(
  #       row = as.character(row),  # Convert 'row' to character
  #       col = as.character(col),  # Convert 'col' to character
  #       combined = paste(pmin(row, col), pmax(row, col), sep = "_")
  #     )
  #   df.bray <- df.bray[!duplicated(df.bray$combined), ]
  #   df.bray <- df.bray[,-4]
  #   
  #   bray_results[prefix,1] <- mean(df.bray$value)
  #   
  #   
  #   
  # }
  # 
  # 
  # results <- data.frame(bray = bray_results$V1, 
  #                       jacc = jaccard_results$V1, 
  #                       turn = turnover_results$V1, 
  #                       nest = nestedness_results$V1)
  # 
  # rownames(results) <- rownames(bray_results)
  # 
  #### Comparaison UDOC ####
  
  df_mean_microhab2 <- data %>% 
    group_by(meta$o_c, meta$Orientation, meta$arms_name) %>% 
    summarise_all(mean, na.rm = TRUE)
  
  set <- df_mean_microhab2$`meta$arms_name`
  
  df_mean_microhab2 <- as.data.frame(df_mean_microhab2)
  rownames(df_mean_microhab2) <- paste0(df_mean_microhab2$`meta$arms_name`, "_", df_mean_microhab2$`meta$Orientation`, df_mean_microhab2$`meta$o_c`)
  df_mean_microhab2 <- df_mean_microhab2[,-c(1,2,3)]
  
  ## Compute index on null model
  
  ## null model microhab ##
  num_iterations <- 1000  
  df_mean_microhab2_t_pa <- vegan::decostand(t(df_mean_microhab2), "pa")
  # null_model <- sim9(df_mean, nReps = num_iterations, metric = "c_score", algo = "sim9")
  null_model <- EcoSimR::sim9(df_mean_microhab2_t_pa, nReps = num_iterations, metric = "c_score", algo = "sim2")
  
  summary(null_model)
  null_model_data <- t(null_model$Randomized.Data)
  
  plot(null_model,type="burn_in")
  plot(null_model,type="hist")
  plot(null_model,type="cooc")
  
  df_mean.rounded <- round(df_mean_microhab2*10000, digits = 0)
  df_mean.rd.bray <- canaper::cpr_rand_comm(df_mean.rounded, "quasiswap_count", 10000)
  df_mean.rd.bray <- df_mean.rd.bray/10000

  
  df_mean.rd.bray <- as.data.frame(vegan::decostand(df_mean.rd.bray, "total"))*100

   
  ## Compute index on null model data- ##
  
  turnover_results <- as.data.frame(matrix(NA, nrow = 6*15, ncol = 1))
  turnover_results$names <- rep(sort(unique(meta$arms_name)), each = 6)
  
  nestedness_results <- as.data.frame(matrix(NA, nrow = 6*15, ncol = 1))
  nestedness_results$names <-  rep(sort(unique(meta$arms_name)), each = 6)
  
  jaccard_results <- as.data.frame(matrix(NA, nrow = 6*15, ncol = 1))
  jaccard_results$names <-  rep(sort(unique(meta$arms_name)), each = 6)
  
  bray_results <- as.data.frame(matrix(NA, nrow = 6*15, ncol = 1))
  bray_results$names <-  rep(sort(unique(meta$arms_name)), each = 6)
  
  
  unique_prefixes <- unique(meta_red$set)
  
  
  # Iterate through unique "prefixe" values
  
  for (prefix in unique_prefixes) {
    #prefix = "CINA1A"
    subset_data_bray <- subset(df_mean.rd.bray, set == prefix)
    mat.bray <- vegan::vegdist(subset_data_bray, "bray")
    
    subset_data <- subset(null_model_data, set == prefix)
    
    B.pair.pa <- betapart::beta.pair(vegan::decostand(subset_data, "pa"), index.family = "jaccard")
    
    mat.turn <- B.pair.pa$beta.jtu
    mat.nest <- B.pair.pa$beta.jne
    mat.jacc <- B.pair.pa$beta.jac
    
    
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
    df.turn <- df.turn[,-4]
    
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
    df.nest <- df.nest[,-4]
    
    nestedness_results$V1[nestedness_results$names == prefix] <- df.nest$value
    
    # Jaccard
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
    df.jacc <- df.jacc[,-4]
    
    jaccard_results$V1[jaccard_results$names == prefix] <- df.jacc$value
    
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
    df.bray <- df.bray[,-4]
    
    bray_results$V1[bray_results$names == prefix] <- df.bray$value
    
  }
  
  null_results <- data.frame(jacc = jaccard_results$V1, 
                             turn = turnover_results$V1, 
                             nest = nestedness_results$V1,
                             bray = bray_results$V1)
  
  rownames(null_results) <- rownames(jaccard_results)
  
  #### compute index on observed data ####
  
  turnover_results <- as.data.frame(matrix(NA, nrow = 6*15, ncol = 1))
  turnover_results$names <- rep(sort(unique(meta$arms_name)), each = 6)
  
  nestedness_results <- as.data.frame(matrix(NA, nrow = 6*15, ncol = 1))
  nestedness_results$names <-  rep(sort(unique(meta$arms_name)), each = 6)
  
  jaccard_results <- as.data.frame(matrix(NA, nrow = 6*15, ncol = 1))
  jaccard_results$names <-  rep(sort(unique(meta$arms_name)), each = 6)
  
  bray_results <- as.data.frame(matrix(NA, nrow = 6*15, ncol = 1))
  bray_results$names <-  rep(sort(unique(meta$arms_name)), each = 6)
  
  # Iterate through unique "prefixe" values
  
  for (prefix in unique_prefixes) {
    #prefix = "CINA1A"
    subset_data <- subset(df_mean_microhab2, set == prefix)
    mat.bray <- vegan::vegdist(subset_data, "bray")
    B.pair.pa <- betapart::beta.pair(vegan::decostand(subset_data, "pa"), index.family = "jaccard") 
    
    mat.turn <- B.pair.pa$beta.jtu
    mat.nest <- B.pair.pa$beta.jne
    mat.jacc <- B.pair.pa$beta.jac
    
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
    df.turn <- df.turn[,-4]
    
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
    df.nest <- df.nest[,-4]
    
    nestedness_results$V1[nestedness_results$names == prefix] <- df.nest$value
    
    # Jaccard
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
    df.jacc <- df.jacc[,-4]
    
    jaccard_results$V1[jaccard_results$names == prefix] <- df.jacc$value
    
    
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
    df.bray <- df.bray[,-4]
    
    bray_results$V1[bray_results$names == prefix] <- df.bray$value
    
    
  }
  
  
  results <- data.frame(bray = bray_results$V1, 
                        jacc = jaccard_results$V1, 
                        turn = turnover_results$V1, 
                        nest = nestedness_results$V1)
  
  rownames(results) <- rownames(bray_results)
  
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
  SES_data$imm_time <- c("6 months","6 months","6 months", 
                         "1 year", "1 year", "1 year", 
                         "6 months","6 months","6 months",
                         "1 year", "1 year", "1 year", 
                         "2 years","2 years","2 years")
  
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
  
  
  # turnover
  
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
  
  path_to_boxplot_SES_UDOC <- paste0("outputs/null_model/boxplot_betadiv_microhab_UDOC_SES_15_10_2024.pdf")
  ggsave(filename =  path_to_boxplot_SES_UDOC, plot = fin, width = 18.5, height = 7.5)
  # 
  # ### 6 moi 1y 2 ans 
  # intra_imm_time = c("between different \n microhab. of  the ARMS \n immersed 6 months", 
  #                    "between different \n microhab. of  the ARMS \n immersed 1 year",
  #                    "between different \n microhab. of  the ARMS \n immersed 2 years")
  # ll2 <- ggplot(SES_data, aes(x = fct_relevel(imm_time, "6 months", "1 year", "2 years"), y = ses_jacc)) +
  #   geom_boxplot(fill =  c("#CC66CC","#1B9E77","#FF7F00") ) +
  #   labs(title = "Jaccard dissimilarity",
  #        x = "Comparisons",
  #        y = "Standardized Effect Size (SES)") +
  #   theme(legend.position = "none") +
  #   scale_x_discrete(labels=intra_imm_time) +
  #   theme_classic() +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.title.x = element_blank(), axis.title.y = element_text(size=12)) +
  #   geom_hline(yintercept = -1.96, colour = "red")+
  #   geom_hline(yintercept = 1.96, colour = "red") +
  #   geom_hline(yintercept = 0, colour = "darkgrey") +
  #   ylim(min = -3, max = 3)
  # # turnover
  # mm2 <- ggplot(SES_data, aes(x = fct_relevel(imm_time, "6 months", "1 year", "2 years"), y = ses_turn)) +
  #   geom_boxplot(fill =  c("#CC66CC","#1B9E77","#FF7F00")) +
  #   labs(title = "Turnover component",
  #        x = "Comparisons",
  #        y = "Standardized Effect Size (SES)") +
  #   theme(legend.position = "none") +
  #   scale_x_discrete(labels=intra_imm_time) +
  #   theme_classic() +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.title.x = element_blank(), axis.title.y = element_text(size=12))+
  #   geom_hline(yintercept = -1.96, colour = "red")+
  #   geom_hline(yintercept = 1.96, colour = "red") +
  #   geom_hline(yintercept = 0, colour = "darkgrey") +
  #   ylim(min = -3, max = 3)
  # # nest
  # nn2 <- ggplot(SES_data, aes(x = fct_relevel(imm_time, "6 months", "1 year", "2 years"), y = ses_nest)) +
  #   geom_boxplot(fill =  c("#CC66CC","#1B9E77","#FF7F00")) +
  #   labs(title = "Nestedness component",
  #        x = "Comparisons",
  #        y = "Standardized Effect Size (SES)") +
  #   theme(legend.position = "none") +
  #   scale_x_discrete(labels=intra_imm_time) +
  #   theme_classic() +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.title.x = element_blank(), axis.title.y = element_text(size=12))+
  #   geom_hline(yintercept = -1.96, colour = "red")+
  #   geom_hline(yintercept = 1.96, colour = "red") +
  #   geom_hline(yintercept = 0, colour = "darkgrey") +
  #   ylim(min = -3, max = 3)
  # 
  # 
  # fin <- cowplot::plot_grid(ll2, mm2, nn2,
  #                           ncol = 3,
  #                           nrow = 1)
  # 
  # path_to_boxplot_SES_UDOC <- paste0("outputs/null_model/boxplot_betadiv_microhab_UDOC_SES_immtime.pdf")
  # ggsave(filename =  path_to_boxplot_SES_UDOC, plot = fin, width = 16, height = 7)
  # 
  # ###
  # 
  # #### Compute boxplots (not SES values) ####
  # 
  # #Jaccard
  # intra = c("between all microhabitat \n of the CINA1 set", 
  #           "between all microhabitat \n of the CINA3 set", 
  #           "between all microhabitat \n of the CINA2 set", 
  #           "between all microhabitat \n of the CINA4 set", 
  #           "between all microhabitat \n of the RUNA2 set")
  # 
  # set_name <- sort(rep(unique(substr(meta$arms_name,1,5)),3))
  # 
  # ww <- ggplot(results, aes(x = fct_relevel(set_name, "CINA1", "CINA3", "CINA2", "CINA4", "RUNA2"), y = jacc)) +
  #   geom_boxplot(fill =  c("#CC66CC","#CC66CC","#1B9E77","#1B9E77","#FF7F00") ) +
  #   labs(title = "",
  #        x = "Comparisons",
  #        y = "Jaccard dissimilarity") +
  #   theme(legend.position = "none") +
  #   scale_x_discrete(labels=intra) +
  #   theme_classic() +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.title.x = element_blank(), axis.title.y = element_text(size=12))
  # 
  # xx <- ggplot(results, aes(x = fct_relevel(set_name, "CINA1", "CINA3", "CINA2", "CINA4", "RUNA2"), y = turn)) +
  #   geom_boxplot(fill =  c("#CC66CC","#CC66CC","#1B9E77","#1B9E77","#FF7F00") ) +
  #   labs(title = "",
  #        x = "Comparisons",
  #        y = "Turnover component") +
  #   theme(legend.position = "none") +
  #   scale_x_discrete(labels=intra) +
  #   theme_classic() +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.title.x = element_blank(), axis.title.y = element_text(size=12))
  # 
  # yy <- ggplot(results, aes(x = fct_relevel(set_name, "CINA1", "CINA3", "CINA2", "CINA4", "RUNA2"), y = nest)) +
  #   geom_boxplot(fill =  c("#CC66CC","#CC66CC","#1B9E77","#1B9E77","#FF7F00") ) +
  #   labs(title = "",
  #        x = "Comparisons",
  #        y = "Nestedness component") +
  #   theme(legend.position = "none") +
  #   scale_x_discrete(labels=intra) +
  #   theme_classic() +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.title.x = element_blank(), axis.title.y = element_text(size=12))
  # 
  # zz <- ggplot(results, aes(x = fct_relevel(set_name, "CINA1", "CINA3", "CINA2", "CINA4", "RUNA2"), y = bray)) +
  #   geom_boxplot(fill =  c("#CC66CC","#CC66CC","#1B9E77","#1B9E77","#FF7F00") ) +
  #   labs(title = "",
  #        x = "Comparisons",
  #        y = "Bray-Curtis index") +
  #   theme(legend.position = "none") +
  #   scale_x_discrete(labels=intra) +
  #   theme_classic() +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.title.x = element_blank(), axis.title.y = element_text(size=12))
  # 
  # fin <- cowplot::plot_grid(ww, xx, yy, zz,
  #                           ncol = 2,
  #                           nrow = 2)
  # 
  # path_to_boxplot_unique_UDOC <- paste0("outputs/beta/boxplot_betadiv_microhab_UDOC_27_10.pdf")
  # ggsave(filename =  path_to_boxplot_unique_UDOC, plot = fin, width = 13.1, height = 12.9)
  # 
  # ### 6 moi 1y 2 ans (not SES)
  # intra_imm_time = c("between different \n microhab. of  the ARMS \n immersed 6 months", 
  #                    "between different \n microhab. of  the ARMS \n immersed 1 year",
  #                    "between different \n microhab. of  the ARMS \n immersed 2 years")
  # 
  # results$imm_time <- c("6 months","6 months","6 months", 
  #                        "1 year", "1 year", "1 year", 
  #                        "6 months","6 months","6 months",
  #                        "1 year", "1 year", "1 year", 
  #                        "2 years","2 years","2 years")
  # 
  # ww2 <- ggplot(results, aes(x = fct_relevel(imm_time, "6 months", "1 year", "2 years"), y = jacc)) +
  #   geom_boxplot(fill =  c("#CC66CC","#1B9E77","#FF7F00") ) +
  #   labs(title = "",
  #        x = "Comparisons",
  #        y = "Jaccard dissimilarity") +
  #   theme(legend.position = "none") +
  #   scale_x_discrete(labels=intra_imm_time) +
  #   theme_classic() +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.title.x = element_blank(), axis.title.y = element_text(size=12))
  # 
  # xx2 <- ggplot(results, aes(x = fct_relevel(imm_time, "6 months", "1 year", "2 years"), y = turn)) +
  #   geom_boxplot(fill =  c("#CC66CC","#1B9E77","#FF7F00") ) +
  #   labs(title = "",
  #        x = "Comparisons",
  #        y = "Turnover component") +
  #   theme(legend.position = "none") +
  #   scale_x_discrete(labels=intra_imm_time) +
  #   theme_classic() +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.title.x = element_blank(), axis.title.y = element_text(size=12))
  # 
  # yy2 <- ggplot(results, aes(x = fct_relevel(imm_time, "6 months", "1 year", "2 years"), y = nest)) +
  #   geom_boxplot(fill =  c("#CC66CC","#1B9E77","#FF7F00") ) +
  #   labs(title = "",
  #        x = "Comparisons",
  #        y = "Nestedness component") +
  #   theme(legend.position = "none") +
  #   scale_x_discrete(labels=intra_imm_time) +
  #   theme_classic() +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.title.x = element_blank(), axis.title.y = element_text(size=12))
  # 
  # zz2 <- ggplot(results, aes(x = fct_relevel(imm_time, "6 months", "1 year", "2 years"), y = bray)) +
  #   geom_boxplot(fill =  c("#CC66CC","#1B9E77","#FF7F00") ) +
  #   labs(title = "",
  #        x = "Comparisons",
  #        y = "Bray-Curtis dissimilarity") +
  #   theme(legend.position = "none") +
  #   scale_x_discrete(labels=intra_imm_time) +
  #   theme_classic() +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.title.x = element_blank(), axis.title.y = element_text(size=12))
  # 
  # fin <- cowplot::plot_grid(ww2, xx2, yy2, zz2,
  #                           ncol = 2,
  #                           nrow = 2)
  # 
  # path_to_boxplot_unique_UDOC <- paste0("outputs/beta/boxplot_betadiv_microhab_UDOC_27_10_imm_time.pdf")
  # ggsave(filename =  path_to_boxplot_unique_UDOC, plot = fin, width = 13.1, height = 12.9)
  #### comparaison UD --> Unique ####
###############################################################################
  # Create empty table to fill with B div values
  # df_mean_microhab <- data %>% 
  #   group_by(meta$Orientation,meta$arms_name) %>% 
  #   summarise_all(mean, na.rm = TRUE)
  # 
  # set <- substr(df_mean_microhab$`meta$arms_name`,1,6)
  # 
  # df_mean_microhab <- as.data.frame(df_mean_microhab)
  # rownames(df_mean_microhab) <- paste0(df_mean_microhab$`meta$arms_name`, "_", df_mean_microhab$`meta$Orientation`)
  # df_mean_microhab <- df_mean_microhab[,-c(1,2)]
  # unique_prefixes <- unique(set)
  # 
  # turnover_results <- as.data.frame(matrix(NA, nrow = 15, ncol = 1))
  # rownames(turnover_results) <- sort(unique(meta$arms_name))
  # 
  # nestedness_results <- as.data.frame(matrix(NA, nrow = 15, ncol = 1))
  # rownames(nestedness_results) <- sort(unique(meta$arms_name))
  # 
  # jaccard_results <- as.data.frame(matrix(NA, nrow = 15, ncol = 1))
  # rownames(jaccard_results) <- sort(unique(meta$arms_name))
  # 
  # bray_results <- as.data.frame(matrix(NA, nrow = 15, ncol = 1))
  # rownames(bray_results) <- sort(unique(meta$arms_name))
  # 
  # # Iterate through unique "prefixe" values
  # 
  # for (prefix in unique_prefixes) {
  #   #prefix = "CINA1A"
  #   subset_data <- subset(df_mean_microhab, set == prefix)
  #   mat.bray <- vegan::vegdist(subset_data, "bray")
  #   B.pair.pa <- betapart::beta.pair(vegan::decostand(subset_data, "pa"), index.family = "jaccard") 
  #   
  #   mat.turn <- B.pair.pa$beta.jtu
  #   mat.nest <- B.pair.pa$beta.jne
  #   mat.jacc <- B.pair.pa$beta.jac
  #   
  #   #Turnover
  #   
  #   df.turn <- melt(as.matrix(mat.turn), varnames = c("row", "col"))
  #   df.turn$row <- as.character(df.turn$row)
  #   df.turn$col <- as.character(df.turn$col)
  #   df.turn <- subset(df.turn, row != col)
  #   df.turn <- df.turn %>%
  #     mutate(
  #       row = as.character(row),  # Convert 'row' to character
  #       col = as.character(col),  # Convert 'col' to character
  #       combined = paste(pmin(row, col), pmax(row, col), sep = "_")
  #     )
  #   df.turn <- df.turn[!duplicated(df.turn$combined), ]
  #   df.turn <- df.turn[,-4]
  #   
  #   turnover_results[prefix,1] <- df.turn$value
  #   
  #   # Nestedness
  #   df.nest <- melt(as.matrix(mat.nest), varnames = c("row", "col"))
  #   df.nest$row <- as.character(df.nest$row)
  #   df.nest$col <- as.character(df.nest$col)
  #   df.nest <- subset(df.nest, row != col)
  #   df.nest <- df.nest %>%
  #     mutate(
  #       row = as.character(row),  # Convert 'row' to character
  #       col = as.character(col),  # Convert 'col' to character
  #       combined = paste(pmin(row, col), pmax(row, col), sep = "_")
  #     )
  #   df.nest <- df.nest[!duplicated(df.nest$combined), ]
  #   df.nest <- df.nest[,-4]
  #   
  #   nestedness_results[prefix,1] <- df.nest$value
  #   
  #   # Jaccard
  #   df.jacc <- melt(as.matrix(mat.jacc), varnames = c("row", "col"))
  #   df.jacc$row <- as.character(df.jacc$row)
  #   df.jacc$col <- as.character(df.jacc$col)
  #   df.jacc <- subset(df.jacc, row != col)
  #   df.jacc <- df.jacc %>%
  #     mutate(
  #       row = as.character(row),  # Convert 'row' to character
  #       col = as.character(col),  # Convert 'col' to character
  #       combined = paste(pmin(row, col), pmax(row, col), sep = "_")
  #     )
  #   df.jacc <- df.jacc[!duplicated(df.jacc$combined), ]
  #   df.jacc <- df.jacc[,-4]
  #   
  #   jaccard_results[prefix,1] <- df.jacc$value
  #   
  #   # Bray
  #   df.bray <- melt(as.matrix(mat.bray), varnames = c("row", "col"))
  #   df.bray$row <- as.character(df.bray$row)
  #   df.bray$col <- as.character(df.bray$col)
  #   df.bray <- subset(df.bray, row != col)
  #   df.bray <- df.bray %>%
  #     mutate(
  #       row = as.character(row),  # Convert 'row' to character
  #       col = as.character(col),  # Convert 'col' to character
  #       combined = paste(pmin(row, col), pmax(row, col), sep = "_")
  #     )
  #   df.bray <- df.bray[!duplicated(df.bray$combined), ]
  #   df.bray <- df.bray[,-4]
  #   
  #   bray_results[prefix,1] <- df.bray$value
  #   
  #   
  #   
  # }
  # 
  # 
  # results <- data.frame(bray = bray_results$V1, 
  #                       jacc = jaccard_results$V1, 
  #                       turn = turnover_results$V1, 
  #                       nest = nestedness_results$V1)
  # 
  # rownames(results) <- rownames(bray_results)
  # 
  # #Jaccard
  # intra = c("between upward \n and downward faces \n of the CINA1 set", 
  #           "between upward \n and downward faces \n of the CINA3 set", 
  #           "between upward \n and downward faces \n of the CINA2 set", 
  #           "between upward \n and downward faces \n of the CINA4 set", 
  #           "between upward \n and downward faces \n of the RUNA2 set")
  # 
  # set_name <- sort(rep(unique(substr(meta$arms_name,1,5)),3))
  # 
  # aa <- ggplot(results, aes(x = fct_relevel(set_name, "CINA1", "CINA3", "CINA2", "CINA4", "RUNA2"), y = jacc)) +
  #   geom_boxplot(fill =  c("#CC66CC","#CC66CC","#1B9E77","#1B9E77","#FF7F00") ) +
  #   labs(title = "",
  #        x = "Comparisons",
  #        y = "Jaccard dissimilarity") +
  #   theme(legend.position = "none") +
  #   scale_x_discrete(labels=intra) +
  #   theme_classic() +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.title.x = element_blank(), axis.title.y = element_text(size=12))
  # 
  # bb <- ggplot(results, aes(x = fct_relevel(set_name, "CINA1", "CINA3", "CINA2", "CINA4", "RUNA2"), y = turn)) +
  #   geom_boxplot(fill =  c("#CC66CC","#CC66CC","#1B9E77","#1B9E77","#FF7F00") ) +
  #   labs(title = "",
  #        x = "Comparisons",
  #        y = "Turnover component") +
  #   theme(legend.position = "none") +
  #   scale_x_discrete(labels=intra) +
  #   theme_classic() +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.title.x = element_blank(), axis.title.y = element_text(size=12))
  # 
  # cc <- ggplot(results, aes(x = fct_relevel(set_name, "CINA1", "CINA3", "CINA2", "CINA4", "RUNA2"), y = nest)) +
  #   geom_boxplot(fill =  c("#CC66CC","#CC66CC","#1B9E77","#1B9E77","#FF7F00") ) +
  #   labs(title = "",
  #        x = "Comparisons",
  #        y = "Nestedness component") +
  #   theme(legend.position = "none") +
  #   scale_x_discrete(labels=intra) +
  #   theme_classic() +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.title.x = element_blank(), axis.title.y = element_text(size=12))
  # 
  # dd <- ggplot(results, aes(x = fct_relevel(set_name, "CINA1", "CINA3", "CINA2", "CINA4", "RUNA2"), y = bray)) +
  #   geom_boxplot(fill =  c("#CC66CC","#CC66CC","#1B9E77","#1B9E77","#FF7F00") ) +
  #   labs(title = "",
  #        x = "Comparisons",
  #        y = "Bray-Curtis index") +
  #   theme(legend.position = "none") +
  #   scale_x_discrete(labels=intra) +
  #   theme_classic() +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.title.x = element_blank(), axis.title.y = element_text(size=12))
  # 
  # fin <- cowplot::plot_grid(aa, bb, cc, dd,
  #                           ncol = 2,
  #                           nrow = 2)
  # 
  # path_to_boxplot_unique <- paste0("outputs/beta/boxplot_betadiv_microhab_27_10.pdf")
  # ggsave(filename =  path_to_boxplot_unique, plot = fin, width = 12, height = 14)
  # 
  #### Comparaison UD ####
  #############################################################################
  # Initialize empty lists to store results
  
  # bray_results <- list()
  # jaccard_results <- list()
  # turnover_results <- list()
  # nestedness_results <- list()
  # 
  # # Iterate through unique "prefixe" values
  # unique_prefixes <- unique(meta_red$arms)
  # 
  # for (prefix in unique_prefixes) {
  #   # Subset the data for the current prefix
  #   #prefix = "CINA4A"
  #   subset_data <-  data[meta_red$arms == prefix, ]
  #   subset_meta <-  meta_red[meta_red$arms == prefix, ]
  #   row.names(subset_data) <- paste0(rownames(subset_data), subset_meta$Orientation)
  #   mat.bray <- vegan::vegdist(subset_data, "bray")
  # 
  #   B.pair.pa <- betapart::beta.pair(vegan::decostand(subset_data, "pa"), index.family = "jaccard")
  #   
  #   mat.turn <- B.pair.pa$beta.jtu
  #   mat.nest <- B.pair.pa$beta.jne
  #   mat.jacc <- B.pair.pa$beta.jac
  #   
  #   #Turnover
  #   df.turn <- melt(as.matrix(mat.turn), varnames = c("row", "col"))
  #   df.turn$row <- as.character(df.turn$row)
  #   df.turn$col <- as.character(df.turn$col)
  #   df.turn <- subset(df.turn, row != col)
  #   df.turn <- df.turn %>%
  #     mutate(
  #       row = as.character(row),  # Convert 'row' to character
  #       col = as.character(col),  # Convert 'col' to character
  #       combined = paste(pmin(row, col), pmax(row, col), sep = "")
  #     )
  #   df.turn <- df.turn[!duplicated(df.turn$combined), ]
  #   
  #   df.turn$row <- substr(df.turn$row, nchar(df.turn$row), nchar(df.turn$row))
  #   df.turn$col <- substr(df.turn$col, nchar(df.turn$col),nchar(df.turn$col))
  #   df.turn$same_value <- ifelse(df.turn$row == df.turn$col, "Yes", "No")
  #   
  #   #Nestedness
  #   df.nest <- melt(as.matrix(mat.nest), varnames = c("row", "col"))
  #   df.nest$row <- as.character(df.nest$row)
  #   df.nest$col <- as.character(df.nest$col)
  #   df.nest <- subset(df.nest, row != col)
  #   df.nest <- df.nest %>%
  #     mutate(
  #       row = as.character(row),  # Convert 'row' to character
  #       col = as.character(col),  # Convert 'col' to character
  #       combined = paste(pmin(row, col), pmax(row, col), sep = "")
  #     )
  #   df.nest <- df.nest[!duplicated(df.nest$combined), ]
  #   
  #   df.nest$row <- substr(df.nest$row, nchar(df.nest$row), nchar(df.nest$row))
  #   df.nest$col <- substr(df.nest$col, nchar(df.nest$col),nchar(df.nest$col))
  #   df.nest$same_value <- ifelse(df.nest$row == df.nest$col, "Yes", "No")
  #   
  #   #Jaccard
  #   df.jacc <- melt(as.matrix(mat.jacc), varnames = c("row", "col"))
  #   df.jacc$row <- as.character(df.jacc$row)
  #   df.jacc$col <- as.character(df.jacc$col)
  #   df.jacc <- subset(df.jacc, row != col)
  #   df.jacc <- df.jacc %>%
  #     mutate(
  #       row = as.character(row),  # Convert 'row' to character
  #       col = as.character(col),  # Convert 'col' to character
  #       combined = paste(pmin(row, col), pmax(row, col), sep = "")
  #     )
  #   df.jacc <- df.jacc[!duplicated(df.jacc$combined), ]
  #   
  #   df.jacc$row <- substr(df.jacc$row, nchar(df.jacc$row), nchar(df.jacc$row))
  #   df.jacc$col <- substr(df.jacc$col, nchar(df.jacc$col),nchar(df.jacc$col))
  #   df.jacc$same_value <- ifelse(df.jacc$row == df.jacc$col, "Yes", "No")
  #   
  #   # Bray Curtis
  #   df.bray <- melt(as.matrix(mat.bray), varnames = c("row", "col"))
  #   df.bray$row <- as.character(df.bray$row)
  #   df.bray$col <- as.character(df.bray$col)
  #   df.bray <- subset(df.bray, row != col)
  #   df.bray <- df.bray %>%
  #     mutate(
  #       row = as.character(row),  # Convert 'row' to character
  #       col = as.character(col),  # Convert 'col' to character
  #       combined = paste(pmin(row, col), pmax(row, col), sep = "")
  #     )
  #   df.bray <- df.bray[!duplicated(df.bray$combined), ]
  #   
  #   df.bray$row <- substr(df.bray$row, nchar(df.bray$row), nchar(df.bray$row))
  #   df.bray$col <- substr(df.bray$col, nchar(df.bray$col),nchar(df.bray$col))
  #   df.bray$same_value <- ifelse(df.bray$row == df.bray$col, "Yes", "No")
  #   
  #   # Store the results in the respective lists
  #   bray_results[[prefix]] <- df.bray
  #   jaccard_results[[prefix]] <- df.jacc
  #   turnover_results[[prefix]] <- df.turn
  #   nestedness_results[[prefix]] <-  df.nest
  # }
  # 
  # jaccard_results
  # turnover_results
  # nestedness_results
  # bray_results
  # 
  # df_jaccard <- data.frame(Same_value = jaccard_results$CINA1A$same_value,
  #                          CINA1A = jaccard_results$CINA1A$value, 
  #                          CINA1B = jaccard_results$CINA1B$value, 
  #                          CINA1C = jaccard_results$CINA1C$value,
  #                          CINA3A = jaccard_results$CINA3A$value, 
  #                          CINA3B = jaccard_results$CINA3B$value,
  #                          CINA3C = jaccard_results$CINA3C$value,
  #                          CINA2A = jaccard_results$CINA2A$value,
  #                          CINA2B = jaccard_results$CINA2B$value,
  #                          CINA2C = jaccard_results$CINA2C$value,
  #                          CINA4A = jaccard_results$CINA4A$value,
  #                          CINA4B = jaccard_results$CINA4B$value,
  #                          CINA4C = jaccard_results$CINA4C$value,
  #                          RUNA2A = jaccard_results$RUNA2A$value,
  #                          RUNA2B = jaccard_results$RUNA2B$value,
  #                          RUNA2C = jaccard_results$RUNA2C$value)
  # 
  # df_jaccard_No <- subset(df_jaccard, Same_value == "No")
  # df_jaccard_No <- melt(as.matrix(df_jaccard_No[,-1]))
  # colnames(df_jaccard_No) <- c("num", "set", "value")
  # df_jaccard_No$set <- substr(df_jaccard_No$set, 1, 5)
  # # df_jaccard_No <- tapply(df_jaccard_No$value, substr(df_jaccard_No$set, 1, 6), mean)
  # 
  # df_turnover <- data.frame(Same_value = turnover_results$CINA1A$same_value,
  #                          CINA1A = turnover_results$CINA1A$value, 
  #                          CINA1B = turnover_results$CINA1B$value, 
  #                          CINA1C = turnover_results$CINA1C$value,
  #                          CINA3A = turnover_results$CINA3A$value, 
  #                          CINA3B = turnover_results$CINA3B$value,
  #                          CINA3C = turnover_results$CINA3C$value,
  #                          CINA2A = turnover_results$CINA2A$value,
  #                          CINA2B = turnover_results$CINA2B$value,
  #                          CINA2C = turnover_results$CINA2C$value,
  #                          CINA4A = turnover_results$CINA4A$value,
  #                          CINA4B = turnover_results$CINA4B$value,
  #                          CINA4C = turnover_results$CINA4C$value,
  #                          RUNA2A = turnover_results$RUNA2A$value,
  #                          RUNA2B = turnover_results$RUNA2B$value,
  #                          RUNA2C = turnover_results$RUNA2C$value)
  # 
  # df_turnover_No <- subset(df_turnover, Same_value == "No")
  # df_turnover_No <- melt(as.matrix(df_turnover_No[,-1]))
  # colnames(df_turnover_No) <- c("num", "set", "value")
  # df_turnover_No$set <- substr(df_turnover_No$set, 1, 5)
  # 
  # 
  # df_nestedness <- data.frame(Same_value = nestedness_results$CINA1A$same_value,
  #                          CINA1A = nestedness_results$CINA1A$value, 
  #                          CINA1B = nestedness_results$CINA1B$value, 
  #                          CINA1C = nestedness_results$CINA1C$value,
  #                          CINA3A = nestedness_results$CINA3A$value, 
  #                          CINA3B = nestedness_results$CINA3B$value,
  #                          CINA3C = nestedness_results$CINA3C$value,
  #                          CINA2A = nestedness_results$CINA2A$value,
  #                          CINA2B = nestedness_results$CINA2B$value,
  #                          CINA2C = nestedness_results$CINA2C$value,
  #                          CINA4A = nestedness_results$CINA4A$value,
  #                          CINA4B = nestedness_results$CINA4B$value,
  #                          CINA4C = nestedness_results$CINA4C$value,
  #                          RUNA2A = nestedness_results$RUNA2A$value,
  #                          RUNA2B = nestedness_results$RUNA2B$value,
  #                          RUNA2C = nestedness_results$RUNA2C$value)
  # 
  # df_nestedness_No <- subset(df_nestedness, Same_value == "No")
  # df_nestedness_No <- melt(as.matrix(df_nestedness_No[,-1]))
  # colnames(df_nestedness_No) <- c("num", "set", "value")
  # df_nestedness_No$set <- substr(df_nestedness_No$set, 1, 5)
  # 
  # df_bray <- data.frame(Same_value = bray_results$CINA1A$same_value,
  #                          CINA1A = bray_results$CINA1A$value, 
  #                          CINA1B = bray_results$CINA1B$value, 
  #                          CINA1C = bray_results$CINA1C$value,
  #                          CINA3A = bray_results$CINA3A$value, 
  #                          CINA3B = bray_results$CINA3B$value,
  #                          CINA3C = bray_results$CINA3C$value,
  #                          CINA2A = bray_results$CINA2A$value,
  #                          CINA2B = bray_results$CINA2B$value,
  #                          CINA2C = bray_results$CINA2C$value,
  #                          CINA4A = bray_results$CINA4A$value,
  #                          CINA4B = bray_results$CINA4B$value,
  #                          CINA4C = bray_results$CINA4C$value,
  #                          RUNA2A = bray_results$RUNA2A$value,
  #                          RUNA2B = bray_results$RUNA2B$value,
  #                          RUNA2C = bray_results$RUNA2C$value)
  # 
  # df_bray_No <- subset(df_bray, Same_value == "No")
  # df_bray_No <- melt(as.matrix(df_bray_No[,-1]))
  # colnames(df_bray_No) <- c("num", "set", "value")
  # df_bray_No$set <- substr(df_bray_No$set, 1, 5)
  # 
  # #Jaccard
  # intra = c("between upward \n and downward faces \n of the CINA1 set", 
  #           "between upward \n and downward faces \n of the CINA3 set", 
  #           "between upward \n and downward faces \n of the CINA2 set", 
  #           "between upward \n and downward faces \n of the CINA4 set", 
  #           "between upward \n and downward faces \n of the RUNA2 set")
  # 
  # kk <- ggplot(df_jaccard_No, aes(x = fct_relevel(set, "CINA1", "CINA3", "CINA2", "CINA4", "RUNA2"), y = value)) +
  #   geom_boxplot(fill =  c("#CC66CC","#CC66CC","#1B9E77","#1B9E77","#FF7F00") ) +
  #   labs(title = "",
  #        x = "Comparisons",
  #        y = "Jaccard dissimilarity") +
  #   theme(legend.position = "none") +
  #   scale_x_discrete(labels=intra) +
  #   theme_classic() +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.title.x = element_blank(), axis.title.y = element_text(size=12))
  # 
  # gg <- ggplot(df_bray_No, aes(x = fct_relevel(set, "CINA1", "CINA3", "CINA2", "CINA4", "RUNA2"), y = value)) +
  #   geom_boxplot(fill =  c("#CC66CC","#CC66CC","#1B9E77","#1B9E77","#FF7F00") ) +
  #   labs(title = "",
  #        x = "Comparisons",
  #        y = "bray component") +
  #   theme(legend.position = "none") +
  #   scale_x_discrete(labels=intra) +
  #   theme_classic() +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.title.x = element_blank(), axis.title.y = element_text(size=12))
  # 
  # dd <- ggplot(df_turnover_No, aes(x = fct_relevel(set, "CINA1", "CINA3", "CINA2", "CINA4", "RUNA2"), y = value)) +
  #   geom_boxplot(fill =  c("#CC66CC","#CC66CC","#1B9E77","#1B9E77","#FF7F00") ) +
  #   labs(title = "",
  #        x = "Comparisons",
  #        y = "Turnover component") +
  #   theme(legend.position = "none") +
  #   scale_x_discrete(labels=intra) +
  #   theme_classic() +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.title.x = element_blank(), axis.title.y = element_text(size=12))
  # 
  # ss <- ggplot(df_nestedness_No, aes(x = fct_relevel(set, "CINA1", "CINA3", "CINA2", "CINA4", "RUNA2"), y = value)) +
  #   geom_boxplot(fill =  c("#CC66CC","#CC66CC","#1B9E77","#1B9E77","#FF7F00") ) +
  #   labs(title = "",
  #        x = "Comparisons",
  #        y = "Nestedness component") +
  #   theme(legend.position = "none") +
  #   scale_x_discrete(labels=intra) +
  #   theme_classic() +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.title.x = element_blank(), axis.title.y = element_text(size=12))
  #   
  # fin <- cowplot::plot_grid(kk, gg, dd, ss,
  #                           ncol = 2,
  #                           nrow = 2)
  # 
  # path_to_boxplot <- paste0("outputs/beta/boxplot_betadiv_microhab.pdf")
  # ggsave(filename =  path_to_boxplot, plot = fin, width = 12, height = 14)
  # 
  # 
  # 
  

  
  return (path_to_boxplot_SES_UDOC)
  
}