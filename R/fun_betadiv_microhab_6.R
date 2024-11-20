#' Null model
#'
#' @param metadata_data_mean data
#' @return the path to the subseted raw data file
#' @export

fun_beta_microhab <- function(meta_data){
  
  # Import data and packages ####
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
  library(ggplot2)
  library(reshape2)
  library(ggpubr)
  library(dplyr)
  library(effectsize)

  meta = read.csv(meta_data[grepl("metadata", meta_data)], header = TRUE)
  data = read.csv(meta_data[!grepl("metadata", meta_data)], header = TRUE)
  data <- data[ , colSums(data) != 0]
  meta <- meta %>%
    mutate(o_c = recode(o_c, "F" = "C"))
  meta_red <- data.frame(prefixe = substr(meta$arms_name,1,5), 
                         set = meta$arms_name, 
                         Orientation = meta$Orientation,
                         Open_Close = meta$o_c)
  
  
  
  rownames(data) <- paste0(meta$arms_name,meta$Orientation,meta$o_c,meta$plate_number)
  
  unique_prefixes <- meta_red$set
  
  # Compute series of observed index ####
  ## Computing beta index for comparisons btw =/= microhab ####
  jaccard_results <- as.data.frame(matrix(NA, nrow = (96/6)*15, ncol = 7))
  colnames(jaccard_results) <- c("names","DC_DO", "DC_UC", "DC_UO", "DO_UC", "DO_UO", "UC_UO")
  jaccard_results$names <-  rep(sort(unique(meta$arms_name)), each = 16)
  
  turnover_results <- as.data.frame(matrix(NA, nrow = (96/6)*15, ncol = 7))
  colnames(turnover_results) <- c("names","DC_DO", "DC_UC", "DC_UO", "DO_UC", "DO_UO", "UC_UO")
  turnover_results$names <-  rep(sort(unique(meta$arms_name)), each = 16)
  
  
  nestedness_results <- as.data.frame(matrix(NA, nrow = (96/6)*15, ncol = 7))
  colnames(nestedness_results) <- c("names","DC_DO", "DC_UC", "DC_UO", "DO_UC", "DO_UO", "UC_UO")
  nestedness_results$names <-  rep(sort(unique(meta$arms_name)), each = 16)
 
  
  bray_results <- as.data.frame(matrix(NA, nrow = (96/6)*15, ncol = 7))
  colnames(bray_results) <- c("names","DC_DO", "DC_UC", "DC_UO", "DO_UC", "DO_UO", "UC_UO")
  bray_results$names <-  rep(sort(unique(meta$arms_name)), each = 16)
  
  
  unique_prefixes <- unique(meta$arms_name)
  unique_comp <- c("DC_DO", "DC_UC", "DC_UO", "DO_UC", "DO_UO", "UC_UO")
  
  
  for (prefix in unique_prefixes) {
    #prefix = "CINA3A"
    for (comp in unique_comp) {
      #comp = "DC_DO"
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
    df.jacc <- df.jacc %>%
      mutate(
        row = as.character(row),  # Convert 'row' to character
        col = as.character(col),  # Convert 'col' to character
        combined = paste(pmin(row, col), pmax(row, col), sep = "_")
      )
    
    unique_comp = c(levels(as.factor(df.jacc$combined)))
    
    df <- subset(df.jacc, combined == comp)
    jaccard_results[,comp][jaccard_results$names == prefix] <- df$value 
    
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
    df.turn <- df.turn %>%
      mutate(
        row = as.character(row),  # Convert 'row' to character
        col = as.character(col),  # Convert 'col' to character
        combined = paste(pmin(row, col), pmax(row, col), sep = "_")
      )
    
    unique_comp = c(levels(as.factor(df.turn$combined)))
    
    df <- subset(df.turn, combined == comp)
    turnover_results[,comp][turnover_results$names == prefix] <- df$value 
    
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
    df.nest <- df.nest %>%
      mutate(
        row = as.character(row),  # Convert 'row' to character
        col = as.character(col),  # Convert 'col' to character
        combined = paste(pmin(row, col), pmax(row, col), sep = "_")
      )
    
    unique_comp = c(levels(as.factor(df.nest$combined)))
    
    df <- subset(df.nest, combined == comp)
    nestedness_results[,comp][nestedness_results$names == prefix] <- df$value 
    
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
    df.bray <- df.bray %>%
      mutate(
        row = as.character(row),  # Convert 'row' to character
        col = as.character(col),  # Convert 'col' to character
        combined = paste(pmin(row, col), pmax(row, col), sep = "_")
      )
    
    unique_comp = c(levels(as.factor(df.bray$combined)))
    
    df <- subset(df.bray, combined == comp)
    bray_results[,comp][bray_results$names == prefix] <- df$value 
    
    }
    
  }

  
  jaccard_long <- melt(jaccard_results, 
                       id.vars = "names", 
                       variable.name = "Comparison", 
                       value.name = "Value")
  
  turnover_long <- melt(turnover_results, 
                       id.vars = "names", 
                       variable.name = "Comparison", 
                       value.name = "Value")
  
  nestedness_long <- melt(nestedness_results, 
                       id.vars = "names", 
                       variable.name = "Comparison", 
                       value.name = "Value")
  
  bray_long <- melt(bray_results, 
                       id.vars = "names", 
                       variable.name = "Comparison", 
                       value.name = "Value")
  
  ## Computing beta index for comparisons btw == microhab ####
  jaccard_results_same <- as.data.frame(matrix(NA, nrow = 48*15, ncol = 3))
  colnames(jaccard_results_same) <- c("names", "Comparison", "Value")
  jaccard_results_same$names <-  rep(sort(unique(meta$arms_name)), each = 48)
  jaccard_results_same$Comparison <- rep("same", 48*15)
  
  turnover_results_same <- as.data.frame(matrix(NA, nrow = 48*15, ncol = 3))
  colnames(turnover_results_same) <- c("names","Comparison", "Value")
  turnover_results_same$names <-  rep(sort(unique(meta$arms_name)), each = 48)
  turnover_results_same$Comparison <- rep("same", 48*15)
  
  nestedness_results_same <- as.data.frame(matrix(NA, nrow = 48*15, ncol = 3))
  colnames(nestedness_results_same) <- c("names","Comparison", "Value")
  nestedness_results_same$names <-  rep(sort(unique(meta$arms_name)), each = 48)
  nestedness_results_same$Comparison <- rep("same", 48*15)
  
  bray_results_same <- as.data.frame(matrix(NA, nrow = 48*15, ncol = 3))
  colnames(bray_results_same) <- c("names", "Comparison", "Value")
  bray_results_same$names <-  rep(sort(unique(meta$arms_name)), each = 48)
  bray_results_same$Comparison <- rep("same", 48*15)
  
  
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
    
    jaccard_results_same$Value[jaccard_results_same$names == prefix] <- df.jacc$value
    
    
    #Turnover
    
    df.turn <- melt(as.matrix(mat.turn), varnames = c("row", "col"))
    df.turn$row <- as.character(df.turn$row)
    df.turn$col <- as.character(df.turn$col)
    df.turn <- subset(df.turn, row != col)
    df.turn$row <- substr(df.turn$row, 7, 8)
    df.turn$col <- substr(df.turn$col, 7, 8)
    df.turn <- subset(df.turn, row == col)
    
    turnover_results_same$Value[turnover_results_same$names == prefix] <- df.turn$value
    
    # Nestedness
    df.nest <- melt(as.matrix(mat.nest), varnames = c("row", "col"))
    df.nest$row <- as.character(df.nest$row)
    df.nest$col <- as.character(df.nest$col)
    df.nest <- subset(df.nest, row != col)
    df.nest$row <- substr(df.nest$row, 7, 8)
    df.nest$col <- substr(df.nest$col, 7, 8)
    df.nest <- subset(df.nest, row == col)
    
    nestedness_results_same$Value[nestedness_results_same$names == prefix] <- df.nest$value
    
    # Bray
    df.bray <- melt(as.matrix(mat.bray), varnames = c("row", "col"))
    df.bray$row <- as.character(df.bray$row)
    df.bray$col <- as.character(df.bray$col)
    df.bray <- subset(df.bray, row != col)
    df.bray$row <- substr(df.bray$row, 7, 8)
    df.bray$col <- substr(df.bray$col, 7, 8)
    df.bray <- subset(df.bray, row == col)
    
    bray_results_same$Value[bray_results_same$names == prefix] <- df.bray$value
    
  }
  

  jaccard_long <- as.data.frame(rbind(jaccard_long, jaccard_results_same))
  jaccard_long$Value <- as.numeric(jaccard_long$Value)
  jaccard_long$triplicat <- substr(jaccard_long$names, 1,5)
  jaccard_long$triplicat <- factor(jaccard_long$triplicat, levels = c("CINA1", "CINA3","RUNA2", "CINA2", "CINA4"))
  
  box_colors <- c("CINA1" = "#CC66CC", "CINA3" = "#CC66CC", "CINA2" = "#1B9E77", "CINA4" = "#1B9E77", "RUNA2" = "#FF7F00")
  
  jaccard_long_filtered <- jaccard_long %>%
    filter(Comparison == "DC_UO")
  
  jaccard_long_filtered$triplicat <- factor(jaccard_long_filtered$triplicat, 
                                             levels = c("CINA1", "CINA3", "CINA2", "CINA4", "RUNA2"))
  
  ggplot(jaccard_long_filtered, aes(x = triplicat, y = Value, fill = triplicat)) +
    geom_boxplot() +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3, color = "red", fill = "red") +
    scale_y_continuous(limits = c(0, 1)) +
    scale_fill_manual(values = box_colors) +
    labs(title = "Boxplot of DC_UO Values by Triplicat", y = "Jaccard dissimilarity", x = "Triplicat") +
    theme_minimal() +
    theme(legend.position = "none")  # Remove legend if not needed
  
  
  turnover_long <- as.data.frame(rbind(turnover_long, turnover_results_same))
  turnover_long$Value <- as.numeric(turnover_long$Value)
  turnover_long$triplicat <- substr(turnover_long$names, 1,5)
  turnover_long$triplicat <- factor(turnover_long$triplicat, levels = c("CINA1", "CINA3","RUNA2", "CINA2", "CINA4"))
    
  turnover_long_filtered <- turnover_long %>%
    filter(Comparison == "DC_UO")
  turnover_long_filtered$triplicat <- factor(turnover_long_filtered$triplicat, 
                                            levels = c("CINA1", "CINA3", "CINA2", "CINA4", "RUNA2"))
  ggplot(turnover_long_filtered, aes(x = triplicat, y = Value, fill = triplicat)) +
    geom_boxplot() +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3, color = "red", fill = "red") +
    scale_y_continuous(limits = c(0, 1)) +
    scale_fill_manual(values = box_colors) +
    labs(title = "Boxplot of DC_UO Values by Triplicat", y = "Turnover component", x = "Triplicat") +
    theme_minimal() +
    theme(legend.position = "none")  # Remove legend if not needed
  
  nestedness_long <- as.data.frame(rbind(nestedness_long, nestedness_results_same))
  nestedness_long$Value <- as.numeric(nestedness_long$Value)
  nestedness_long$triplicat <- substr(nestedness_long$names, 1,5)
  nestedness_long$triplicat <- factor(nestedness_long$triplicat, levels = c("CINA1", "CINA3","RUNA2", "CINA2", "CINA4"))
  
  nestedness_long_filtered <- nestedness_long %>%
    filter(Comparison == "DC_UO")
  nestedness_long_filtered$triplicat <- factor(nestedness_long_filtered$triplicat, 
                                            levels = c("CINA1", "CINA3", "CINA2", "CINA4", "RUNA2"))
  ggplot(nestedness_long_filtered, aes(x = triplicat, y = Value, fill = triplicat)) +
    geom_boxplot() +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3, color = "red", fill = "red") +
    scale_y_continuous(limits = c(0, 1)) +
    scale_fill_manual(values = box_colors) +
    labs(title = "Boxplot of DC_UO Values by Triplicat", y = "Nestedness component", x = "Triplicat") +
    theme_minimal() +
    theme(legend.position = "none")  # Remove legend if not needed
  
  bray_long <- as.data.frame(rbind(bray_long, bray_results_same))
  bray_long$Value <- as.numeric(bray_long$Value)
  bray_long$triplicat <- substr(bray_long$names, 1,5)
  bray_long$triplicat <- factor(bray_long$triplicat, levels = c("CINA1", "CINA3","RUNA2", "CINA2", "CINA4"))
  
  bray_long_filtered <- bray_long %>%
    filter(Comparison == "DC_UO")
  bray_long_filtered$triplicat <- factor(bray_long_filtered$triplicat, 
                                            levels = c("CINA1", "CINA3", "CINA2", "CINA4", "RUNA2"))
  ggplot(bray_long_filtered, aes(x = triplicat, y = Value, fill = triplicat)) +
    geom_boxplot() +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3, color = "red", fill = "red") +
    scale_y_continuous(limits = c(0, 1)) +
    scale_fill_manual(values = box_colors) +
    labs(title = "Boxplot of DC_UO Values by Triplicat", y = "Bray-Curtis dissimilarity", x = "Triplicat") +
    theme_minimal() +
    theme(legend.position = "none")  # Remove legend if not needed
  
  ## Testing effect of immersion time on difference between microhab with PERMANOVA ####
  meta$microhab <- as.factor(paste0(meta$o_c, "_", meta$Orientation))
  meta$triplicat <- as.factor(paste0(meta$imm_time, "_", meta$immersion_season))
  meta$microhab <- as.factor(meta$microhab)
  meta$imm_time <- as.factor(meta$imm_time)
  meta$immersion_season <- as.factor(meta$immersion_season)

  # perm.b <- vegan::adonis2(data ~ microhab * imm_time, data = meta, method = "bray", strata = meta$arms_name, permutations = 999, by = "terms")
  # 
  # perm.j <- vegan::adonis2(data ~ microhab * imm_time, data = meta, method = "jaccard", strata = meta$arms_name, permutations = 999, by = "terms")
  
  ## Plot graph of beta index comparisons btw same and differents ####
  
  # ### jaccard ####
  # #### Testing normality of data for Jaccard ####
  # 
  # step = 0.05*(max(jaccard_long$Value)-min(jaccard_long$Value))
  # brek = seq(min(jaccard_long$Value),max(jaccard_long$Value),  step)
  # 
  # j.norm = ggplot(jaccard_long, aes(x = Value)) +
  #   geom_histogram(breaks = brek, fill = "coral", color = "black", alpha = 0.7) +
  #   facet_wrap(~ triplicat, scales = "free_x") +
  #   labs(title = paste0("jacc ", prefix), x = "Beta-div index", y = "Frequency") +
  #   theme_minimal()
  # 
  # path_to_boxplot <- paste0("outputs/beta/boxplot_norm_J_30_10_24.pdf")
  # ggsave(filename =  path_to_boxplot, plot = j.norm, width = 12, height = 9)
  # 
  # #### Testing mean difference between comparisons ####
  # stat_test <- jaccard_long %>%
  #   group_by(triplicat) %>%
  #   rstatix::pairwise_t_test(Value ~ Comparison, p.adjust.method = "bonferroni") %>%
  #   filter(group2 == "same")  #
  # 
  # stat_test <- stat_test %>%
  #   mutate(y.position = rep(c(0.8, 0.835, 0.870, 0.905, 0.940, 0.975), 5))  
  # 
  # print(stat_test, n = 44) 
  # 
  # #### Plot the boxplot of beta div ####
  # j.beta <- ggplot(jaccard_long, aes(x = Comparison, y = Value)) +
  #           geom_boxplot(outlier.shape = NA) +  # Masquer les outliers pour ne pas interférer avec la moyenne
  #           stat_summary(fun = mean, geom = "point", shape = 18, color = "red", size = 3) +  # Ajouter la moyenne avec un losange rouge
  #           facet_wrap(~ triplicat, scales = "free_x") +  # Garder seulement l'échelle libre sur l'axe des x
  #           labs(title = "Boxplot of Jaccard Index Values by Comparison and Names",
  #                x = "Comparison", y = "Jaccard Index Value") +
  #           theme_minimal() +
  #           theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #           ylim(0, 1) +  # Définir l'échelle de l'axe des y entre 0 et 1
  #           stat_pvalue_manual(
  #             stat_test, 
  #             label = "p.adj.signif", 
  #             y.position = "y.position",  # Utiliser la colonne y.position de stat_test pour placer les segments
  #             bracket.size = 0.3
  #           )
  # path_to_boxplot <- paste0("outputs/beta/boxplot_betadiv_microhab_J_30_10_24.pdf")
  # ggsave(filename =  path_to_boxplot, plot = j.beta, width = 16, height = 11)
  # 
  # #### Compute the cohen's d ####
  # results <-jaccard_long %>%
  #           group_by(triplicat) %>%
  #           rstatix::pairwise_t_test(Value ~ Comparison, p.adjust.method = "bonferroni") %>%
  #           filter(group2 == "same")  #
  #         
  # 
  # cohen_results <- data.frame()
  # 
  # # Initialiser un dataframe pour stocker les résultats
  # d_cohen_results <- jaccard_long %>%
  #     group_by(triplicat) %>%
  #     summarise(
  #       Comparison_Pairs = list(combn(unique(Comparison), 2, simplify = FALSE)),
  #       .groups = 'drop'
  #     )
  # 
  # 
  # for (i in 1:nrow(d_cohen_results)) {
  #   triplicat_name <- d_cohen_results$triplicat[i]
  #   pairs <- d_cohen_results$Comparison_Pairs[[i]]
  #   
  #   for (pair in pairs) {
  #     values1 <- jaccard_long$Value[jaccard_long$Comparison == pair[1] & jaccard_long$triplicat == triplicat_name]
  #     values2 <- jaccard_long$Value[jaccard_long$Comparison == pair[2] & jaccard_long$triplicat == triplicat_name]
  #     
  #     # Calculer le d de Cohen
  #     if (length(values1) > 1 & length(values2) > 1) {
  #       d_cohen_value <- effectsize::cohens_d(values1, values2)
  #       cohen_results <- rbind(cohen_results, data.frame(triplicat = triplicat_name, Comparison1 = pair[1], Comparison2 = pair[2], d_cohen = d_cohen_value))
  #     }
  #   }
  # }
  # 
  # colnames(results) <-  c("triplicat", ".y.", "Comparison1", "Comparison2", "n1", "n2",  "p", "p.signif", "p.adj", "p.adj.signif")
  # 
  # 
  # final_results <- results %>%
  #   left_join(cohen_results, by = c("triplicat", "Comparison1", "Comparison2")) %>%
  #   filter(Comparison2 == "same") %>%
  #   filter(p.adj < 0.05)
  # 
  # #### plot the cohen's d ####
  # j.coh <- ggplot(final_results,aes(x = triplicat, y = d_cohen.Cohens_d)) +
  #             geom_boxplot(fill = "grey", outlier.colour = "black", outlier.shape = 16, outlier.size = 2) +
  #             stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "red", position = position_dodge(width = 0.75)) +
  #             labs(title = "",
  #                  x = "Comparisons",
  #                  y = "Cohen's d") +
  #             theme_minimal() +
  #             theme(axis.text.x = element_text(angle = 45, hjust = 1))
  #           
  # path_to_boxplot <- paste0("outputs/beta/boxplot_cohen_J_30_10_24.pdf")
  # ggsave(filename =  path_to_boxplot, plot = j.coh, width = 5.5, height = 4.5)
  # 
  # ### turnover ####
  # #### Testing normality of data for turnover ####
  # 
  # step = 0.05*(max(turnover_long$Value)-min(turnover_long$Value))
  # brek = seq(min(turnover_long$Value),max(turnover_long$Value),  step)
  # 
  # t.norm = ggplot(turnover_long, aes(x = Value)) +
  #   geom_histogram(breaks = brek, fill = "coral", color = "black", alpha = 0.7) +
  #   facet_wrap(~ triplicat, scales = "free_x") +
  #   labs(title = paste0("turn", prefix), x = "Beta-div index", y = "Frequency") +
  #   theme_minimal()
  # 
  # path_to_boxplot <- paste0("outputs/beta/boxplot_norm_T_30_10_24.pdf")
  # ggsave(filename =  path_to_boxplot, plot = t.norm, width = 12, height = 9)
  # 
  # #### Testing mean difference between comparisons ####
  # stat_test <- turnover_long %>%
  #   group_by(triplicat) %>%
  #   rstatix::pairwise_t_test(Value ~ Comparison, p.adjust.method = "bonferroni") %>%
  #   filter(group2 == "same")  #
  # 
  # stat_test <- stat_test %>%
  #   mutate(y.position = rep(c(0.8, 0.835, 0.870, 0.905, 0.940, 0.975), 5))  
  # 
  # print(stat_test, n = 44) 
  # 
  # #### Plot the boxplot of beta div ####
  # t.beta <- ggplot(turnover_long, aes(x = Comparison, y = Value)) +
  #   geom_boxplot(outlier.shape = NA) +  # Masquer les outliers pour ne pas interférer avec la moyenne
  #   stat_summary(fun = mean, geom = "point", shape = 18, color = "red", size = 3) +  # Ajouter la moyenne avec un losange rouge
  #   facet_wrap(~ triplicat, scales = "free_x") +  # Garder seulement l'échelle libre sur l'axe des x
  #   labs(title = "Boxplot of turnover Index Values by Comparison and Names",
  #        x = "Comparison", y = "turnover Index Value") +
  #   theme_minimal() +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #   ylim(0, 1) +  # Définir l'échelle de l'axe des y entre 0 et 1
  #   stat_pvalue_manual(
  #     stat_test, 
  #     label = "p.adj.signif", 
  #     y.position = "y.position",  # Utiliser la colonne y.position de stat_test pour placer les segments
  #     bracket.size = 0.3
  #   )
  # 
  # path_to_boxplot <- paste0("outputs/beta/boxplot_betadiv_microhab_T_30_10_24.pdf")
  # ggsave(filename =  path_to_boxplot, plot = t.beta, width = 16, height = 11)
  # 
  # #### Compute the cohen's d ####
  # results <-turnover_long %>%
  #   group_by(triplicat) %>%
  #   rstatix::pairwise_t_test(Value ~ Comparison, p.adjust.method = "bonferroni") %>%
  #   filter(group2 == "same")  #
  # 
  # 
  # cohen_results <- data.frame()
  # 
  # # Initialiser un dataframe pour stocker les résultats
  # d_cohen_results <- turnover_long %>%
  #   group_by(triplicat) %>%
  #   summarise(
  #     Comparison_Pairs = list(combn(unique(Comparison), 2, simplify = FALSE)),
  #     .groups = 'drop'
  #   )
  # 
  # 
  # for (i in 1:nrow(d_cohen_results)) {
  #   triplicat_name <- d_cohen_results$triplicat[i]
  #   pairs <- d_cohen_results$Comparison_Pairs[[i]]
  #   
  #   for (pair in pairs) {
  #     values1 <- turnover_long$Value[turnover_long$Comparison == pair[1] & turnover_long$triplicat == triplicat_name]
  #     values2 <- turnover_long$Value[turnover_long$Comparison == pair[2] & turnover_long$triplicat == triplicat_name]
  #     
  #     # Calculer le d de Cohen
  #     if (length(values1) > 1 & length(values2) > 1) {
  #       d_cohen_value <- effectsize::cohens_d(values1, values2)
  #       cohen_results <- rbind(cohen_results, data.frame(triplicat = triplicat_name, Comparison1 = pair[1], Comparison2 = pair[2], d_cohen = d_cohen_value))
  #     }
  #   }
  # }
  # 
  # colnames(results) <-  c("triplicat", ".y.", "Comparison1", "Comparison2", "n1", "n2",  "p", "p.signif", "p.adj", "p.adj.signif")
  # 
  # 
  # final_results <- results %>%
  #   left_join(cohen_results, by = c("triplicat", "Comparison1", "Comparison2")) %>%
  #   filter(Comparison2 == "same") %>%
  #   filter(p.adj < 0.05)
  # 
  # #### plot the cohen's d ####
  # t.coh <- ggplot(final_results,aes(x = triplicat, y = d_cohen.Cohens_d)) +
  #   geom_boxplot(fill = "grey", outlier.colour = "black", outlier.shape = 16, outlier.size = 2) +
  #   stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "red", position = position_dodge(width = 0.75)) +
  #   labs(title = "",
  #        x = "Comparisons",
  #        y = "Cohen's d") +
  #   theme_minimal() +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1))
  # 
  # 
  # path_to_boxplot <- paste0("outputs/beta/boxplot_cohen_T_30_10_24.pdf")
  # ggsave(filename =  path_to_boxplot, plot = t.coh, width = 5.5, height = 4.5)
  # 
  # 
  # ### nestedness ####
  # #### Testing normality of data for nestedness ####
  # 
  # step = 0.05*(max(nestedness_long$Value)-min(nestedness_long$Value))
  # brek = seq(min(nestedness_long$Value),max(nestedness_long$Value),  step)
  # 
  # n.norm = ggplot(nestedness_long, aes(x = Value)) +
  #   geom_histogram(breaks = brek, fill = "coral", color = "black", alpha = 0.7) +
  #   facet_wrap(~ triplicat, scales = "free_x") +
  #   labs(title = paste0("nest ", prefix), x = "Beta-div index", y = "Frequency") +
  #   theme_minimal()
  # 
  # path_to_boxplot <- paste0("outputs/beta/boxplot_norm_N_30_10_24.pdf")
  # ggsave(filename =  path_to_boxplot, plot = n.norm, width = 12, height = 9)
  # 
  # #### Testing mean difference between comparisons ####
  # stat_test <- nestedness_long %>%
  #   group_by(triplicat) %>%
  #   rstatix::pairwise_wilcox_test(Value ~ Comparison, p.adjust.method = "bonferroni") %>%
  #   filter(group2 == "same")  #
  # 
  # stat_test <- stat_test %>%
  #   mutate(y.position = rep(c(0.4, 0.435, 0.470, 0.505, 0.540, 0.575), 5))  
  # 
  # print(stat_test, n = 44) 
  # 
  # #### Plot the boxplot of beta div ####
  # n.beta <- ggplot(nestedness_long, aes(x = Comparison, y = Value)) +
  #   geom_boxplot(outlier.shape = NA) +  # Masquer les outliers pour ne pas interférer avec la moyenne
  #   stat_summary(fun = mean, geom = "point", shape = 18, color = "red", size = 3) +  # Ajouter la moyenne avec un losange rouge
  #   facet_wrap(~ triplicat, scales = "free_x") +  # Garder seulement l'échelle libre sur l'axe des x
  #   labs(title = "",
  #        x = "Comparison", y = "nestedness Index Value") +
  #   theme_minimal() +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #   ylim(0, 1) +  # Définir l'échelle de l'axe des y entre 0 et 1
  #   stat_pvalue_manual(
  #     stat_test, 
  #     label = "p.adj.signif", 
  #     y.position = "y.position",  # Utiliser la colonne y.position de stat_test pour placer les segments
  #     bracket.size = 0.3
  #   )
  # 
  # path_to_boxplot <- paste0("outputs/beta/boxplot_betadiv_microhab_N_30_10_24.pdf")
  # ggsave(filename =  path_to_boxplot, plot = n.beta, width = 16, height = 11)
  # 
  # #### Compute the r effsize (non parametric) ####
  # results <- nestedness_long %>%
  #   group_by(triplicat) %>%
  #   rstatix::pairwise_wilcox_test(Value ~ Comparison, p.adjust.method = "bonferroni") %>%
  #   filter(group2 == "same")  #
  # 
  # 
  # wilcox_results <- data.frame()
  # 
  # # Initialiser un dataframe pour stocker les résultats
  # r_wilcox_results <- nestedness_long %>%
  #   group_by(triplicat) %>%
  #   summarise(
  #     Comparison_Pairs = list(combn(unique(Comparison), 2, simplify = FALSE)),
  #     .groups = 'drop'
  #   )
  # 
  # 
  # for (i in 1:nrow(r_wilcox_results)) {
  #   #i = 1
  #   triplicat_name <- r_wilcox_results$triplicat[i]
  #   pairs <- r_wilcox_results$Comparison_Pairs[[i]]
  #   
  #   for (pair in pairs) {
  #     #pair = pairs[[1]]
  #     values1 <- nestedness_long$Value[nestedness_long$Comparison == pair[1] & nestedness_long$triplicat == triplicat_name]
  #     values2 <- nestedness_long$Value[nestedness_long$Comparison == pair[2] & nestedness_long$triplicat == triplicat_name]
  #     df <- as.data.frame(cbind(values1 = values1,
  #                               values2 = values2))
  #     df <- melt(df, value.name = "nest")
  #     r_value <- rstatix::wilcox_effsize(df, nest ~ variable)
  #     wilcox_results <- rbind(wilcox_results, data.frame(triplicat = triplicat_name, Comparison1 = pair[1], Comparison2 = pair[2], r = r_value))
  #     }
  #   }
  # 
  # colnames(results) <-  c("triplicat", ".y.", "Comparison1", "Comparison2", "n1", "n2",  "p", "p.signif", "p.adj", "p.adj.signif")
  # 
  # 
  # final_results <- results %>%
  #   left_join(wilcox_results, by = c("triplicat", "Comparison1", "Comparison2")) %>%
  #   filter(Comparison2 == "same") %>%
  #   filter(p.adj < 0.05)
  # 
  # #### plot the r eff size ####
  # n.r <- ggplot(final_results,aes(x = triplicat, y = r.effsize)) +
  #   geom_boxplot(fill = "grey", outlier.colour = "black", outlier.shape = 16, outlier.size = 2) +
  #   stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "red", position = position_dodge(width = 0.75)) +
  #   labs(title = "",
  #        x = "Comparisons",
  #        y = "r W-U effect size") +
  #   theme_minimal() +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1))
  # 
  # path_to_boxplot <- paste0("outputs/beta/boxplot_cohen_N_30_10_24.pdf")
  # ggsave(filename =  path_to_boxplot, plot = n.r, width = 5.5, height = 4.5)
  # 
  ### bray ####
  #### Testing normality of data for bray curtis ####
  
  step = 0.05*(max(bray_long$Value)-min(bray_long$Value))
  brek = seq(min(bray_long$Value),max(bray_long$Value),  step)
  
  b.norm = ggplot(bray_long, aes(x = Value)) +
    geom_histogram(breaks = brek, fill = "coral", color = "black", alpha = 0.7) +
    facet_wrap(~ triplicat, scales = "free_x") +
    labs(title = paste0("bray", prefix), x = "Beta-div index", y = "Frequency") +
    theme_minimal()
  
  path_to_boxplot <- paste0("outputs/beta/boxplot_norm_B_30_10_24.pdf")
  ggsave(filename =  path_to_boxplot, plot = b.norm, width = 12, height = 9)
  
  #### Testing mean difference between comparisons ####

  
  bray_long <- bray_long %>%
    mutate(imm_time = case_when(
      triplicat %in% c("CINA1", "CINA3") ~ "6m",  # Si "CINA1" ou "CINA3", assigner "6m"
      triplicat %in% c("CINA2", "CINA4") ~ "1y",  # Si "CINA2" ou "CINA4", assigner "1y"
      triplicat == "RUNA2" ~ "2y",                # Si "RUNA2", assigner "2y"
      TRUE ~ NA_character_                        # Si aucune des conditions précédentes, NA
    ))
  

  
  stat_test <- bray_long %>%
    group_by(imm_time) %>%
    rstatix::pairwise_t_test(Value ~ Comparison, p.adjust.method = "bonferroni") %>%
    filter(group2 == "same")  #
  
    
  stat_test <- stat_test %>%
    mutate(y.position = rep(c(0.8, 0.835, 0.870, 0.905, 0.940, 0.975), 3))

  print(stat_test, n = 100)
  
  #### Plot the boxplot of beta div ####
  library(forcats)
  
  
  # Modifier l'ordre des niveaux de 'imm_time' en utilisant fct_relevel
  
  bray_long$imm_time <- factor(bray_long$imm_time, levels = c("6m", "1y", "2y"))
  
  # Then, arrange the table by this column
  bray_long <- bray_long %>%
    arrange(imm_time)

  library(plot3D)
  
  
  matrice <- with(bray_long, tapply(Value, list(Comparison, imm_time), mean, na.rm = TRUE))
  
  # Calculer la somme des valeurs pour chaque ligne
  row_sums <- rowSums(matrice)
  
  # Trier les lignes par ordre croissant des sommes
  sorted_matrice <- matrice[order(row_sums), ]
  
  # Trier les lignes par ordre croissant des sommes
  sorted_matrice <- sorted_matrice[,c("2y","1y","6m")]
  
  color_palette <- rev(c("#CC66CC", "#1B9E77", "#FF7F00"))
  colvar <- matrix(c(rep(1, 7), rep(2, 7), rep(3, 7)), nrow = 7, ncol = 3)
  
  
 
  hist3D(x = 1:7, y = 1:3, z = matrice)
  

  x_labels <- c("same", "UC_UO", "DC_DO", "DO_UO", "DC_UC", "DO_UC", "DC_UO")
  y_labels <- rev(c("6m", "1y", "2y"))
  z_labels <- c("0", "0.7") 
  

  # Ajouter les labels pour les axes X et Y
  axis(1, at = 1:7, labels = x_labels, las = 2)  # Label de l'axe X
  axis(2, at = 1:3, labels = y_labels, las = 1)  # Label de l'axe Y
  library(plot3D)
  # pdf(file =  "outputs/hist3D", width = 8, height = 7)
  hist3D(x = 1:7, y = 1:3, z = sorted_matrice,
         colvar = colvar,
         xlab = "", 
         ylab = "", 
         zlab = "",
         ticks = TRUE,
         col = color_palette,
         border = "black",
         bty = "g", 
         phi = 40,  
         theta = -50,
         space=0.500,
         shade = 0.8,
         zlim = c(0, 0.7),
         cex.axis = 1e-9)
  
  # Ajouter les étiquettes pour l'axe X
  text3D(x = 1:7, y = rep(0.3, 7), z = rep(0, 7),add = TRUE, labels = x_labels, col = "black", cex = 1)
  
  # Ajouter les étiquettes pour l'axe Y
  text3D(x = rep(-0.3, 3), y = 1:3, z = rep(0, 3), add = TRUE, labels = y_labels, col = "black", cex = 1)
                
  # text3D(x = rep(0, 2), y = rep(0, 2), z = c(0, 0.7), add = TRUE, labels = z_labels, col = "black", cex = 1)
  
  
  b.beta <- ggplot(bray_long, aes(x = Comparison, y = Value, fill = imm_time)) +
    geom_boxplot(outlier.shape = NA) +  # Hide outliers to avoid interference with the mean
    stat_summary(fun = mean, geom = "point", shape = 18, color = "grey", size = 3) +  # Add red diamond for the mean
    facet_wrap(~ imm_time, ncol = 3) +  # Keep free x-scale only
    labs(title = "",
         x = "Comparison", y = "Bray Index Value") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylim(0, 1) +   # Set y-axis scale between 0 and 1
    scale_fill_manual(values = c("#CC66CC", "#1B9E77", "#FF7F00")) +  # Set colors for each value of imm_time
    stat_pvalue_manual(
      stat_test,
      label = "p.adj.signif",
      y.position = "y.position",  # Use y.position column in stat_test to place segments
      bracket.size = 0.3
    )
  
  path_to_boxplot <- paste0("outputs/beta/boxplot_betadiv_microhab_20_11_24.png")
  ggsave(filename =  path_to_boxplot, plot = b.beta, width = 13, height = 4.5)

  # install.packages("ggridges")
  # library(ggridges)
  # library(ggplot2)
  # 
  # ggplot(bray_long, aes(x = Value, y = Comparison)) +
  #   geom_density_ridges(stat="binline") +
  #   theme_ridges() + 
  #   theme(legend.position = "none")
  # 
  
  
  
  #### Compute the cohen's d ####
  results <-bray_long %>%
    group_by(triplicat) %>%
    rstatix::pairwise_t_test(Value ~ Comparison, p.adjust.method = "bonferroni") %>%
    filter(group2 == "same")  #
  
  
  cohen_results <- data.frame()
  
  # Initialiser un dataframe pour stocker les résultats
  d_cohen_results <- bray_long %>%
    group_by(triplicat) %>%
    summarise(
      Comparison_Pairs = list(combn(unique(Comparison), 2, simplify = FALSE)),
      .groups = 'drop'
    )
  
  
  for (i in 1:nrow(d_cohen_results)) {
    triplicat_name <- d_cohen_results$triplicat[i]
    pairs <- d_cohen_results$Comparison_Pairs[[i]]
    
    for (pair in pairs) {
      values1 <- bray_long$Value[bray_long$Comparison == pair[1] & bray_long$triplicat == triplicat_name]
      values2 <- bray_long$Value[bray_long$Comparison == pair[2] & bray_long$triplicat == triplicat_name]
      
      # Calculer le d de Cohen
      if (length(values1) > 1 & length(values2) > 1) {
        d_cohen_value <- effectsize::cohens_d(values1, values2)
        cohen_results <- rbind(cohen_results, data.frame(triplicat = triplicat_name, Comparison1 = pair[1], Comparison2 = pair[2], d_cohen = d_cohen_value))
      }
    }
  }
  
  colnames(results) <-  c("triplicat", ".y.", "Comparison1", "Comparison2", "n1", "n2",  "p", "p.signif", "p.adj", "p.adj.signif")
  
  
  final_results <- results %>%
    left_join(cohen_results, by = c("triplicat", "Comparison1", "Comparison2")) %>%
    filter(Comparison2 == "same") %>%
    filter(p.adj < 0.05)
  print(final_results, n = 100)
  
  #### plot the cohen's d ####
  b.coh <- ggplot(final_results,aes(x = triplicat, y = d_cohen.Cohens_d)) +
    geom_boxplot(fill = "grey", outlier.colour = "black", outlier.shape = 16, outlier.size = 2) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "red", position = position_dodge(width = 0.75)) +
    labs(title = "Boxplot des Valeurs de d de Cohen",
         x = "Comparisons",
         y = "Cohen's d") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  path_to_boxplot <- paste0("outputs/beta/boxplot_cohen_B_30_10_24.pdf")
  ggsave(filename =  path_to_boxplot, plot = b.coh, width = 5.5, height = 4.5)
  
  
  return(NULL)
}  
