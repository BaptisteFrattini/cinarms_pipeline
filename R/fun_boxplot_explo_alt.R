#' boxplot season separatively
#'
#' @param meta_data the path to the raw data file
#' @param data_full_pool the path to the raw data file (species pool)
#' @return 
#' @export
#' 

boxplot_explo_alt <- function(data_full_pool, meta_data){
  
  # data_full_pool = targets::tar_read(data_pool)
  # meta_data = targets::tar_read(metadata_data)
  
  data_pool <- read.csv(data_full_pool, header = TRUE)
  
  sort(colSums(data_pool))
  ordre <- sort(colSums(data_pool))
  data_pool <- data_pool[, names(ordre)]
  
  meta <- read.csv(meta_data[grepl("metadata", meta_data)], header = TRUE)
  
  names(data_pool)
  library(ggpubr)
  library(forcats)
  library(ggplot2)
  time <- c("6 month", "1 year", "2 years")
  col1 <- rev(c("darkolivegreen","darkolivegreen3","darkolivegreen1"))
  data_pool$imm_time <- meta$imm_time
  data_pool$immersion_season <- meta$immersion_season
  data_pool$recovery_season <- meta$recovery_season
  
  #### Selecting comparisons ####
  
  data_pool$set <- paste0(data_pool$imm_time,"_", data_pool$immersion_season)
  data_pool$set2 <- paste0(data_pool$imm_time,"_", data_pool$recovery_season)
  data_pool
  levels(as.factor(data_pool$set))
  time_imm_seas <- c("6m_imm_cold", "6m_imm_hot", "1y_imm_cold", "1y_imm_hot", "2y_imm_hot")
  levels(as.factor(data_pool$set2))
  time_rec_seas <- c("1y_rec_cold", "1y_rec_hot", "2y_rec_hot", "6m_rec_cold", "6m_rec_hot")
  #testing on bare plate 
  res.aov <- rstatix::kruskal_test(data_pool, bare_plate ~ set)
  p.sed <- rstatix::wilcox_test(data_pool, bare_plate ~ set, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  
  
  x1_bare_plate <- ggplot(data_pool, aes(x = fct_relevel(set,  "6m_imm_cold", "6m_imm_hot", "1y_imm_cold", "1y_imm_hot", "2y_imm_hot"), y = bare_plate)) +
    geom_boxplot(fill =  c("grey"),
                 color = c("dodgerblue", "firebrick3","dodgerblue", "firebrick3", "firebrick3")) +
    labs(title = "",
         x = "Immersion time",
         y = "Percentage cover of bare_plate") +
    scale_x_discrete(labels=time_imm_seas) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  res.aov <- rstatix::kruskal_test(data_pool, bare_plate ~ set2)
  p.sed <- rstatix::wilcox_test(data_pool, bare_plate ~ set2, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  
  x2_bare_plate <- ggplot(data_pool, aes(x = fct_relevel(set2,  "6m_rec_cold", "6m_rec_hot", "1y_rec_cold", "1y_rec_hot","2y_rec_hot"), y = bare_plate)) +
    geom_boxplot(fill =  c("grey"),
                 color = c("dodgerblue", "firebrick3","dodgerblue", "firebrick3", "firebrick3")) +
    labs(title = "",
         x = "Immersion time",
         y = "Percentage cover of bare_plate") +
    scale_x_discrete(labels=time_rec_seas) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  ### OU bien : 
  
  #### Immersion season
  
  data_pool_hot <- subset(data_pool, data_pool$immersion_season == "imm_hot")
  data_pool_cool <- subset(data_pool, data_pool$immersion_season == "imm_cold")
  
  data_pool_six <- subset(data_pool, data_pool$imm_time == "6m")
  data_pool_one <- subset(data_pool, data_pool$imm_time == "1y")
  
  res.aov <- rstatix::kruskal_test(data_pool_hot, bare_plate ~ imm_time)
  p.sed <- rstatix::wilcox_test(data_pool_hot, bare_plate ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y1_bare_plate <- ggplot(data_pool_hot, aes(x = fct_relevel(data_pool_hot$imm_time, "6m", "1y", "2y"), y = bare_plate)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3","darkolivegreen") ) +
    labs(title = "Immersion season : Hot",
         x = "Immersion time",
         y = "Percentage cover of bare_plate") +
    scale_x_discrete(labels=time) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  res.aov <- rstatix::kruskal_test(data_pool_cool, bare_plate ~ imm_time)
  p.sed <- rstatix::wilcox_test(data_pool_cool, bare_plate ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y2_bare_plate <- ggplot(data_pool_cool, aes(x = fct_relevel(data_pool_cool$imm_time, "6m", "1y"), y = bare_plate)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3")) +
    labs(title = "Immersion season : Cool",
         x = "Immersion time",
         y = "Percentage cover of bare_plate") +
    scale_x_discrete(labels=c("6 month", "1 year")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  res.aov <- rstatix::kruskal_test(data_pool_six, bare_plate ~ immersion_season)
  p.sed <- rstatix::wilcox_test(data_pool_six, bare_plate ~ immersion_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y3_bare_plate <- ggplot(data_pool_six, aes(x = fct_relevel(data_pool_six$immersion_season, "imm_cold", "imm_hot"), y = bare_plate)) +
    geom_boxplot(fill =  c("dodgerblue","firebrick3")) +
    labs(title = "Immersion time : 6 months",
         x = "Immersion season",
         y = "Percentage cover of bare_plate") +
    scale_x_discrete(labels=c("cool", "hot")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  res.aov <- rstatix::kruskal_test(data_pool_one, bare_plate ~ immersion_season)
  p.sed <- rstatix::wilcox_test(data_pool_one, bare_plate ~ immersion_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y4_bare_plate <- ggplot(data_pool_one, aes(x = fct_relevel(data_pool_one$immersion_season, "imm_cold", "imm_hot"), y = bare_plate)) +
    geom_boxplot(fill =  c("dodgerblue","firebrick3")) +
    labs(title = "Immersion time : one year",
         x = "Immersion season",
         y = "Percentage cover of bare_plate") +
    scale_x_discrete(labels=c("cool", "hot")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  #### Retrieval season
  
  data_pool_r_hot <- subset(data_pool, data_pool$recovery_season == "rec_hot")
  data_pool_r_cool <- subset(data_pool, data_pool$recovery_season == "rec_cold")
  
  data_pool_six <- subset(data_pool, data_pool$imm_time == "6m")
  data_pool_one <- subset(data_pool, data_pool$imm_time == "1y")
  
  res.aov <- rstatix::kruskal_test(data_pool_r_hot, bare_plate ~ imm_time)
  p.sed <- rstatix::wilcox_test(data_pool_r_hot, bare_plate ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y5_bare_plate <- ggplot(data_pool_r_hot, aes(x = fct_relevel(data_pool_r_hot$imm_time, "6m", "1y", "2y"), y = bare_plate)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3","darkolivegreen") ) +
    labs(title = "Retrieval season : Hot",
         x = "Immersion time",
         y = "Percentage cover of bare_plate") +
    scale_x_discrete(labels=time) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  res.aov <- rstatix::kruskal_test(data_pool_r_cool, bare_plate ~ imm_time)
  p.sed <- rstatix::wilcox_test(data_pool_r_cool, bare_plate ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y6_bare_plate <- ggplot(data_pool_r_cool, aes(x = fct_relevel(imm_time, "6m", "1y"), y = bare_plate)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3")) +
    labs(title = "Rerieval season : Cool",
         x = "Immersion time",
         y = "Percentage cover of bare_plate") +
    scale_x_discrete(labels=c("6 month", "1 year")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  res.aov <- rstatix::kruskal_test(data_pool_six, bare_plate ~ recovery_season)
  p.sed <- rstatix::wilcox_test(data_pool_six, bare_plate ~ recovery_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y7_bare_plate <- ggplot(data_pool_six, aes(x = fct_relevel(data_pool_six$recovery_season, "rec_cold", "rec_hot"), y = bare_plate)) +
    geom_boxplot(fill =  c("dodgerblue","firebrick3")) +
    labs(title = "Immersion time : 6 months",
         x = "retrieval season",
         y = "Percentage cover of bare_plate") +
    scale_x_discrete(labels=c("cool", "hot")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  res.aov <- rstatix::kruskal_test(data_pool_one, bare_plate ~ recovery_season)
  p.sed <- rstatix::wilcox_test(data_pool_one, bare_plate ~ recovery_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y8_bare_plate <- ggplot(data_pool_one, aes(x = fct_relevel(data_pool_one$recovery_season, "rec_cold", "rec_hot"), y = bare_plate)) +
    geom_boxplot(fill =  c("dodgerblue","firebrick3")) +
    labs(title = "Immersion time : one year",
         x = "retrieval season",
         y = "Percentage cover of bare_plate") +
    scale_x_discrete(labels=c("cool", "hot")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  
  
  first_sol <- cowplot::plot_grid(x1_bare_plate,x2_bare_plate ,  labels = c("imm_seas","rec_seas"))
  second_sol <- cowplot::plot_grid(y1_bare_plate, y2_bare_plate, y3_bare_plate, y4_bare_plate,
                                   y5_bare_plate, y6_bare_plate, y7_bare_plate, y8_bare_plate,
                                   labels = c(" "," "," "," "," "," "," "," "),
                                   ncol = 4,
                                   nrow = 2)
  
  
  path_to_boxplot_alt1 <- paste0("outputs/boxplot_pool_alt1.pdf")
  ggsave(filename =  path_to_boxplot_alt1, plot = first_sol, width = 12, height = 8)  
  
  path_to_boxplot_alt2 <- paste0("outputs/boxplot_pool_alt2.pdf")
  ggsave(filename =  path_to_boxplot_alt2, plot = second_sol, width = 12, height = 8)
  

return(c(path_to_boxplot_alt1, path_to_boxplot_alt2))

}