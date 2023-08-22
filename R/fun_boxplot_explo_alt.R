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
  data_pool$Ascidiacea <- rowSums(data.frame(data_pool$ascidiacea_c, data_pool$ascidiacea_s))
  data_pool <- data_pool[,-c(3,4)]
  
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
  
  #### set comparisons ####
  
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
    stat_pvalue_manual(p.sed)+
    annotate(geom="text", x=1, y=-0.1, label = paste0("N = ",length(data_pool$set[grepl("6m_imm_cold", data_pool$set)])),
             color="black")+
    annotate(geom="text", x=2, y=-0.1, label = paste0("N = ",length(data_pool$set[grepl("6m_imm_hot", data_pool$set)])),
             color="black")+
    annotate(geom="text", x=3, y=-0.1, label = paste0("N = ",length(data_pool$set[grepl("1y_imm_cold", data_pool$set)])),
             color="black")+
    annotate(geom="text", x=4, y=-0.1, label = paste0("N = ",length(data_pool$set[grepl("1y_imm_hot", data_pool$set)])),
             color="black")+
    annotate(geom="text", x=5, y=-0.1, label = paste0("N = ",length(data_pool$set[grepl("2y_imm_hot", data_pool$set)])),
             color="black")
    
  
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
    stat_pvalue_manual(p.sed)+
    annotate(geom="text", x=1, y=-0.1, label = paste0("N = ",length(data_pool$set2[grepl("6m_rec_cold", data_pool$set2)])),
             color="black")+
    annotate(geom="text", x=2, y=-0.1, label = paste0("N = ",length(data_pool$set2[grepl("6m_rec_hot", data_pool$set2)])),
             color="black")+
    annotate(geom="text", x=3, y=-0.1, label = paste0("N = ",length(data_pool$set2[grepl("1y_rec_cold", data_pool$set2)])),
             color="black")+
    annotate(geom="text", x=4, y=-0.1, label = paste0("N = ",length(data_pool$set2[grepl("1y_rec_hot", data_pool$set2)])),
             color="black")+
    annotate(geom="text", x=5, y=-0.1, label = paste0("N = ",length(data_pool$set2[grepl("2y_rec_hot", data_pool$set2)])),
             color="black")
  
  #### Other solution  ####
  #### Deployment season
  
  data_pool_hot <- subset(data_pool, data_pool$immersion_season == "imm_hot")
  data_pool_cool <- subset(data_pool, data_pool$immersion_season == "imm_cold")
  
  data_pool_six <- subset(data_pool, data_pool$imm_time == "6m")
  data_pool_one <- subset(data_pool, data_pool$imm_time == "1y")
  
  #### bare plate ####
  
  ggplot(data_pool, aes(x=bare_plate)) + 
    geom_density()                        #Data not normal
  ggpubr::ggqqplot(data_pool$bare_plate)  #Data not normal
  shapiro.test(data_pool$bare_plate)      #Data not normal
  res.aov2 <- aov(bare_plate ~ imm_time, data=data_pool)
  residuals <- resid(res.aov2)
  plot(data_pool$bare_plate, residuals, xlab="Measurement", ylab="Residuals") # really bad
  abline(0,0)
  ggplot() +
    geom_qq(aes(sample = rstandard(res.aov2))) +
    geom_abline(color = "red") +
    coord_fixed() # QQplot on residuals --> not OK
  
  res.aov <- rstatix::kruskal_test(data_pool_hot, bare_plate ~ imm_time)
  p.sed <- rstatix::wilcox_test(data_pool_hot, bare_plate ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  
  y1_bare_plate <- ggplot(data_pool_hot, aes(x = fct_relevel(data_pool_hot$imm_time, "6m", "1y", "2y"), y = bare_plate)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3","darkolivegreen") ) +
    labs(title = "Deployment season : Hot",
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
    labs(title = "Deployment season : Cool",
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
         x = "Deployment season",
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
         x = "Deployment season",
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
  
  
  path_to_boxplot_alt2 <- paste0("outputs/box/boxplot_bare_plate.pdf")
  ggsave(filename =  path_to_boxplot_alt2, plot = second_sol, width = 12, height = 8)
  
  
  #### Sediments ####
  #### Deployment season
  
  ggplot(data_pool, aes(x=sediment)) + 
    geom_density()                        #Data not normal
  ggpubr::ggqqplot(data_pool$sediment)  #Data not normal
  shapiro.test(data_pool$sediment)      #Data not normal
  res.aov2 <- aov(sediment ~ imm_time*immersion_season, data=data_pool)
  residuals <- resid(res.aov2)
  plot(data_pool$sediment, residuals, xlab="Measurement", ylab="Residuals") # really bad
  abline(0,0)
  ggplot() +
    geom_qq(aes(sample = rstandard(res.aov2))) +
    geom_abline(color = "red") +
    coord_fixed() # QQplot on residuals --> not OK
  
  res.aov <- rstatix::kruskal_test(data_pool_hot, sediment ~ imm_time)
  p.sed <- rstatix::wilcox_test(data_pool_hot, sediment ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  
  y1_sediment <- ggplot(data_pool_hot, aes(x = fct_relevel(data_pool_hot$imm_time, "6m", "1y", "2y"), y = sediment)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3","darkolivegreen") ) +
    labs(title = "Deployment season : Hot",
         x = "Immersion time",
         y = "Percentage cover of sediment") +
    scale_x_discrete(labels=time) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  res.aov <- rstatix::kruskal_test(data_pool_cool, sediment ~ imm_time)
  p.sed <- rstatix::wilcox_test(data_pool_cool, sediment ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y2_sediment <- ggplot(data_pool_cool, aes(x = fct_relevel(data_pool_cool$imm_time, "6m", "1y"), y = sediment)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3")) +
    labs(title = "Deployment season : Cool",
         x = "Immersion time",
         y = "Percentage cover of sediment") +
    scale_x_discrete(labels=c("6 month", "1 year")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed) 
  
  res.aov <- rstatix::kruskal_test(data_pool_six, sediment ~ immersion_season)
  p.sed <- rstatix::wilcox_test(data_pool_six, sediment ~ immersion_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y3_sediment <- ggplot(data_pool_six, aes(x = fct_relevel(data_pool_six$immersion_season, "imm_cold", "imm_hot"), y = sediment)) +
    geom_boxplot(fill =  c("dodgerblue","firebrick3")) +
    labs(title = "Immersion time : 6 months",
         x = "Deployment season",
         y = "Percentage cover of sediment") +
    scale_x_discrete(labels=c("cool", "hot")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  res.aov <- rstatix::kruskal_test(data_pool_one, sediment ~ immersion_season)
  p.sed <- rstatix::wilcox_test(data_pool_one, sediment ~ immersion_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y4_sediment <- ggplot(data_pool_one, aes(x = fct_relevel(data_pool_one$immersion_season, "imm_cold", "imm_hot"), y = sediment)) +
    geom_boxplot(fill =  c("dodgerblue","firebrick3")) +
    labs(title = "Immersion time : one year",
         x = "Deployment season",
         y = "Percentage cover of sediment") +
    scale_x_discrete(labels=c("cool", "hot")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  #### Retrieval season
  
  res.aov <- rstatix::kruskal_test(data_pool_r_hot, sediment ~ imm_time)
  p.sed <- rstatix::wilcox_test(data_pool_r_hot, sediment ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y5_sediment <- ggplot(data_pool_r_hot, aes(x = fct_relevel(data_pool_r_hot$imm_time, "6m", "1y", "2y"), y = sediment)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3","darkolivegreen") ) +
    labs(title = "Retrieval season : Hot",
         x = "Immersion time",
         y = "Percentage cover of sediment") +
    scale_x_discrete(labels=time) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  res.aov <- rstatix::kruskal_test(data_pool_r_cool, sediment ~ imm_time)
  p.sed <- rstatix::wilcox_test(data_pool_r_cool, sediment ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y6_sediment <- ggplot(data_pool_r_cool, aes(x = fct_relevel(imm_time, "6m", "1y"), y = sediment)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3")) +
    labs(title = "Rerieval season : Cool",
         x = "Immersion time",
         y = "Percentage cover of sediment") +
    scale_x_discrete(labels=c("6 month", "1 year")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)  
  
  res.aov <- rstatix::kruskal_test(data_pool_six, sediment ~ recovery_season)
  p.sed <- rstatix::wilcox_test(data_pool_six, sediment ~ recovery_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y7_sediment <- ggplot(data_pool_six, aes(x = fct_relevel(data_pool_six$recovery_season, "rec_cold", "rec_hot"), y = sediment)) +
    geom_boxplot(fill =  c("dodgerblue","firebrick3")) +
    labs(title = "Immersion time : 6 months",
         x = "retrieval season",
         y = "Percentage cover of sediment") +
    scale_x_discrete(labels=c("cool", "hot")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  res.aov <- rstatix::kruskal_test(data_pool_one, sediment ~ recovery_season)
  p.sed <- rstatix::wilcox_test(data_pool_one, sediment ~ recovery_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y8_sediment <- ggplot(data_pool_one, aes(x = fct_relevel(data_pool_one$recovery_season, "rec_cold", "rec_hot"), y = sediment)) +
    geom_boxplot(fill =  c("dodgerblue","firebrick3")) +
    labs(title = "Immersion time : one year",
         x = "retrieval season",
         y = "Percentage cover of sediment") +
    scale_x_discrete(labels=c("cool", "hot")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)  
  
  
  
  second_sol <- cowplot::plot_grid(y1_sediment, y2_sediment, y3_sediment, y4_sediment,
                                   y5_sediment, y6_sediment, y7_sediment, y8_sediment,
                                   labels = c(" "," "," "," "," "," "," "," "),
                                   ncol = 4,
                                   nrow = 2)
  
  
  path_to_boxplot_alt2 <- paste0("outputs/box/boxplot_sediment.pdf")
  ggsave(filename =  path_to_boxplot_alt2, plot = second_sol, width = 12, height = 8)
  
  #### CCA ####
  #### Deployment season
  
  ggplot(data_pool, aes(x=CCA)) + 
    geom_density()                        #Data not normal
  ggpubr::ggqqplot(data_pool$CCA)  #Data not normal
  shapiro.test(data_pool$CCA)      #Data not normal
  res.aov2 <- aov(CCA ~ imm_time*immersion_season, data=data_pool)
  residuals <- resid(res.aov2)
  plot(data_pool$CCA, residuals, xlab="Measurement", ylab="Residuals") # really bad
  abline(0,0)
  ggplot() +
    geom_qq(aes(sample = rstandard(res.aov2))) +
    geom_abline(color = "red") +
    coord_fixed() # QQplot on residuals --> not OK
  
  res.aov <- rstatix::kruskal_test(data_pool_hot, CCA ~ imm_time)
  p.sed <- rstatix::wilcox_test(data_pool_hot, CCA ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  
  y1_CCA <- ggplot(data_pool_hot, aes(x = fct_relevel(data_pool_hot$imm_time, "6m", "1y", "2y"), y = CCA)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3","darkolivegreen") ) +
    labs(title = "Deployment season : Hot",
         x = "Immersion time",
         y = "Percentage cover of CCA") +
    scale_x_discrete(labels=time) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed) 
  
  res.aov <- rstatix::kruskal_test(data_pool_cool, CCA ~ imm_time)
  p.sed <- rstatix::wilcox_test(data_pool_cool, CCA ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y2_CCA <- ggplot(data_pool_cool, aes(x = fct_relevel(data_pool_cool$imm_time, "6m", "1y"), y = CCA)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3")) +
    labs(title = "Deployment season : Cool",
         x = "Immersion time",
         y = "Percentage cover of CCA") +
    scale_x_discrete(labels=c("6 month", "1 year")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed) 
  
  res.aov <- rstatix::kruskal_test(data_pool_six, CCA ~ immersion_season)
  p.sed <- rstatix::wilcox_test(data_pool_six, CCA ~ immersion_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y3_CCA <- ggplot(data_pool_six, aes(x = fct_relevel(data_pool_six$immersion_season, "imm_cold", "imm_hot"), y = CCA)) +
    geom_boxplot(fill =  c("dodgerblue","firebrick3")) +
    labs(title = "Immersion time : 6 months",
         x = "Deployment season",
         y = "Percentage cover of CCA") +
    scale_x_discrete(labels=c("cool", "hot")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  res.aov <- rstatix::kruskal_test(data_pool_one, CCA ~ immersion_season)
  p.sed <- rstatix::wilcox_test(data_pool_one, CCA ~ immersion_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y4_CCA <- ggplot(data_pool_one, aes(x = fct_relevel(data_pool_one$immersion_season, "imm_cold", "imm_hot"), y = CCA)) +
    geom_boxplot(fill =  c("dodgerblue","firebrick3")) +
    labs(title = "Immersion time : one year",
         x = "Deployment season",
         y = "Percentage cover of CCA") +
    scale_x_discrete(labels=c("cool", "hot")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  #### Retrieval season
  
  res.aov <- rstatix::kruskal_test(data_pool_r_hot, CCA ~ imm_time)
  p.sed <- rstatix::wilcox_test(data_pool_r_hot, CCA ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y5_CCA <- ggplot(data_pool_r_hot, aes(x = fct_relevel(data_pool_r_hot$imm_time, "6m", "1y", "2y"), y = CCA)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3","darkolivegreen") ) +
    labs(title = "Retrieval season : Hot",
         x = "Immersion time",
         y = "Percentage cover of CCA") +
    scale_x_discrete(labels=time) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  res.aov <- rstatix::kruskal_test(data_pool_r_cool, CCA ~ imm_time)
  p.sed <- rstatix::wilcox_test(data_pool_r_cool, CCA ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y6_CCA <- ggplot(data_pool_r_cool, aes(x = fct_relevel(imm_time, "6m", "1y"), y = CCA)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3")) +
    labs(title = "Rerieval season : Cool",
         x = "Immersion time",
         y = "Percentage cover of CCA") +
    scale_x_discrete(labels=c("6 month", "1 year")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed) 
  
  res.aov <- rstatix::kruskal_test(data_pool_six, CCA ~ recovery_season)
  
  p.sed <- rstatix::wilcox_test(data_pool_six, CCA ~ recovery_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y7_CCA <- ggplot(data_pool_six, aes(x = fct_relevel(data_pool_six$recovery_season, "rec_cold", "rec_hot"), y = CCA)) +
    geom_boxplot(fill =  c("dodgerblue","firebrick3")) +
    labs(title = "Immersion time : 6 months",
         x = "retrieval season",
         y = "Percentage cover of CCA") +
    scale_x_discrete(labels=c("cool", "hot")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  res.aov <- rstatix::kruskal_test(data_pool_one, CCA ~ recovery_season)
  p.sed <- rstatix::wilcox_test(data_pool_one, CCA ~ recovery_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y8_CCA <- ggplot(data_pool_one, aes(x = fct_relevel(data_pool_one$recovery_season, "rec_cold", "rec_hot"), y = CCA)) +
    geom_boxplot(fill =  c("dodgerblue","firebrick3")) +
    labs(title = "Immersion time : one year",
         x = "retrieval season",
         y = "Percentage cover of CCA") +
    scale_x_discrete(labels=c("cool", "hot")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed) 
  
  
  second_sol <- cowplot::plot_grid(y1_CCA, y2_CCA, y3_CCA, y4_CCA,
                                   y5_CCA, y6_CCA, y7_CCA, y8_CCA,
                                   labels = c(" "," "," "," "," "," "," "," "),
                                   ncol = 4,
                                   nrow = 2)
  
  
  path_to_boxplot_alt2 <- paste0("outputs/box/boxplot_CCA.pdf")
  ggsave(filename =  path_to_boxplot_alt2, plot = second_sol, width = 12, height = 8)
  
  #### Porifera ####
  
  ggplot(data_pool, aes(x=porifera)) + 
    geom_density()                        #Data not normal
  ggpubr::ggqqplot(data_pool$porifera)  #Data not normal
  shapiro.test(data_pool$porifera)      #Data not normal
  res.aov2 <- aov(porifera ~ imm_time*immersion_season, data=data_pool)
  residuals <- resid(res.aov2)
  plot(data_pool$porifera, residuals, xlab="Measurement", ylab="Residuals") # really bad
  abline(0,0)
  ggplot() +
    geom_qq(aes(sample = rstandard(res.aov2))) +
    geom_abline(color = "red") +
    coord_fixed() # QQplot on residuals --> not OK
  
  res.aov <- rstatix::kruskal_test(data_pool_hot, porifera ~ imm_time)
  p.sed <- rstatix::wilcox_test(data_pool_hot, porifera ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  
  y1_porifera <- ggplot(data_pool_hot, aes(x = fct_relevel(data_pool_hot$imm_time, "6m", "1y", "2y"), y = porifera)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3","darkolivegreen") ) +
    labs(title = "Deployment season : Hot",
         x = "Immersion time",
         y = "Percentage cover of porifera") +
    scale_x_discrete(labels=time) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed) 
  
  res.aov <- rstatix::kruskal_test(data_pool_cool, porifera ~ imm_time)
  p.sed <- rstatix::wilcox_test(data_pool_cool, porifera ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y2_porifera <- ggplot(data_pool_cool, aes(x = fct_relevel(data_pool_cool$imm_time, "6m", "1y"), y = porifera)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3")) +
    labs(title = "Deployment season : Cool",
         x = "Immersion time",
         y = "Percentage cover of porifera") +
    scale_x_discrete(labels=c("6 month", "1 year")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  res.aov <- rstatix::kruskal_test(data_pool_six, porifera ~ immersion_season)
  p.sed <- rstatix::wilcox_test(data_pool_six, porifera ~ immersion_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y3_porifera <- ggplot(data_pool_six, aes(x = fct_relevel(data_pool_six$immersion_season, "imm_cold", "imm_hot"), y = porifera)) +
    geom_boxplot(fill =  c("dodgerblue","firebrick3")) +
    labs(title = "Immersion time : 6 months",
         x = "Deployment season",
         y = "Percentage cover of porifera") +
    scale_x_discrete(labels=c("cool", "hot")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  res.aov <- rstatix::kruskal_test(data_pool_one, porifera ~ immersion_season)
  p.sed <- rstatix::wilcox_test(data_pool_one, porifera ~ immersion_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y4_porifera <- ggplot(data_pool_one, aes(x = fct_relevel(data_pool_one$immersion_season, "imm_cold", "imm_hot"), y = porifera)) +
    geom_boxplot(fill =  c("dodgerblue","firebrick3")) +
    labs(title = "Immersion time : one year",
         x = "Deployment season",
         y = "Percentage cover of porifera") +
    scale_x_discrete(labels=c("cool", "hot")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  #### Retrieval season
  
  res.aov <- rstatix::kruskal_test(data_pool_r_hot, porifera ~ imm_time)
  p.sed <- rstatix::wilcox_test(data_pool_r_hot, porifera ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y5_porifera <- ggplot(data_pool_r_hot, aes(x = fct_relevel(data_pool_r_hot$imm_time, "6m", "1y", "2y"), y = porifera)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3","darkolivegreen") ) +
    labs(title = "Retrieval season : Hot",
         x = "Immersion time",
         y = "Percentage cover of porifera") +
    scale_x_discrete(labels=time) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  res.aov <- rstatix::kruskal_test(data_pool_r_cool, porifera ~ imm_time)
  p.sed <- rstatix::wilcox_test(data_pool_r_cool, porifera ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y6_porifera <- ggplot(data_pool_r_cool, aes(x = fct_relevel(imm_time, "6m", "1y"), y = porifera)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3")) +
    labs(title = "Rerieval season : Cool",
         x = "Immersion time",
         y = "Percentage cover of porifera") +
    scale_x_discrete(labels=c("6 month", "1 year")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed) 
  res.aov <- rstatix::kruskal_test(data_pool_six, porifera ~ recovery_season)
  
  p.sed <- rstatix::wilcox_test(data_pool_six, porifera ~ recovery_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y7_porifera <- ggplot(data_pool_six, aes(x = fct_relevel(data_pool_six$recovery_season, "rec_cold", "rec_hot"), y = porifera)) +
    geom_boxplot(fill =  c("dodgerblue","firebrick3")) +
    labs(title = "Immersion time : 6 months",
         x = "retrieval season",
         y = "Percentage cover of porifera") +
    scale_x_discrete(labels=c("cool", "hot")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  res.aov <- rstatix::kruskal_test(data_pool_one, porifera ~ recovery_season)
  p.sed <- rstatix::wilcox_test(data_pool_one, porifera ~ recovery_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y8_porifera <- ggplot(data_pool_one, aes(x = fct_relevel(data_pool_one$recovery_season, "rec_cold", "rec_hot"), y = porifera)) +
    geom_boxplot(fill =  c("dodgerblue","firebrick3")) +
    labs(title = "Immersion time : one year",
         x = "retrieval season",
         y = "Percentage cover of porifera") +
    scale_x_discrete(labels=c("cool", "hot")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed) 
  
  
  second_sol <- cowplot::plot_grid(y1_porifera, y2_porifera, y3_porifera, y4_porifera,
                                   y5_porifera, y6_porifera, y7_porifera, y8_porifera,
                                   labels = c(" "," "," "," "," "," "," "," "),
                                   ncol = 4,
                                   nrow = 2)
  
  
  path_to_boxplot_alt2 <- paste0("outputs/box/boxplot_porifera.pdf")
  ggsave(filename =  path_to_boxplot_alt2, plot = second_sol, width = 12, height = 8)
  
  #### Ascidiacea ####
  
  ggplot(data_pool, aes(x=Ascidiacea)) + 
    geom_density()                        #Data not normal
  ggpubr::ggqqplot(data_pool$Ascidiacea)  #Data not normal
  shapiro.test(data_pool$Ascidiacea)      #Data not normal
  res.aov2 <- aov(Ascidiacea ~ imm_time*immersion_season, data=data_pool)
  residuals <- resid(res.aov2)
  plot(data_pool$Ascidiacea, residuals, xlab="Measurement", ylab="Residuals") # really bad
  abline(0,0)
  ggplot() +
    geom_qq(aes(sample = rstandard(res.aov2))) +
    geom_abline(color = "red") +
    coord_fixed() # QQplot on residuals --> not OK
  
  res.aov <- rstatix::kruskal_test(data_pool_hot, Ascidiacea ~ imm_time)
  p.sed <- rstatix::wilcox_test(data_pool_hot, Ascidiacea ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  
  y1_Ascidiacea <- ggplot(data_pool_hot, aes(x = fct_relevel(data_pool_hot$imm_time, "6m", "1y", "2y"), y = Ascidiacea)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3","darkolivegreen") ) +
    labs(title = "Deployment season : Hot",
         x = "Immersion time",
         y = "Percentage cover of Ascidiacea") +
    scale_x_discrete(labels=time) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed) 
  
  res.aov <- rstatix::kruskal_test(data_pool_cool, Ascidiacea ~ imm_time)
  p.sed <- rstatix::wilcox_test(data_pool_cool, Ascidiacea ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y2_Ascidiacea <- ggplot(data_pool_cool, aes(x = fct_relevel(data_pool_cool$imm_time, "6m", "1y"), y = Ascidiacea)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3")) +
    labs(title = "Deployment season : Cool",
         x = "Immersion time",
         y = "Percentage cover of Ascidiacea") +
    scale_x_discrete(labels=c("6 month", "1 year")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed) 
  
  res.aov <- rstatix::kruskal_test(data_pool_six, Ascidiacea ~ immersion_season)
  p.sed <- rstatix::wilcox_test(data_pool_six, Ascidiacea ~ immersion_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y3_Ascidiacea <- ggplot(data_pool_six, aes(x = fct_relevel(data_pool_six$immersion_season, "imm_cold", "imm_hot"), y = Ascidiacea)) +
    geom_boxplot(fill =  c("dodgerblue","firebrick3")) +
    labs(title = "Immersion time : 6 months",
         x = "Deployment season",
         y = "Percentage cover of Ascidiacea") +
    scale_x_discrete(labels=c("cool", "hot")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  res.aov <- rstatix::kruskal_test(data_pool_one, Ascidiacea ~ immersion_season)
  p.sed <- rstatix::wilcox_test(data_pool_one, Ascidiacea ~ immersion_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y4_Ascidiacea <- ggplot(data_pool_one, aes(x = fct_relevel(data_pool_one$immersion_season, "imm_cold", "imm_hot"), y = Ascidiacea)) +
    geom_boxplot(fill =  c("dodgerblue","firebrick3")) +
    labs(title = "Immersion time : one year",
         x = "Deployment season",
         y = "Percentage cover of Ascidiacea") +
    scale_x_discrete(labels=c("cool", "hot")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  #### Retrieval season
  
  res.aov <- rstatix::kruskal_test(data_pool_r_hot, Ascidiacea ~ imm_time)
  p.sed <- rstatix::wilcox_test(data_pool_r_hot, Ascidiacea ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y5_Ascidiacea <- ggplot(data_pool_r_hot, aes(x = fct_relevel(data_pool_r_hot$imm_time, "6m", "1y", "2y"), y = Ascidiacea)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3","darkolivegreen") ) +
    labs(title = "Retrieval season : Hot",
         x = "Immersion time",
         y = "Percentage cover of Ascidiacea") +
    scale_x_discrete(labels=time) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  res.aov <- rstatix::kruskal_test(data_pool_r_cool, Ascidiacea ~ imm_time)
  p.sed <- rstatix::wilcox_test(data_pool_r_cool, Ascidiacea ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y6_Ascidiacea <- ggplot(data_pool_r_cool, aes(x = fct_relevel(imm_time, "6m", "1y"), y = Ascidiacea)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3")) +
    labs(title = "Rerieval season : Cool",
         x = "Immersion time",
         y = "Percentage cover of Ascidiacea") +
    scale_x_discrete(labels=c("6 month", "1 year")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed) 
  res.aov <- rstatix::kruskal_test(data_pool_six, Ascidiacea ~ recovery_season)
  
  p.sed <- rstatix::wilcox_test(data_pool_six, Ascidiacea ~ recovery_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y7_Ascidiacea <- ggplot(data_pool_six, aes(x = fct_relevel(data_pool_six$recovery_season, "rec_cold", "rec_hot"), y = Ascidiacea)) +
    geom_boxplot(fill =  c("dodgerblue","firebrick3")) +
    labs(title = "Immersion time : 6 months",
         x = "retrieval season",
         y = "Percentage cover of Ascidiacea") +
    scale_x_discrete(labels=c("cool", "hot")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  res.aov <- rstatix::kruskal_test(data_pool_one, Ascidiacea ~ recovery_season)
  p.sed <- rstatix::wilcox_test(data_pool_one, Ascidiacea ~ recovery_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y8_Ascidiacea <- ggplot(data_pool_one, aes(x = fct_relevel(data_pool_one$recovery_season, "rec_cold", "rec_hot"), y = Ascidiacea)) +
    geom_boxplot(fill =  c("dodgerblue","firebrick3")) +
    labs(title = "Immersion time : one year",
         x = "retrieval season",
         y = "Percentage cover of Ascidiacea") +
    scale_x_discrete(labels=c("cool", "hot")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed) 
  
  
  second_sol <- cowplot::plot_grid(y1_Ascidiacea, y2_Ascidiacea, y3_Ascidiacea, y4_Ascidiacea,
                                   y5_Ascidiacea, y6_Ascidiacea, y7_Ascidiacea, y8_Ascidiacea,
                                   labels = c(" "," "," "," "," "," "," "," "),
                                   ncol = 4,
                                   nrow = 2)
  
  
  path_to_boxplot_alt2 <- paste0("outputs/box/boxplot_Ascidiacea.pdf")
  ggsave(filename =  path_to_boxplot_alt2, plot = second_sol, width = 12, height = 8)
  
  #### Foraminifera ####
  
  ggplot(data_pool, aes(x=foraminifera)) + 
    geom_density()                        #Data not normal
  ggpubr::ggqqplot(data_pool$foraminifera)  #Data not normal
  shapiro.test(data_pool$foraminifera)      #Data not normal
  res.aov2 <- aov(foraminifera ~ imm_time*immersion_season, data=data_pool)
  residuals <- resid(res.aov2)
  plot(data_pool$foraminifera, residuals, xlab="Measurement", ylab="Residuals") # really bad
  abline(0,0)
  ggplot() +
    geom_qq(aes(sample = rstandard(res.aov2))) +
    geom_abline(color = "red") +
    coord_fixed() # QQplot on residuals --> not OK
  
  res.aov <- rstatix::kruskal_test(data_pool_hot, foraminifera ~ imm_time)
  p.sed <- rstatix::wilcox_test(data_pool_hot, foraminifera ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  
  y1_foraminifera <- ggplot(data_pool_hot, aes(x = fct_relevel(data_pool_hot$imm_time, "6m", "1y", "2y"), y = foraminifera)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3","darkolivegreen") ) +
    labs(title = "Deployment season : Hot",
         x = "Immersion time",
         y = "Percentage cover of foraminifera") +
    scale_x_discrete(labels=time) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  res.aov <- rstatix::kruskal_test(data_pool_cool, foraminifera ~ imm_time)
  p.sed <- rstatix::wilcox_test(data_pool_cool, foraminifera ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y2_foraminifera <- ggplot(data_pool_cool, aes(x = fct_relevel(data_pool_cool$imm_time, "6m", "1y"), y = foraminifera)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3")) +
    labs(title = "Deployment season : Cool",
         x = "Immersion time",
         y = "Percentage cover of foraminifera") +
    scale_x_discrete(labels=c("6 month", "1 year")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed) 
  
  res.aov <- rstatix::kruskal_test(data_pool_six, foraminifera ~ immersion_season)
  p.sed <- rstatix::wilcox_test(data_pool_six, foraminifera ~ immersion_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y3_foraminifera <- ggplot(data_pool_six, aes(x = fct_relevel(data_pool_six$immersion_season, "imm_cold", "imm_hot"), y = foraminifera)) +
    geom_boxplot(fill =  c("dodgerblue","firebrick3")) +
    labs(title = "Immersion time : 6 months",
         x = "Deployment season",
         y = "Percentage cover of foraminifera") +
    scale_x_discrete(labels=c("cool", "hot")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  res.aov <- rstatix::kruskal_test(data_pool_one, foraminifera ~ immersion_season)
  p.sed <- rstatix::wilcox_test(data_pool_one, foraminifera ~ immersion_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y4_foraminifera <- ggplot(data_pool_one, aes(x = fct_relevel(data_pool_one$immersion_season, "imm_cold", "imm_hot"), y = foraminifera)) +
    geom_boxplot(fill =  c("dodgerblue","firebrick3")) +
    labs(title = "Immersion time : one year",
         x = "Deployment season",
         y = "Percentage cover of foraminifera") +
    scale_x_discrete(labels=c("cool", "hot")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  #### Retrieval season
  
  res.aov <- rstatix::kruskal_test(data_pool_r_hot, foraminifera ~ imm_time)
  p.sed <- rstatix::wilcox_test(data_pool_r_hot, foraminifera ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y5_foraminifera <- ggplot(data_pool_r_hot, aes(x = fct_relevel(data_pool_r_hot$imm_time, "6m", "1y", "2y"), y = foraminifera)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3","darkolivegreen") ) +
    labs(title = "Retrieval season : Hot",
         x = "Immersion time",
         y = "Percentage cover of foraminifera") +
    scale_x_discrete(labels=time) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  res.aov <- rstatix::kruskal_test(data_pool_r_cool, foraminifera ~ imm_time)
  p.sed <- rstatix::wilcox_test(data_pool_r_cool, foraminifera ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y6_foraminifera <- ggplot(data_pool_r_cool, aes(x = fct_relevel(imm_time, "6m", "1y"), y = foraminifera)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3")) +
    labs(title = "Rerieval season : Cool",
         x = "Immersion time",
         y = "Percentage cover of foraminifera") +
    scale_x_discrete(labels=c("6 month", "1 year")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed) 
  
  res.aov <- rstatix::kruskal_test(data_pool_six, foraminifera ~ recovery_season)
  
  p.sed <- rstatix::wilcox_test(data_pool_six, foraminifera ~ recovery_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y7_foraminifera <- ggplot(data_pool_six, aes(x = fct_relevel(data_pool_six$recovery_season, "rec_cold", "rec_hot"), y = foraminifera)) +
    geom_boxplot(fill =  c("dodgerblue","firebrick3")) +
    labs(title = "Immersion time : 6 months",
         x = "retrieval season",
         y = "Percentage cover of foraminifera") +
    scale_x_discrete(labels=c("cool", "hot")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed) 
  
  res.aov <- rstatix::kruskal_test(data_pool_one, foraminifera ~ recovery_season)
  p.sed <- rstatix::wilcox_test(data_pool_one, foraminifera ~ recovery_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y8_foraminifera <- ggplot(data_pool_one, aes(x = fct_relevel(data_pool_one$recovery_season, "rec_cold", "rec_hot"), y = foraminifera)) +
    geom_boxplot(fill =  c("dodgerblue","firebrick3")) +
    labs(title = "Immersion time : one year",
         x = "retrieval season",
         y = "Percentage cover of foraminifera") +
    scale_x_discrete(labels=c("cool", "hot")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed) 
  
  
  
  second_sol <- cowplot::plot_grid(y1_foraminifera, y2_foraminifera, y3_foraminifera, y4_foraminifera,
                                   y5_foraminifera, y6_foraminifera, y7_foraminifera, y8_foraminifera,
                                   labels = c(" "," "," "," "," "," "," "," "),
                                   ncol = 4,
                                   nrow = 2)
  
  
  path_to_boxplot_alt2 <- paste0("outputs/box/boxplot_foraminifera.pdf")
  ggsave(filename =  path_to_boxplot_alt2, plot = second_sol, width = 12, height = 8)
  
  #### Hydrozoa ####
  ggplot(data_pool, aes(x=Hydrozoa)) + 
    geom_density()                        #Data not normal
  ggpubr::ggqqplot(data_pool$Hydrozoa)  #Data not normal
  shapiro.test(data_pool$Hydrozoa)      #Data not normal
  res.aov2 <- aov(Hydrozoa ~ imm_time*immersion_season, data=data_pool)
  residuals <- resid(res.aov2)
  plot(data_pool$Hydrozoa, residuals, xlab="Measurement", ylab="Residuals") # really bad
  abline(0,0)
  ggplot() +
    geom_qq(aes(sample = rstandard(res.aov2))) +
    geom_abline(color = "red") +
    coord_fixed() # QQplot on residuals --> not OK
  
  res.aov <- rstatix::kruskal_test(data_pool_hot, Hydrozoa ~ imm_time)
  p.sed <- rstatix::wilcox_test(data_pool_hot, Hydrozoa ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  
  y1_Hydrozoa <- ggplot(data_pool_hot, aes(x = fct_relevel(data_pool_hot$imm_time, "6m", "1y", "2y"), y = Hydrozoa)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3","darkolivegreen") ) +
    labs(title = "Deployment season : Hot",
         x = "Immersion time",
         y = "Percentage cover of Hydrozoa") +
    scale_x_discrete(labels=time) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed) 
  
  res.aov <- rstatix::kruskal_test(data_pool_cool, Hydrozoa ~ imm_time)
  p.sed <- rstatix::wilcox_test(data_pool_cool, Hydrozoa ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y2_Hydrozoa <- ggplot(data_pool_cool, aes(x = fct_relevel(data_pool_cool$imm_time, "6m", "1y"), y = Hydrozoa)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3")) +
    labs(title = "Deployment season : Cool",
         x = "Immersion time",
         y = "Percentage cover of Hydrozoa") +
    scale_x_discrete(labels=c("6 month", "1 year")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed) 
  
  res.aov <- rstatix::kruskal_test(data_pool_six, Hydrozoa ~ immersion_season)
  p.sed <- rstatix::wilcox_test(data_pool_six, Hydrozoa ~ immersion_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y3_Hydrozoa <- ggplot(data_pool_six, aes(x = fct_relevel(data_pool_six$immersion_season, "imm_cold", "imm_hot"), y = Hydrozoa)) +
    geom_boxplot(fill =  c("dodgerblue","firebrick3")) +
    labs(title = "Immersion time : 6 months",
         x = "Deployment season",
         y = "Percentage cover of Hydrozoa") +
    scale_x_discrete(labels=c("cool", "hot")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  res.aov <- rstatix::kruskal_test(data_pool_one, Hydrozoa ~ immersion_season)
  p.sed <- rstatix::wilcox_test(data_pool_one, Hydrozoa ~ immersion_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y4_Hydrozoa <- ggplot(data_pool_one, aes(x = fct_relevel(data_pool_one$immersion_season, "imm_cold", "imm_hot"), y = Hydrozoa)) +
    geom_boxplot(fill =  c("dodgerblue","firebrick3")) +
    labs(title = "Immersion time : one year",
         x = "Deployment season",
         y = "Percentage cover of Hydrozoa") +
    scale_x_discrete(labels=c("cool", "hot")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  #### Retrieval season
  
  res.aov <- rstatix::kruskal_test(data_pool_r_hot, Hydrozoa ~ imm_time)
  p.sed <- rstatix::wilcox_test(data_pool_r_hot, Hydrozoa ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y5_Hydrozoa <- ggplot(data_pool_r_hot, aes(x = fct_relevel(data_pool_r_hot$imm_time, "6m", "1y", "2y"), y = Hydrozoa)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3","darkolivegreen") ) +
    labs(title = "Retrieval season : Hot",
         x = "Immersion time",
         y = "Percentage cover of Hydrozoa") +
    scale_x_discrete(labels=time) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  res.aov <- rstatix::kruskal_test(data_pool_r_cool, Hydrozoa ~ imm_time)
  p.sed <- rstatix::wilcox_test(data_pool_r_cool, Hydrozoa ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y6_Hydrozoa <- ggplot(data_pool_r_cool, aes(x = fct_relevel(imm_time, "6m", "1y"), y = Hydrozoa)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3")) +
    labs(title = "Rerieval season : Cool",
         x = "Immersion time",
         y = "Percentage cover of Hydrozoa") +
    scale_x_discrete(labels=c("6 month", "1 year")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)  
  
  res.aov <- rstatix::kruskal_test(data_pool_six, Hydrozoa ~ recovery_season)
  
  p.sed <- rstatix::wilcox_test(data_pool_six, Hydrozoa ~ recovery_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y7_Hydrozoa <- ggplot(data_pool_six, aes(x = fct_relevel(data_pool_six$recovery_season, "rec_cold", "rec_hot"), y = Hydrozoa)) +
    geom_boxplot(fill =  c("dodgerblue","firebrick3")) +
    labs(title = "Immersion time : 6 months",
         x = "retrieval season",
         y = "Percentage cover of Hydrozoa") +
    scale_x_discrete(labels=c("cool", "hot")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  res.aov <- rstatix::kruskal_test(data_pool_one, Hydrozoa ~ recovery_season)
  p.sed <- rstatix::wilcox_test(data_pool_one, Hydrozoa ~ recovery_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y8_Hydrozoa <- ggplot(data_pool_one, aes(x = fct_relevel(data_pool_one$recovery_season, "rec_cold", "rec_hot"), y = Hydrozoa)) +
    geom_boxplot(fill =  c("dodgerblue","firebrick3")) +
    labs(title = "Immersion time : one year",
         x = "retrieval season",
         y = "Percentage cover of Hydrozoa") +
    scale_x_discrete(labels=c("cool", "hot")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed) 
  
  
  
  second_sol <- cowplot::plot_grid(y1_Hydrozoa, y2_Hydrozoa, y3_Hydrozoa, y4_Hydrozoa,
                                   y5_Hydrozoa, y6_Hydrozoa, y7_Hydrozoa, y8_Hydrozoa,
                                   labels = c(" "," "," "," "," "," "," "," "),
                                   ncol = 4,
                                   nrow = 2)
  
  
  path_to_boxplot_alt2 <- paste0("outputs/box/boxplot_Hydrozoa.pdf")
  ggsave(filename =  path_to_boxplot_alt2, plot = second_sol, width = 12, height = 8)
  
  #### Annelida ####
  
  ggplot(data_pool, aes(x=annelida)) + 
    geom_density()                        #Data not normal
  ggpubr::ggqqplot(data_pool$annelida)  #Data not normal
  shapiro.test(data_pool$annelida)      #Data not normal
  res.aov2 <- aov(annelida ~ imm_time*immersion_season, data=data_pool)
  residuals <- resid(res.aov2)
  plot(data_pool$annelida, residuals, xlab="Measurement", ylab="Residuals") # really bad
  abline(0,0)
  ggplot() +
    geom_qq(aes(sample = rstandard(res.aov2))) +
    geom_abline(color = "red") +
    coord_fixed() # QQplot on residuals --> not OK
  
  res.aov <- rstatix::kruskal_test(data_pool_hot, annelida ~ imm_time)
  p.sed <- rstatix::wilcox_test(data_pool_hot, annelida ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  
  y1_annelida <- ggplot(data_pool_hot, aes(x = fct_relevel(data_pool_hot$imm_time, "6m", "1y", "2y"), y = annelida)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3","darkolivegreen") ) +
    labs(title = "Deployment season : Hot",
         x = "Immersion time",
         y = "Percentage cover of annelida") +
    scale_x_discrete(labels=time) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed) 
  
  res.aov <- rstatix::kruskal_test(data_pool_cool, annelida ~ imm_time)
  p.sed <- rstatix::wilcox_test(data_pool_cool, annelida ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y2_annelida <- ggplot(data_pool_cool, aes(x = fct_relevel(data_pool_cool$imm_time, "6m", "1y"), y = annelida)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3")) +
    labs(title = "Deployment season : Cool",
         x = "Immersion time",
         y = "Percentage cover of annelida") +
    scale_x_discrete(labels=c("6 month", "1 year")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed) 
  
  res.aov <- rstatix::kruskal_test(data_pool_six, annelida ~ immersion_season)
  p.sed <- rstatix::wilcox_test(data_pool_six, annelida ~ immersion_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y3_annelida <- ggplot(data_pool_six, aes(x = fct_relevel(data_pool_six$immersion_season, "imm_cold", "imm_hot"), y = annelida)) +
    geom_boxplot(fill =  c("dodgerblue","firebrick3")) +
    labs(title = "Immersion time : 6 months",
         x = "Deployment season",
         y = "Percentage cover of annelida") +
    scale_x_discrete(labels=c("cool", "hot")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  
  
  res.aov <- rstatix::kruskal_test(data_pool_one, annelida ~ immersion_season)
  p.sed <- rstatix::wilcox_test(data_pool_one, annelida ~ immersion_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y4_annelida <- ggplot(data_pool_one, aes(x = fct_relevel(data_pool_one$immersion_season, "imm_cold", "imm_hot"), y = annelida)) +
    geom_boxplot(fill =  c("dodgerblue","firebrick3")) +
    labs(title = "Immersion time : one year",
         x = "Deployment season",
         y = "Percentage cover of annelida") +
    scale_x_discrete(labels=c("cool", "hot")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  #### Retrieval season
  
  res.aov <- rstatix::kruskal_test(data_pool_r_hot, annelida ~ imm_time)
  p.sed <- rstatix::wilcox_test(data_pool_r_hot, annelida ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y5_annelida <- ggplot(data_pool_r_hot, aes(x = fct_relevel(data_pool_r_hot$imm_time, "6m", "1y", "2y"), y = annelida)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3","darkolivegreen") ) +
    labs(title = "Retrieval season : Hot",
         x = "Immersion time",
         y = "Percentage cover of annelida") +
    scale_x_discrete(labels=time) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  res.aov <- rstatix::kruskal_test(data_pool_r_cool, annelida ~ imm_time)
  p.sed <- rstatix::wilcox_test(data_pool_r_cool, annelida ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y6_annelida <- ggplot(data_pool_r_cool, aes(x = fct_relevel(imm_time, "6m", "1y"), y = annelida)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3")) +
    labs(title = "Rerieval season : Cool",
         x = "Immersion time",
         y = "Percentage cover of annelida") +
    scale_x_discrete(labels=c("6 month", "1 year")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed) 
  
  res.aov <- rstatix::kruskal_test(data_pool_six, annelida ~ recovery_season)
  
  p.sed <- rstatix::wilcox_test(data_pool_six, annelida ~ recovery_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y7_annelida <- ggplot(data_pool_six, aes(x = fct_relevel(data_pool_six$recovery_season, "rec_cold", "rec_hot"), y = annelida)) +
    geom_boxplot(fill =  c("dodgerblue","firebrick3")) +
    labs(title = "Immersion time : 6 months",
         x = "retrieval season",
         y = "Percentage cover of annelida") +
    scale_x_discrete(labels=c("cool", "hot")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  res.aov <- rstatix::kruskal_test(data_pool_one, annelida ~ recovery_season)
  p.sed <- rstatix::wilcox_test(data_pool_one, annelida ~ recovery_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y8_annelida <- ggplot(data_pool_one, aes(x = fct_relevel(data_pool_one$recovery_season, "rec_cold", "rec_hot"), y = annelida)) +
    geom_boxplot(fill =  c("dodgerblue","firebrick3")) +
    labs(title = "Immersion time : one year",
         x = "retrieval season",
         y = "Percentage cover of annelida") +
    scale_x_discrete(labels=c("cool", "hot")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed) 
  
  
  
  second_sol <- cowplot::plot_grid(y1_annelida, y2_annelida, y3_annelida, y4_annelida,
                                   y5_annelida, y6_annelida, y7_annelida, y8_annelida,
                                   labels = c(" "," "," "," "," "," "," "," "),
                                   ncol = 4,
                                   nrow = 2)
  
  
  path_to_boxplot_alt2 <- paste0("outputs/box/boxplot_annelida.pdf")
  ggsave(filename =  path_to_boxplot_alt2, plot = second_sol, width = 12, height = 8)
  
  #### Bryozoa ####
  ggplot(data_pool, aes(x=bryozoa)) + 
    geom_density()                        #Data not normal
  ggpubr::ggqqplot(data_pool$bryozoa)  #Data not normal
  shapiro.test(data_pool$bryozoa)      #Data not normal
  res.aov2 <- aov(bryozoa ~ imm_time*immersion_season, data=data_pool)
  residuals <- resid(res.aov2)
  plot(data_pool$bryozoa, residuals, xlab="Measurement", ylab="Residuals") # really bad
  abline(0,0)
  ggplot() +
    geom_qq(aes(sample = rstandard(res.aov2))) +
    geom_abline(color = "red") +
    coord_fixed() # QQplot on residuals --> not OK
  
  res.aov <- rstatix::kruskal_test(data_pool_hot, bryozoa ~ imm_time)
  p.sed <- rstatix::wilcox_test(data_pool_hot, bryozoa ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  
  y1_bryozoa <- ggplot(data_pool_hot, aes(x = fct_relevel(data_pool_hot$imm_time, "6m", "1y", "2y"), y = bryozoa)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3","darkolivegreen") ) +
    labs(title = "Deployment season : Hot",
         x = "Immersion time",
         y = "Percentage cover of bryozoa") +
    scale_x_discrete(labels=time) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed) 
  
  res.aov <- rstatix::kruskal_test(data_pool_cool, bryozoa ~ imm_time)
  p.sed <- rstatix::wilcox_test(data_pool_cool, bryozoa ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y2_bryozoa <- ggplot(data_pool_cool, aes(x = fct_relevel(data_pool_cool$imm_time, "6m", "1y"), y = bryozoa)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3")) +
    labs(title = "Deployment season : Cool",
         x = "Immersion time",
         y = "Percentage cover of bryozoa") +
    scale_x_discrete(labels=c("6 month", "1 year")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)  
  
  res.aov <- rstatix::kruskal_test(data_pool_six, bryozoa ~ immersion_season)
  p.sed <- rstatix::wilcox_test(data_pool_six, bryozoa ~ immersion_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y3_bryozoa <- ggplot(data_pool_six, aes(x = fct_relevel(data_pool_six$immersion_season, "imm_cold", "imm_hot"), y = bryozoa)) +
    geom_boxplot(fill =  c("dodgerblue","firebrick3")) +
    labs(title = "Immersion time : 6 months",
         x = "Deployment season",
         y = "Percentage cover of bryozoa") +
    scale_x_discrete(labels=c("cool", "hot")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed) 
  
  res.aov <- rstatix::kruskal_test(data_pool_one, bryozoa ~ immersion_season)
  p.sed <- rstatix::wilcox_test(data_pool_one, bryozoa ~ immersion_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y4_bryozoa <- ggplot(data_pool_one, aes(x = fct_relevel(data_pool_one$immersion_season, "imm_cold", "imm_hot"), y = bryozoa)) +
    geom_boxplot(fill =  c("dodgerblue","firebrick3")) +
    labs(title = "Immersion time : one year",
         x = "Deployment season",
         y = "Percentage cover of bryozoa") +
    scale_x_discrete(labels=c("cool", "hot")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  #### Retrieval season
  
  res.aov <- rstatix::kruskal_test(data_pool_r_hot, bryozoa ~ imm_time)
  p.sed <- rstatix::wilcox_test(data_pool_r_hot, bryozoa ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y5_bryozoa <- ggplot(data_pool_r_hot, aes(x = fct_relevel(data_pool_r_hot$imm_time, "6m", "1y", "2y"), y = bryozoa)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3","darkolivegreen") ) +
    labs(title = "Retrieval season : Hot",
         x = "Immersion time",
         y = "Percentage cover of bryozoa") +
    scale_x_discrete(labels=time) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  res.aov <- rstatix::kruskal_test(data_pool_r_cool, bryozoa ~ imm_time)
  p.sed <- rstatix::wilcox_test(data_pool_r_cool, bryozoa ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y6_bryozoa <- ggplot(data_pool_r_cool, aes(x = fct_relevel(imm_time, "6m", "1y"), y = bryozoa)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3")) +
    labs(title = "Rerieval season : Cool",
         x = "Immersion time",
         y = "Percentage cover of bryozoa") +
    scale_x_discrete(labels=c("6 month", "1 year")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)  
  
  res.aov <- rstatix::kruskal_test(data_pool_six, bryozoa ~ recovery_season)
  
  p.sed <- rstatix::wilcox_test(data_pool_six, bryozoa ~ recovery_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y7_bryozoa <- ggplot(data_pool_six, aes(x = fct_relevel(data_pool_six$recovery_season, "rec_cold", "rec_hot"), y = bryozoa)) +
    geom_boxplot(fill =  c("dodgerblue","firebrick3")) +
    labs(title = "Immersion time : 6 months",
         x = "retrieval season",
         y = "Percentage cover of bryozoa") +
    scale_x_discrete(labels=c("cool", "hot")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed) 
  
  res.aov <- rstatix::kruskal_test(data_pool_one, bryozoa ~ recovery_season)
  p.sed <- rstatix::wilcox_test(data_pool_one, bryozoa ~ recovery_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y8_bryozoa <- ggplot(data_pool_one, aes(x = fct_relevel(data_pool_one$recovery_season, "rec_cold", "rec_hot"), y = bryozoa)) +
    geom_boxplot(fill =  c("dodgerblue","firebrick3")) +
    labs(title = "Immersion time : one year",
         x = "retrieval season",
         y = "Percentage cover of bryozoa") +
    scale_x_discrete(labels=c("cool", "hot")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed) 
  
  
  
  second_sol <- cowplot::plot_grid(y1_bryozoa, y2_bryozoa, y3_bryozoa, y4_bryozoa,
                                   y5_bryozoa, y6_bryozoa, y7_bryozoa, y8_bryozoa,
                                   labels = c(" "," "," "," "," "," "," "," "),
                                   ncol = 4,
                                   nrow = 2)
  
  
  path_to_boxplot_alt2 <- paste0("outputs/box/boxplot_bryozoa.pdf")
  ggsave(filename =  path_to_boxplot_alt2, plot = second_sol, width = 12, height = 8)
  
  #### Bivalvia ####
  
  ggplot(data_pool, aes(x=Bivalvia)) + 
    geom_density()                        #Data not normal
  ggpubr::ggqqplot(data_pool$Bivalvia)  #Data not normal
  shapiro.test(data_pool$Bivalvia)      #Data not normal
  res.aov2 <- aov(Bivalvia ~ imm_time*immersion_season, data=data_pool)
  residuals <- resid(res.aov2)
  plot(data_pool$Bivalvia, residuals, xlab="Measurement", ylab="Residuals") # really bad
  abline(0,0)
  ggplot() +
    geom_qq(aes(sample = rstandard(res.aov2))) +
    geom_abline(color = "red") +
    coord_fixed() # QQplot on residuals --> not OK
  
  res.aov <- rstatix::kruskal_test(data_pool_hot, Bivalvia ~ imm_time)
  p.sed <- rstatix::wilcox_test(data_pool_hot, Bivalvia ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  
  y1_Bivalvia <- ggplot(data_pool_hot, aes(x = fct_relevel(data_pool_hot$imm_time, "6m", "1y", "2y"), y = Bivalvia)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3","darkolivegreen") ) +
    labs(title = "Deployment season : Hot",
         x = "Immersion time",
         y = "Percentage cover of Bivalvia") +
    scale_x_discrete(labels=time) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed) 
  
  res.aov <- rstatix::kruskal_test(data_pool_cool, Bivalvia ~ imm_time)
  p.sed <- rstatix::wilcox_test(data_pool_cool, Bivalvia ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y2_Bivalvia <- ggplot(data_pool_cool, aes(x = fct_relevel(data_pool_cool$imm_time, "6m", "1y"), y = Bivalvia)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3")) +
    labs(title = "Deployment season : Cool",
         x = "Immersion time",
         y = "Percentage cover of Bivalvia") +
    scale_x_discrete(labels=c("6 month", "1 year")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)  
  
  res.aov <- rstatix::kruskal_test(data_pool_six, Bivalvia ~ immersion_season)
  p.sed <- rstatix::wilcox_test(data_pool_six, Bivalvia ~ immersion_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y3_Bivalvia <- ggplot(data_pool_six, aes(x = fct_relevel(data_pool_six$immersion_season, "imm_cold", "imm_hot"), y = Bivalvia)) +
    geom_boxplot(fill =  c("dodgerblue","firebrick3")) +
    labs(title = "Immersion time : 6 months",
         x = "Deployment season",
         y = "Percentage cover of Bivalvia") +
    scale_x_discrete(labels=c("cool", "hot")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  res.aov <- rstatix::kruskal_test(data_pool_one, Bivalvia ~ immersion_season)
  p.sed <- rstatix::wilcox_test(data_pool_one, Bivalvia ~ immersion_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y4_Bivalvia <- ggplot(data_pool_one, aes(x = fct_relevel(data_pool_one$immersion_season, "imm_cold", "imm_hot"), y = Bivalvia)) +
    geom_boxplot(fill =  c("dodgerblue","firebrick3")) +
    labs(title = "Immersion time : one year",
         x = "Deployment season",
         y = "Percentage cover of Bivalvia") +
    scale_x_discrete(labels=c("cool", "hot")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed) 
  
  #### Retrieval season
  
  res.aov <- rstatix::kruskal_test(data_pool_r_hot, Bivalvia ~ imm_time)
  p.sed <- rstatix::wilcox_test(data_pool_r_hot, Bivalvia ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y5_Bivalvia <- ggplot(data_pool_r_hot, aes(x = fct_relevel(data_pool_r_hot$imm_time, "6m", "1y", "2y"), y = Bivalvia)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3","darkolivegreen") ) +
    labs(title = "Retrieval season : Hot",
         x = "Immersion time",
         y = "Percentage cover of Bivalvia") +
    scale_x_discrete(labels=time) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  res.aov <- rstatix::kruskal_test(data_pool_r_cool, Bivalvia ~ imm_time)
  p.sed <- rstatix::wilcox_test(data_pool_r_cool, Bivalvia ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y6_Bivalvia <- ggplot(data_pool_r_cool, aes(x = fct_relevel(imm_time, "6m", "1y"), y = Bivalvia)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3")) +
    labs(title = "Rerieval season : Cool",
         x = "Immersion time",
         y = "Percentage cover of Bivalvia") +
    scale_x_discrete(labels=c("6 month", "1 year")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)  
  
  res.aov <- rstatix::kruskal_test(data_pool_six, Bivalvia ~ recovery_season)
  
  p.sed <- rstatix::wilcox_test(data_pool_six, Bivalvia ~ recovery_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y7_Bivalvia <- ggplot(data_pool_six, aes(x = fct_relevel(data_pool_six$recovery_season, "rec_cold", "rec_hot"), y = Bivalvia)) +
    geom_boxplot(fill =  c("dodgerblue","firebrick3")) +
    labs(title = "Immersion time : 6 months",
         x = "retrieval season",
         y = "Percentage cover of Bivalvia") +
    scale_x_discrete(labels=c("cool", "hot")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed) 
  
  res.aov <- rstatix::kruskal_test(data_pool_one, Bivalvia ~ recovery_season)
  p.sed <- rstatix::wilcox_test(data_pool_one, Bivalvia ~ recovery_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y8_Bivalvia <- ggplot(data_pool_one, aes(x = fct_relevel(data_pool_one$recovery_season, "rec_cold", "rec_hot"), y = Bivalvia)) +
    geom_boxplot(fill =  c("dodgerblue","firebrick3")) +
    labs(title = "Immersion time : one year",
         x = "retrieval season",
         y = "Percentage cover of Bivalvia") +
    scale_x_discrete(labels=c("cool", "hot")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed) 
  
  
  
  second_sol <- cowplot::plot_grid(y1_Bivalvia, y2_Bivalvia, y3_Bivalvia, y4_Bivalvia,
                                   y5_Bivalvia, y6_Bivalvia, y7_Bivalvia, y8_Bivalvia,
                                   labels = c(" "," "," "," "," "," "," "," "),
                                   ncol = 4,
                                   nrow = 2)
  
  
  path_to_boxplot_alt2 <- paste0("outputs/box/boxplot_Bivalvia.pdf")
  ggsave(filename =  path_to_boxplot_alt2, plot = second_sol, width = 12, height = 8)
  
  #### Other algae ####
  
  ggplot(data_pool, aes(x=other_algae)) + 
    geom_density()                        #Data not normal
  ggpubr::ggqqplot(data_pool$other_algae)  #Data not normal
  shapiro.test(data_pool$other_algae)      #Data not normal
  res.aov2 <- aov(other_algae ~ imm_time*immersion_season, data=data_pool)
  residuals <- resid(res.aov2)
  plot(data_pool$other_algae, residuals, xlab="Measurement", ylab="Residuals") # really bad
  abline(0,0)
  ggplot() +
    geom_qq(aes(sample = rstandard(res.aov2))) +
    geom_abline(color = "red") +
    coord_fixed() # QQplot on residuals --> not OK
  
  res.aov <- rstatix::kruskal_test(data_pool_hot, other_algae ~ imm_time)
  p.sed <- rstatix::wilcox_test(data_pool_hot, other_algae ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  
  y1_other_algae <- ggplot(data_pool_hot, aes(x = fct_relevel(data_pool_hot$imm_time, "6m", "1y", "2y"), y = other_algae)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3","darkolivegreen") ) +
    labs(title = "Deployment season : Hot",
         x = "Immersion time",
         y = "Percentage cover of other algae") +
    scale_x_discrete(labels=time) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed) 
  
  res.aov <- rstatix::kruskal_test(data_pool_cool, other_algae ~ imm_time)
  p.sed <- rstatix::wilcox_test(data_pool_cool, other_algae ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y2_other_algae <- ggplot(data_pool_cool, aes(x = fct_relevel(data_pool_cool$imm_time, "6m", "1y"), y = other_algae)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3")) +
    labs(title = "Deployment season : Cool",
         x = "Immersion time",
         y = "Percentage cover of other algae") +
    scale_x_discrete(labels=c("6 month", "1 year")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)  
  
  res.aov <- rstatix::kruskal_test(data_pool_six, other_algae ~ immersion_season)
  p.sed <- rstatix::wilcox_test(data_pool_six, other_algae ~ immersion_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y3_other_algae <- ggplot(data_pool_six, aes(x = fct_relevel(data_pool_six$immersion_season, "imm_cold", "imm_hot"), y = other_algae)) +
    geom_boxplot(fill =  c("dodgerblue","firebrick3")) +
    labs(title = "Immersion time : 6 months",
         x = "Deployment season",
         y = "Percentage cover of other algae") +
    scale_x_discrete(labels=c("cool", "hot")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  res.aov <- rstatix::kruskal_test(data_pool_one, other_algae ~ immersion_season)
  p.sed <- rstatix::wilcox_test(data_pool_one, other_algae ~ immersion_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y4_other_algae <- ggplot(data_pool_one, aes(x = fct_relevel(data_pool_one$immersion_season, "imm_cold", "imm_hot"), y = other_algae)) +
    geom_boxplot(fill =  c("dodgerblue","firebrick3")) +
    labs(title = "Immersion time : one year",
         x = "Deployment season",
         y = "Percentage cover of other algae") +
    scale_x_discrete(labels=c("cool", "hot")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed) 
  
  #### Retrieval season
  
  res.aov <- rstatix::kruskal_test(data_pool_r_hot, other_algae ~ imm_time)
  p.sed <- rstatix::wilcox_test(data_pool_r_hot, other_algae ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y5_other_algae <- ggplot(data_pool_r_hot, aes(x = fct_relevel(data_pool_r_hot$imm_time, "6m", "1y", "2y"), y = other_algae)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3","darkolivegreen") ) +
    labs(title = "Retrieval season : Hot",
         x = "Immersion time",
         y = "Percentage cover of other algae") +
    scale_x_discrete(labels=time) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  res.aov <- rstatix::kruskal_test(data_pool_r_cool, other_algae ~ imm_time)
  p.sed <- rstatix::wilcox_test(data_pool_r_cool, other_algae ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y6_other_algae <- ggplot(data_pool_r_cool, aes(x = fct_relevel(imm_time, "6m", "1y"), y = other_algae)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3")) +
    labs(title = "Rerieval season : Cool",
         x = "Immersion time",
         y = "Percentage cover of other algae") +
    scale_x_discrete(labels=c("6 month", "1 year")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed) 
  
  res.aov <- rstatix::kruskal_test(data_pool_six, other_algae ~ recovery_season)
  
  p.sed <- rstatix::wilcox_test(data_pool_six, other_algae ~ recovery_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y7_other_algae <- ggplot(data_pool_six, aes(x = fct_relevel(data_pool_six$recovery_season, "rec_cold", "rec_hot"), y = other_algae)) +
    geom_boxplot(fill =  c("dodgerblue","firebrick3")) +
    labs(title = "Immersion time : 6 months",
         x = "retrieval season",
         y = "Percentage cover of other algae") +
    scale_x_discrete(labels=c("cool", "hot")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed) 
  
  res.aov <- rstatix::kruskal_test(data_pool_one, other_algae ~ recovery_season)
  p.sed <- rstatix::wilcox_test(data_pool_one, other_algae ~ recovery_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y8_other_algae <- ggplot(data_pool_one, aes(x = fct_relevel(data_pool_one$recovery_season, "rec_cold", "rec_hot"), y = other_algae)) +
    geom_boxplot(fill =  c("dodgerblue","firebrick3")) +
    labs(title = "Immersion time : one year",
         x = "retrieval season",
         y = "Percentage cover of other algae") +
    scale_x_discrete(labels=c("cool", "hot")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed) 
  
  
  
  second_sol <- cowplot::plot_grid(y1_other_algae, y2_other_algae, y3_other_algae, y4_other_algae,
                                   y5_other_algae, y6_other_algae, y7_other_algae, y8_other_algae,
                                   labels = c(" "," "," "," "," "," "," "," "),
                                   ncol = 4,
                                   nrow = 2)
  
  
  path_to_boxplot_alt2 <- paste0("outputs/box/boxplot_other_algae.pdf")
  ggsave(filename =  path_to_boxplot_alt2, plot = second_sol, width = 12, height = 8)
  
  
  #### Prokariotic biotas ####
  ggplot(data_pool, aes(x=prokariotic_biotas)) + 
    geom_density()                        #Data not normal
  ggpubr::ggqqplot(data_pool$prokariotic_biotas)  #Data not normal
  shapiro.test(data_pool$prokariotic_biotas)      #Data not normal
  res.aov2 <- aov(prokariotic_biotas ~ imm_time*immersion_season, data=data_pool)
  residuals <- resid(res.aov2)
  plot(data_pool$prokariotic_biotas, residuals, xlab="Measurement", ylab="Residuals") # really bad
  abline(0,0)
  ggplot() +
    geom_qq(aes(sample = rstandard(res.aov2))) +
    geom_abline(color = "red") +
    coord_fixed() # QQplot on residuals --> not OK
  
  res.aov <- rstatix::kruskal_test(data_pool_hot, prokariotic_biotas ~ imm_time)
  p.sed <- rstatix::wilcox_test(data_pool_hot, prokariotic_biotas ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  
  y1_prokariotic_biotas <- ggplot(data_pool_hot, aes(x = fct_relevel(data_pool_hot$imm_time, "6m", "1y", "2y"), y = prokariotic_biotas)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3","darkolivegreen") ) +
    labs(title = "Deployment season : Hot",
         x = "Immersion time",
         y = "Percentage cover of prokariotic_biotas") +
    scale_x_discrete(labels=time) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed) 
  
  res.aov <- rstatix::kruskal_test(data_pool_cool, prokariotic_biotas ~ imm_time)
  p.sed <- rstatix::wilcox_test(data_pool_cool, prokariotic_biotas ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y2_prokariotic_biotas <- ggplot(data_pool_cool, aes(x = fct_relevel(data_pool_cool$imm_time, "6m", "1y"), y = prokariotic_biotas)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3")) +
    labs(title = "Deployment season : Cool",
         x = "Immersion time",
         y = "Percentage cover of prokariotic_biotas") +
    scale_x_discrete(labels=c("6 month", "1 year")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)  
  
  res.aov <- rstatix::kruskal_test(data_pool_six, prokariotic_biotas ~ immersion_season)
  p.sed <- rstatix::wilcox_test(data_pool_six, prokariotic_biotas ~ immersion_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y3_prokariotic_biotas <- ggplot(data_pool_six, aes(x = fct_relevel(data_pool_six$immersion_season, "imm_cold", "imm_hot"), y = prokariotic_biotas)) +
    geom_boxplot(fill =  c("dodgerblue","firebrick3")) +
    labs(title = "Immersion time : 6 months",
         x = "Deployment season",
         y = "Percentage cover of prokariotic_biotas") +
    scale_x_discrete(labels=c("cool", "hot")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  res.aov <- rstatix::kruskal_test(data_pool_one, prokariotic_biotas ~ immersion_season)
  p.sed <- rstatix::wilcox_test(data_pool_one, prokariotic_biotas ~ immersion_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y4_prokariotic_biotas <- ggplot(data_pool_one, aes(x = fct_relevel(data_pool_one$immersion_season, "imm_cold", "imm_hot"), y = prokariotic_biotas)) +
    geom_boxplot(fill =  c("dodgerblue","firebrick3")) +
    labs(title = "Immersion time : one year",
         x = "Deployment season",
         y = "Percentage cover of prokariotic_biotas") +
    scale_x_discrete(labels=c("cool", "hot")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  #### Retrieval season
  
  res.aov <- rstatix::kruskal_test(data_pool_r_hot, prokariotic_biotas ~ imm_time)
  p.sed <- rstatix::wilcox_test(data_pool_r_hot, prokariotic_biotas ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y5_prokariotic_biotas <- ggplot(data_pool_r_hot, aes(x = fct_relevel(data_pool_r_hot$imm_time, "6m", "1y", "2y"), y = prokariotic_biotas)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3","darkolivegreen") ) +
    labs(title = "Retrieval season : Hot",
         x = "Immersion time",
         y = "Percentage cover of prokariotic_biotas") +
    scale_x_discrete(labels=time) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  res.aov <- rstatix::kruskal_test(data_pool_r_cool, prokariotic_biotas ~ imm_time)
  p.sed <- rstatix::wilcox_test(data_pool_r_cool, prokariotic_biotas ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y6_prokariotic_biotas <- ggplot(data_pool_r_cool, aes(x = fct_relevel(imm_time, "6m", "1y"), y = prokariotic_biotas)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3")) +
    labs(title = "Rerieval season : Cool",
         x = "Immersion time",
         y = "Percentage cover of prokariotic_biotas") +
    scale_x_discrete(labels=c("6 month", "1 year")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed) 
  
  res.aov <- rstatix::kruskal_test(data_pool_six, prokariotic_biotas ~ recovery_season)
  
  p.sed <- rstatix::wilcox_test(data_pool_six, prokariotic_biotas ~ recovery_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y7_prokariotic_biotas <- ggplot(data_pool_six, aes(x = fct_relevel(data_pool_six$recovery_season, "rec_cold", "rec_hot"), y = prokariotic_biotas)) +
    geom_boxplot(fill =  c("dodgerblue","firebrick3")) +
    labs(title = "Immersion time : 6 months",
         x = "retrieval season",
         y = "Percentage cover of prokariotic_biotas") +
    scale_x_discrete(labels=c("cool", "hot")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  res.aov <- rstatix::kruskal_test(data_pool_one, prokariotic_biotas ~ recovery_season)
  p.sed <- rstatix::wilcox_test(data_pool_one, prokariotic_biotas ~ recovery_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.1)
  
  y8_prokariotic_biotas <- ggplot(data_pool_one, aes(x = fct_relevel(data_pool_one$recovery_season, "rec_cold", "rec_hot"), y = prokariotic_biotas)) +
    geom_boxplot(fill =  c("dodgerblue","firebrick3")) +
    labs(title = "Immersion time : one year",
         x = "retrieval season",
         y = "Percentage cover of prokariotic_biotas") +
    scale_x_discrete(labels=c("cool", "hot")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed) 
  
  
  
  second_sol <- cowplot::plot_grid(y1_prokariotic_biotas, y2_prokariotic_biotas, y3_prokariotic_biotas, y4_prokariotic_biotas,
                                   y5_prokariotic_biotas, y6_prokariotic_biotas, y7_prokariotic_biotas, y8_prokariotic_biotas,
                                   labels = c(" "," "," "," "," "," "," "," "),
                                   ncol = 4,
                                   nrow = 2)
  
  
  path_to_boxplot_alt2 <- paste0("outputs/box/boxplot_prokariotic_biotas.pdf")
  ggsave(filename =  path_to_boxplot_alt2, plot = second_sol, width = 12, height = 8)
  
  
  
  #### final plot ####
  
  second_sol1 <- cowplot::plot_grid(y1_porifera, y2_porifera, y3_porifera, y4_porifera,
                                   y1_Ascidiacea, y2_Ascidiacea, y3_Ascidiacea, y4_Ascidiacea,
                                   y1_foraminifera, y2_foraminifera, y3_foraminifera, y4_foraminifera,
                                   y1_Hydrozoa, y2_Hydrozoa, y3_Hydrozoa, y4_Hydrozoa,
                                   y1_annelida, y2_annelida, y3_annelida, y4_annelida,
                                   y1_bryozoa, y2_bryozoa, y3_bryozoa, y4_bryozoa,
                                   ncol = 4,
                                   nrow = 6)
  
  path_to_boxplot_list1 <- paste0("outputs/box/boxplot_list1.pdf")
  ggsave(filename =  path_to_boxplot_list1, plot = second_sol1, width = 16, height = 22)
  
  second_sol2 <- cowplot::plot_grid(y1_Bivalvia, y2_Bivalvia, y3_Bivalvia, y4_Bivalvia,
                                   y1_CCA, y2_CCA, y3_CCA, y4_CCA,
                                   y1_other_algae, y2_other_algae, y3_other_algae, y4_other_algae,
                                   y1_prokariotic_biotas, y2_prokariotic_biotas, y3_prokariotic_biotas, y4_prokariotic_biotas,
                                   y1_sediment, y2_sediment, y3_sediment, y4_sediment,
                                   y1_bare_plate, y2_bare_plate, y3_bare_plate, y4_bare_plate,
                                   ncol = 4,
                                   nrow = 6)
  
  path_to_boxplot_list2 <- paste0("outputs/box/boxplot_list2.pdf")
  ggsave(filename =  path_to_boxplot_list2, plot = second_sol2, width = 16, height = 22)
  
  ####return####
return(c(path_to_boxplot_alt2))

  
  
  
}

