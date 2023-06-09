#' boxplot
#'
#' @param metadata_data_mean the path to the raw data file
#'
#' @return 
#' @export
#' 

boxplot_explo <- function(data_full_pool, meta_data){
 
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
  #### Bare plate ####
  
  ggplot(data_pool, aes(x=bare_plate)) + 
    geom_density()                        #Data not normal
  ggpubr::ggqqplot(data_pool$bare_plate)  #Data not normal
  shapiro.test(data_pool$bare_plate)      #Data not normal
  res.aov2 <- aov(bare_plate ~ imm_time*immersion_season, data=data_pool)
  residuals <- resid(res.aov2)
  plot(data_pool$bare_plate, residuals, xlab="Measurement", ylab="Residuals") # really bad
  abline(0,0)
  ggplot() +
    geom_qq(aes(sample = rstandard(mod))) +
    geom_abline(color = "red") +
    coord_fixed() # QQplot on residuals --> not OK
  
  # # #If distrib was normal 
  # res.aov <- rstatix::anova_test(data_pool, bare_plate ~ imm_time*immersion_season)
  # 
  # #plot the interaction
  # interaction.plot(x.factor = meta$recovery_season,
  #                  trace.factor = meta$imm_time,
  #                  response = data_pool$bare_plate)
  #distrib not normal
  
  library(lmPerm)
  mod <- aovp(bare_plate ~ imm_time*immersion_season,
              data = data_pool,
              perm="Exact") #Anova par permutation
  summary(mod)
  # checking residuals
  residuals <- resid(mod)
  plot(data_pool$bare_plate,
       residuals,
       xlab="Measurement",
       ylab="Residuals") # not OK
  abline(0,0)
  
  ggplot() +
    geom_qq(aes(sample = rstandard(mod))) +
    geom_abline(color = "red") +
    coord_fixed() # QQplot on residuals --> not OK
  
  res.aov <- rstatix::kruskal_test(data_pool, bare_plate ~ meta$imm_time)
  p.sed <- rstatix::wilcox_test(data_pool, bare_plate ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.05)
  
  a1_bare_plate <- ggplot(data_pool, aes(x = fct_relevel(meta$imm_time, "6m", "1y", "2y"), y = bare_plate)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3","darkolivegreen") ) +
    labs(title = "",
         x = "Immersion time",
         y = "Percentage cover of bare_plate") +
    scale_x_discrete(labels=time) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)


  p.sed <- rstatix::wilcox_test(data_pool, bare_plate ~ immersion_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.5)
  
  a2_bare_plate <- ggplot(data_pool, aes(x = meta$immersion_season, y = bare_plate)) +
    geom_boxplot(fill = c("dodgerblue2","firebrick3")) +
    labs(title = "",
         x = "Immersion season",
         y = "Percentage cover of bare_plate") +
    scale_x_discrete(labels = c("Cool", "Hot")) +
    theme(legend.position = "none") +
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  
  p.sed <- rstatix::wilcox_test(data_pool, bare_plate ~ recovery_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.5)

  a3_bare_plate <- ggplot(data_pool, aes(x = meta$recovery_season, y = bare_plate)) +
    geom_boxplot(fill = c("dodgerblue2","firebrick3")) +
    labs(title = "",
         x = "Recovery season",
         y = "Percentage cover of bare_plate") +
    scale_x_discrete(labels = c("Cool", "Hot")) +
    theme(legend.position = "none") +
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  first_row_bare_plate <- cowplot::plot_grid(a1_bare_plate, labels = c("Bare plate"))
  second_row_bare_plate <- cowplot::plot_grid(a2_bare_plate, a3_bare_plate, labels = c(" ", " "))
  gg_all_bare_plate = cowplot::plot_grid(first_row_bare_plate, second_row_bare_plate, labels=c('', ''), ncol=1)
  
  
  #### Porifera ####
  ggplot(data_pool, aes(x=porifera)) + 
    geom_density()
  ggpubr::ggqqplot(data_pool$porifera)
  shapiro.test(data_pool$porifera)
  
  #distrib not normal
  
  res.aov <- rstatix::kruskal_test(data_pool, porifera ~ meta$imm_time)
  p.sed <- rstatix::wilcox_test(data_pool, porifera ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.05)
  
  a1_porifera <- ggplot(data_pool, aes(x = fct_relevel(meta$imm_time, "6m", "1y", "2y"), y = porifera)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3","darkolivegreen") ) +
    labs(title = "",
         x = "Immersion time",
         y = "Percentage cover of porifera") +
    scale_x_discrete(labels=time) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  p.sed <- rstatix::wilcox_test(data_pool, porifera ~ immersion_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.5)
  
  a2_porifera <- ggplot(data_pool, aes(x = meta$immersion_season, y = porifera)) +
    geom_boxplot(fill = c("dodgerblue2","firebrick3")) +
    labs(title = "",
         x = "Immersion season",
         y = "Percentage cover of porifera") +
    scale_x_discrete(labels = c("Cool", "Hot")) +
    theme(legend.position = "none") +
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  
  p.sed <- rstatix::wilcox_test(data_pool, porifera ~ recovery_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.5)
  
  a3_porifera <- ggplot(data_pool, aes(x = meta$recovery_season, y = porifera)) +
    geom_boxplot(fill = c("dodgerblue2","firebrick3")) +
    labs(title = "",
         x = "Recovery season",
         y = "Percentage cover of porifera") +
    scale_x_discrete(labels = c("Cool", "Hot")) +
    theme(legend.position = "none") +
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  first_row_porifera <- cowplot::plot_grid(a1_porifera, labels = c("Porifera"))
  second_row_porifera <- cowplot::plot_grid(a2_porifera, a3_porifera, labels = c(" ", " "))
  gg_all_porifera = cowplot::plot_grid(first_row_porifera, second_row_porifera, labels=c('', ''), ncol=1)
  
  
  #### Ascidiacea colonial####
  ggplot(data_pool, aes(x=ascidiacea_c)) + 
    geom_density()
  ggpubr::ggqqplot(data_pool$ascidiacea_c)
  shapiro.test(data_pool$ascidiacea_c)
  
  #distrib not normal
  ?rstatix::kruskal_test
  res.aov <- rstatix::kruskal_test(data_pool, ascidiacea_c ~ meta$imm_time)
  p.sed <- rstatix::wilcox_test(data_pool, ascidiacea_c ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.05)
  
  a1_ascidiacea_c <- ggplot(data_pool, aes(x = fct_relevel(meta$imm_time, "6m", "1y", "2y"), y = ascidiacea_c)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3","darkolivegreen") ) +
    labs(title = "",
         x = "Immersion time",
         y = "Percentage cover of ascidiacea_c") +
    scale_x_discrete(labels=time) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)

  
  p.sed <- rstatix::wilcox_test(data_pool, ascidiacea_c ~ immersion_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.5)
  
  a2_ascidiacea_c <- ggplot(data_pool, aes(x = meta$immersion_season, y = ascidiacea_c)) +
    geom_boxplot(fill = c("dodgerblue2","firebrick3")) +
    labs(title = "",
         x = "Immersion season",
         y = "Percentage cover of ascidiacea_c") +
    scale_x_discrete(labels = c("Cool", "Hot")) +
    theme(legend.position = "none") +
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  
  p.sed <- rstatix::wilcox_test(data_pool, ascidiacea_c ~ recovery_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.5)
  
  a3_ascidiacea_c <- ggplot(data_pool, aes(x = meta$recovery_season, y = ascidiacea_c)) +
    geom_boxplot(fill = c("dodgerblue2","firebrick3")) +
    labs(title = "",
         x = "Recovery season",
         y = "Percentage cover of ascidiacea_c") +
    scale_x_discrete(labels = c("Cool", "Hot")) +
    theme(legend.position = "none") +
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  first_row_ascidiacea_c <- cowplot::plot_grid(a1_ascidiacea_c, labels = c("ascidiacea_c"))
  second_row_ascidiacea_c <- cowplot::plot_grid(a2_ascidiacea_c, a3_ascidiacea_c, labels = c(" ", " "))
  gg_all_ascidiacea_c = cowplot::plot_grid(first_row_ascidiacea_c, second_row_ascidiacea_c, labels=c('', ''), ncol=1)
  #### Ascidiacea solitary ####
  ggplot(data_pool, aes(x=ascidiacea_s)) + 
    geom_density()
  ggpubr::ggqqplot(data_pool$ascidiacea_s)
  shapiro.test(data_pool$ascidiacea_s)
  
  #distrib not normal
  ?rstatix::kruskal_test
  res.aov <- rstatix::kruskal_test(data_pool, ascidiacea_s ~ meta$imm_time)
  p.sed <- rstatix::wilcox_test(data_pool, ascidiacea_s ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.05)
  
  a1_ascidiacea_s <- ggplot(data_pool, aes(x = fct_relevel(meta$imm_time, "6m", "1y", "2y"), y = ascidiacea_s)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3","darkolivegreen") ) +
    labs(title = "",
         x = "Immersion time",
         y = "Percentage cover of ascidiacea_s") +
    scale_x_discrete(labels=time) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  
  p.sed <- rstatix::wilcox_test(data_pool, ascidiacea_s ~ immersion_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.5)
  
  a2_ascidiacea_s <- ggplot(data_pool, aes(x = meta$immersion_season, y = ascidiacea_s)) +
    geom_boxplot(fill = c("dodgerblue2","firebrick3")) +
    labs(title = "",
         x = "Immersion season",
         y = "Percentage cover of ascidiacea_s") +
    scale_x_discrete(labels = c("Cool", "Hot")) +
    theme(legend.position = "none") +
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  
  p.sed <- rstatix::wilcox_test(data_pool, ascidiacea_s ~ recovery_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.5)
  
  a3_ascidiacea_s <- ggplot(data_pool, aes(x = meta$recovery_season, y = ascidiacea_s)) +
    geom_boxplot(fill = c("dodgerblue2","firebrick3")) +
    labs(title = "",
         x = "Recovery season",
         y = "Percentage cover of ascidiacea_s") +
    scale_x_discrete(labels = c("Cool", "Hot")) +
    theme(legend.position = "none") +
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  first_row_ascidiacea_s <- cowplot::plot_grid(a1_ascidiacea_s, labels = c("ascidiacea_s"))
  second_row_ascidiacea_s <- cowplot::plot_grid(a2_ascidiacea_s, a3_ascidiacea_s, labels = c(" ", " "))
  gg_all_ascidiacea_s = cowplot::plot_grid(first_row_ascidiacea_s, second_row_ascidiacea_s, labels=c('', ''), ncol=1)

  #### Bryozoa ####
  ggplot(data_pool, aes(x=bryozoa)) + 
    geom_density()
  ggpubr::ggqqplot(data_pool$bryozoa)
  shapiro.test(data_pool$bryozoa)
  
  #distrib not normal
  
  res.aov <- rstatix::kruskal_test(data_pool, bryozoa ~ meta$imm_time)
  p.sed <- rstatix::wilcox_test(data_pool, bryozoa ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.05)
  
  a1_bryozoa <- ggplot(data_pool, aes(x = fct_relevel(meta$imm_time, "6m", "1y", "2y"), y = bryozoa)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3","darkolivegreen") ) +
    labs(title = "",
         x = "Immersion time",
         y = "Percentage cover of bryozoa") +
    scale_x_discrete(labels=time) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  
  p.sed <- rstatix::wilcox_test(data_pool, bryozoa ~ immersion_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.5)
  
  a2_bryozoa <- ggplot(data_pool, aes(x = meta$immersion_season, y = bryozoa)) +
    geom_boxplot(fill = c("dodgerblue2","firebrick3")) +
    labs(title = "",
         x = "Immersion season",
         y = "Percentage cover of bryozoa") +
    scale_x_discrete(labels = c("Cool", "Hot")) +
    theme(legend.position = "none") +
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  
  p.sed <- rstatix::wilcox_test(data_pool, bryozoa ~ recovery_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.5)
  
  a3_bryozoa <- ggplot(data_pool, aes(x = meta$recovery_season, y = bryozoa)) +
    geom_boxplot(fill = c("dodgerblue2","firebrick3")) +
    labs(title = "",
         x = "Recovery season",
         y = "Percentage cover of bryozoa") +
    scale_x_discrete(labels = c("Cool", "Hot")) +
    theme(legend.position = "none") +
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  first_row_bryozoa <- cowplot::plot_grid(a1_bryozoa, labels = c("Bryozoa"))
  second_row_bryozoa <- cowplot::plot_grid(a2_bryozoa, a3_bryozoa, labels = c(" ", " "))
  gg_all_bryozoa = cowplot::plot_grid(first_row_bryozoa, second_row_bryozoa, labels=c('', ''), ncol=1)
  
  
  #### Annelida ####
  ggplot(data_pool, aes(x=annelida)) + 
    geom_density()
  ggpubr::ggqqplot(data_pool$annelida)
  shapiro.test(data_pool$annelida)
  
  # #If distrib was normal
  res.aov <- rstatix::anova_test(data_pool, annelida ~ imm_time*immersion_season)
  # #plot the interaction
  interaction.plot(x.factor = meta$immersion_season,
                   trace.factor = meta$imm_time,
                   response = data_pool$annelida)

  
  #distrib not normal
  
  res.aov <- rstatix::kruskal_test(data_pool, annelida ~ meta$imm_time)
  p.sed <- rstatix::wilcox_test(data_pool, annelida ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.05)
  
  a1_annelida <- ggplot(data_pool, aes(x = fct_relevel(meta$imm_time, "6m", "1y", "2y"), y = annelida)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3","darkolivegreen") ) +
    labs(title = "",
         x = "Immersion time",
         y = "Percentage cover of annelida") +
    scale_x_discrete(labels=time) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)

  p.sed <- rstatix::wilcox_test(data_pool, annelida ~ immersion_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.5)
  
  a2_annelida <- ggplot(data_pool, aes(x = meta$immersion_season, y = annelida)) +
    geom_boxplot(fill = c("dodgerblue2","firebrick3")) +
    labs(title = "",
         x = "Immersion season",
         y = "Percentage cover of annelida") +
    scale_x_discrete(labels = c("Cool", "Hot")) +
    theme(legend.position = "none") +
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  
  p.sed <- rstatix::wilcox_test(data_pool, annelida ~ recovery_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.5)
  
  a3_annelida <- ggplot(data_pool, aes(x = meta$recovery_season, y = annelida)) +
    geom_boxplot(fill = c("dodgerblue2","firebrick3")) +
    labs(title = "",
         x = "Recovery season",
         y = "Percentage cover of annelida") +
    scale_x_discrete(labels = c("Cool", "Hot")) +
    theme(legend.position = "none") +
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  first_row_annelida <- cowplot::plot_grid(a1_annelida, labels = c("Annelida"))
  second_row_annelida <- cowplot::plot_grid(a2_annelida, a3_annelida, labels = c(" ", " "))
  gg_all_annelida = cowplot::plot_grid(first_row_annelida, second_row_annelida, labels=c('', ''), ncol=1)
  
  #### foraminifera ####
  ggplot(data_pool, aes(x=foraminifera)) + 
    geom_density()
  ggpubr::ggqqplot(data_pool$foraminifera)
  shapiro.test(data_pool$foraminifera)
  
  #distrib not normal
  
  res.aov <- rstatix::kruskal_test(data_pool, foraminifera ~ meta$imm_time)
  p.sed <- rstatix::wilcox_test(data_pool, foraminifera ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.05)
  
  a1_foraminifera <- ggplot(data_pool, aes(x = fct_relevel(meta$imm_time, "6m", "1y", "2y"), y = foraminifera)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3","darkolivegreen") ) +
    labs(title = "",
         x = "Immersion time",
         y = "Percentage cover of foraminifera") +
    scale_x_discrete(labels=time) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  
  p.sed <- rstatix::wilcox_test(data_pool, foraminifera ~ immersion_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.5)
  
  a2_foraminifera <- ggplot(data_pool, aes(x = meta$immersion_season, y = foraminifera)) +
    geom_boxplot(fill = c("dodgerblue2","firebrick3")) +
    labs(title = "",
         x = "Immersion season",
         y = "Percentage cover of foraminifera") +
    scale_x_discrete(labels = c("Cool", "Hot")) +
    theme(legend.position = "none") +
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  
  p.sed <- rstatix::wilcox_test(data_pool, foraminifera ~ recovery_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.5)
  
  a3_foraminifera <- ggplot(data_pool, aes(x = meta$recovery_season, y = foraminifera)) +
    geom_boxplot(fill = c("dodgerblue2","firebrick3")) +
    labs(title = "",
         x = "Recovery season",
         y = "Percentage cover of foraminifera") +
    scale_x_discrete(labels = c("Cool", "Hot")) +
    theme(legend.position = "none") +
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  first_row_foraminifera <- cowplot::plot_grid(a1_foraminifera, labels = c("Foraminifera"))
  second_row_foraminifera <- cowplot::plot_grid(a2_foraminifera, a3_foraminifera, labels = c(" ", " "))
  gg_all_foraminifera = cowplot::plot_grid(first_row_foraminifera, second_row_foraminifera, labels=c('', ''), ncol=1)
  

  
  
  #### Hydrozoa ####
  ggplot(data_pool, aes(x=Hydrozoa)) + 
    geom_density()
  ggpubr::ggqqplot(data_pool$Hydrozoa)
  shapiro.test(data_pool$Hydrozoa)
  
  # #If distrib was normal 
  # res.aov <- rstatix::anova_test(data_pool, Hydrozoa ~ imm_time*immersion_season)
  # #plot the interaction
  # interaction.plot(x.factor = meta$immersion_season,
  #                  trace.factor = meta$imm_time,
  #                  response = data_pool$Hydrozoa)
  
  
  #distrib not normal
  
  res.aov <- rstatix::kruskal_test(data_pool, Hydrozoa ~ meta$imm_time)
  p.sed <- rstatix::wilcox_test(data_pool, Hydrozoa ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.05)
  
  a1_Hydrozoa <- ggplot(data_pool, aes(x = fct_relevel(meta$imm_time, "6m", "1y", "2y"), y = Hydrozoa)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3","darkolivegreen") ) +
    labs(title = "",
         x = "Immersion time",
         y = "Percentage cover of Hydrozoa") +
    scale_x_discrete(labels=time) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  
  p.sed <- rstatix::wilcox_test(data_pool, Hydrozoa ~ immersion_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.5)
  
  a2_Hydrozoa <- ggplot(data_pool, aes(x = meta$immersion_season, y = Hydrozoa)) +
    geom_boxplot(fill = c("dodgerblue2","firebrick3")) +
    labs(title = "",
         x = "Immersion season",
         y = "Percentage cover of Hydrozoa") +
    scale_x_discrete(labels = c("Cool", "Hot")) +
    theme(legend.position = "none") +
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  
  p.sed <- rstatix::wilcox_test(data_pool, Hydrozoa ~ recovery_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.5)
  
  a3_Hydrozoa <- ggplot(data_pool, aes(x = meta$recovery_season, y = Hydrozoa)) +
    geom_boxplot(fill = c("dodgerblue2","firebrick3")) +
    labs(title = "",
         x = "Recovery season",
         y = "Percentage cover of Hydrozoa") +
    scale_x_discrete(labels = c("Cool", "Hot")) +
    theme(legend.position = "none") +
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  first_row_Hydrozoa <- cowplot::plot_grid(a1_Hydrozoa, labels = c("Hydrozoa"))
  second_row_Hydrozoa <- cowplot::plot_grid(a2_Hydrozoa, a3_Hydrozoa, labels = c(" ", " "))
  gg_all_Hydrozoa = cowplot::plot_grid(first_row_Hydrozoa, second_row_Hydrozoa, labels=c('', ''), ncol=1)
  
  #### CCA ####
  ggplot(data_pool, aes(x=CCA)) + 
    geom_density()
  ggpubr::ggqqplot(data_pool$CCA)
  shapiro.test(data_pool$CCA)
  
  #distrib not normal
  
  res.aov <- rstatix::anova_test(data_pool, CCA ~ imm_time*immersion_season)
  p.sed <- rstatix::pairwise_t_test(data_pool, CCA ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.05)
  
  a1_CCA <- ggplot(data_pool, aes(x = fct_relevel(meta$imm_time, "6m", "1y", "2y"), y = CCA)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3","darkolivegreen") ) +
    labs(title = "",
         x = "Immersion time",
         y = "Percentage cover of CCA") +
    scale_x_discrete(labels=time) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  #plot the interaction
  # ggplot(data_pool, aes(x = meta$imm_time, y = CCA, color = meta$immersion_season)) +
  #   geom_point() +
  #   labs(x = "Factor 1", y = "Response", color = "Factor 2") +
  #   ggtitle("Interaction Plot") +
  #   theme_bw()
  
  p.sed <- rstatix::pairwise_t_test(data_pool, CCA ~ immersion_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.5)
  
  a2_CCA <- ggplot(data_pool, aes(x = meta$immersion_season, y = CCA)) +
    geom_boxplot(fill = c("dodgerblue2","firebrick3")) +
    labs(title = "",
         x = "Immersion season",
         y = "Percentage cover of CCA") +
    scale_x_discrete(labels = c("Cool", "Hot")) +
    theme(legend.position = "none") +
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  
  p.sed <- rstatix::pairwise_t_test(data_pool, CCA ~ recovery_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.5)
  
  a3_CCA <- ggplot(data_pool, aes(x = meta$recovery_season, y = CCA)) +
    geom_boxplot(fill = c("dodgerblue2","firebrick3")) +
    labs(title = "",
         x = "Recovery season",
         y = "Percentage cover of CCA") +
    scale_x_discrete(labels = c("Cool", "Hot")) +
    theme(legend.position = "none") +
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  first_row_CCA <- cowplot::plot_grid(a1_CCA, labels = c("CCA"))
  second_row_CCA <- cowplot::plot_grid(a2_CCA, a3_CCA, labels = c(" ", " "))
  gg_all_CCA = cowplot::plot_grid(first_row_CCA, second_row_CCA, labels=c('', ''), ncol=1)
  
  
  #### Bivalvia ####
  ggplot(data_pool, aes(x=Bivalvia)) + 
    geom_density()
  ggpubr::ggqqplot(data_pool$Bivalvia)
  shapiro.test(data_pool$Bivalvia)
  
  #distrib not normal
  
  res.aov <- rstatix::kruskal_test(data_pool, Bivalvia ~ meta$imm_time)
  p.sed <- rstatix::wilcox_test(data_pool, Bivalvia ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.05)
  
  a1_Bivalvia <- ggplot(data_pool, aes(x = fct_relevel(meta$imm_time, "6m", "1y", "2y"), y = Bivalvia)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3","darkolivegreen") ) +
    labs(title = "",
         x = "Immersion time",
         y = "Percentage cover of Bivalvia") +
    scale_x_discrete(labels=time) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  
  p.sed <- rstatix::wilcox_test(data_pool, Bivalvia ~ immersion_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.5)
  
  a2_Bivalvia <- ggplot(data_pool, aes(x = meta$immersion_season, y = Bivalvia)) +
    geom_boxplot(fill = c("dodgerblue2","firebrick3")) +
    labs(title = "",
         x = "Immersion season",
         y = "Percentage cover of Bivalvia") +
    scale_x_discrete(labels = c("Cool", "Hot")) +
    theme(legend.position = "none") +
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  
  p.sed <- rstatix::wilcox_test(data_pool, Bivalvia ~ recovery_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.5)
  
  a3_Bivalvia <- ggplot(data_pool, aes(x = meta$recovery_season, y = Bivalvia)) +
    geom_boxplot(fill = c("dodgerblue2","firebrick3")) +
    labs(title = "",
         x = "Recovery season",
         y = "Percentage cover of Bivalvia") +
    scale_x_discrete(labels = c("Cool", "Hot")) +
    theme(legend.position = "none") +
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  first_row_Bivalvia <- cowplot::plot_grid(a1_Bivalvia, labels = c("Bivalvia"))
  second_row_Bivalvia <- cowplot::plot_grid(a2_Bivalvia, a3_Bivalvia, labels = c(" ", " "))
  gg_all_Bivalvia = cowplot::plot_grid(first_row_Bivalvia, second_row_Bivalvia, labels=c('', ''), ncol=1)
  #### other_algae ####
  ggplot(data_pool, aes(x=other_algae)) + 
    geom_density()
  ggpubr::ggqqplot(data_pool$other_algae)
  shapiro.test(data_pool$other_algae)
  
  #distrib not normal
  
  res.aov <- rstatix::kruskal_test(data_pool, other_algae ~ meta$imm_time)
  p.sed <- rstatix::wilcox_test(data_pool, other_algae ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.05)
  
  a1_other_algae <- ggplot(data_pool, aes(x = fct_relevel(meta$imm_time, "6m", "1y", "2y"), y = other_algae)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3","darkolivegreen") ) +
    labs(title = "",
         x = "Immersion time",
         y = "Percentage cover of other_algae") +
    scale_x_discrete(labels=time) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  
  p.sed <- rstatix::wilcox_test(data_pool, other_algae ~ immersion_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.5)
  
  a2_other_algae <- ggplot(data_pool, aes(x = meta$immersion_season, y = other_algae)) +
    geom_boxplot(fill = c("dodgerblue2","firebrick3")) +
    labs(title = "",
         x = "Immersion season",
         y = "Percentage cover of other_algae") +
    scale_x_discrete(labels = c("Cool", "Hot")) +
    theme(legend.position = "none") +
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  
  p.sed <- rstatix::wilcox_test(data_pool, other_algae ~ recovery_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.5)
  
  a3_other_algae <- ggplot(data_pool, aes(x = meta$recovery_season, y = other_algae)) +
    geom_boxplot(fill = c("dodgerblue2","firebrick3")) +
    labs(title = "",
         x = "Recovery season",
         y = "Percentage cover of other_algae") +
    scale_x_discrete(labels = c("Cool", "Hot")) +
    theme(legend.position = "none") +
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  first_row_other_algae <- cowplot::plot_grid(a1_other_algae, labels = c("Other algae"))
  second_row_other_algae <- cowplot::plot_grid(a2_other_algae, a3_other_algae, labels = c(" ", " "))
  gg_all_other_algae = cowplot::plot_grid(first_row_other_algae, second_row_other_algae, labels=c('', ''), ncol=1)
  
  #### sediment ####
  ggplot(data_pool, aes(x=sediment)) + 
    geom_density()
  ggpubr::ggqqplot(data_pool$sediment)
  shapiro.test(data_pool$sediment)
  
  #distrib not normal
  
  res.aov <- rstatix::kruskal_test(data_pool, sediment ~ meta$imm_time)
  p.sed <- rstatix::wilcox_test(data_pool, sediment ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.05)
  
  a1_sediment <- ggplot(data_pool, aes(x = fct_relevel(meta$imm_time, "6m", "1y", "2y"), y = sediment)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3","darkolivegreen") ) +
    labs(title = "",
         x = "Immersion time",
         y = "Percentage cover of sediment") +
    scale_x_discrete(labels=time) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  
  p.sed <- rstatix::wilcox_test(data_pool, sediment ~ immersion_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.5)
  
  a2_sediment <- ggplot(data_pool, aes(x = meta$immersion_season, y = sediment)) +
    geom_boxplot(fill = c("dodgerblue2","firebrick3")) +
    labs(title = "",
         x = "Immersion season",
         y = "Percentage cover of sediment") +
    scale_x_discrete(labels = c("Cool", "Hot")) +
    theme(legend.position = "none") +
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  
  p.sed <- rstatix::wilcox_test(data_pool, sediment ~ recovery_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.5)
  
  a3_sediment <- ggplot(data_pool, aes(x = meta$recovery_season, y = sediment)) +
    geom_boxplot(fill = c("dodgerblue2","firebrick3")) +
    labs(title = "",
         x = "Recovery season",
         y = "Percentage cover of sediment") +
    scale_x_discrete(labels = c("Cool", "Hot")) +
    theme(legend.position = "none") +
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  first_row_sediment <- cowplot::plot_grid(a1_sediment, labels = c("Sediments"))
  second_row_sediment <- cowplot::plot_grid(a2_sediment, a3_sediment, labels = c(" ", " "))
  gg_all_sediment = cowplot::plot_grid(first_row_sediment, second_row_sediment, labels=c('', ''), ncol=1)
  
  #### prokariotic_biotas ####
  
  ggplot(data_pool, aes(x=prokariotic_biotas)) + 
    geom_density()
  ggpubr::ggqqplot(data_pool$prokariotic_biotas)
  shapiro.test(data_pool$prokariotic_biotas)
  
  #distrib not normal
  
  res.aov <- rstatix::kruskal_test(data_pool, prokariotic_biotas ~ meta$imm_time)
  p.sed <- rstatix::wilcox_test(data_pool, prokariotic_biotas ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.05)
  
  a1_prokariotic_biotas <- ggplot(data_pool, aes(x = fct_relevel(meta$imm_time, "6m", "1y", "2y"), y = prokariotic_biotas)) +
    geom_boxplot(fill =  c("darkolivegreen1","darkolivegreen3","darkolivegreen") ) +
    labs(title = "",
         x = "Immersion time",
         y = "Percentage cover of prokariotic_biotas") +
    scale_x_discrete(labels=time) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  
  p.sed <- rstatix::wilcox_test(data_pool, prokariotic_biotas ~ immersion_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.5)
  
  a2_prokariotic_biotas <- ggplot(data_pool, aes(x = meta$immersion_season, y = prokariotic_biotas)) +
    geom_boxplot(fill = c("dodgerblue2","firebrick3")) +
    labs(title = "",
         x = "Immersion season",
         y = "Percentage cover of prokariotic_biotas") +
    scale_x_discrete(labels = c("Cool", "Hot")) +
    theme(legend.position = "none") +
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  
  p.sed <- rstatix::wilcox_test(data_pool, prokariotic_biotas ~ recovery_season, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.5)
  
  a3_prokariotic_biotas <- ggplot(data_pool, aes(x = meta$recovery_season, y = prokariotic_biotas)) +
    geom_boxplot(fill = c("dodgerblue2","firebrick3")) +
    labs(title = "",
         x = "Recovery season",
         y = "Percentage cover of prokariotic_biotas") +
    scale_x_discrete(labels = c("Cool", "Hot")) +
    theme(legend.position = "none") +
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  first_row_prokariotic_biotas <- cowplot::plot_grid(a1_prokariotic_biotas, labels = c("Prokariotic biotas"))
  second_row_prokariotic_biotas <- cowplot::plot_grid(a2_prokariotic_biotas, a3_prokariotic_biotas, labels = c(" ", " "))
  gg_all_prokariotic_biotas = cowplot::plot_grid(first_row_prokariotic_biotas, second_row_prokariotic_biotas, labels=c('', ''), ncol=1)
  #### plot all ####
  
  
  fin <- cowplot::plot_grid(gg_all_CCA,
                            gg_all_foraminifera,
                            gg_all_annelida,
                            gg_all_bryozoa,
                            gg_all_ascidiacea_c,
                            gg_all_ascidiacea_s,
                            gg_all_Hydrozoa,
                            gg_all_Bivalvia,
                            gg_all_porifera,
                            gg_all_prokariotic_biotas,
                            gg_all_bare_plate,
                            gg_all_sediment,
                            ncol = 6,
                            nrow = 2)
  
  
  path_to_boxplot <- paste0("outputs/boxplot_pool.pdf")
  ggsave(filename =  path_to_boxplot , width = 27, height = 17)
  
  return(path_to_boxplot)
}


