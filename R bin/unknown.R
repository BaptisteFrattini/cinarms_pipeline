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
         x = "Deployment season",
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
         x = "Deployment season",
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

  return(NULL)  
}