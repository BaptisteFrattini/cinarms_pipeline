#' Diversity analysis
#'
#' @param metadata_data_mean the path to the raw data file
#'
#' @return 
#' @export
#' 

diversity_explo <- function(metadata_data_mean){
  
  # metadata_data_mean = targets::tar_read(mean_metadata_data)
  
  #### Load data and meta data ####
  
  df_mean <- read.csv(metadata_data_mean[!grepl("metadata", metadata_data_mean)], header = TRUE)
  meta_mean <- read.csv(metadata_data_mean[grepl("metadata", metadata_data_mean)], header = TRUE)
  
  
  S <- vegan::specnumber(df_mean)
  
  data_div <- data.frame(s = S,
                         name = meta_mean$arms,
                         imm_time = meta_mean$imm_time,
                         imm_seas = meta_mean$immersion_season,
                         imm_ret = meta_mean$recovery_season)
  
  tapply(S, substr(meta_mean$arms, 1, 5), mean)
  tapply(S, meta_mean$imm_time, mean)
  time <- c("6 month", "1 year", "2 years")
  
  library(vctrs)
  library(ggpubr)
  library(forcats)
  library(ggplot2)
  
  ggplot(data_div, aes(x=s)) + 
    geom_density() #Normality OK
  ggpubr::ggqqplot(data_div$s) #Normality ok
  shapiro.test(data_div$s) #Normality OK
  
  #distrib not normal
  ?rstatix::kruskal_test
  # res.seas <-  p.sed <- rstatix::pairwise_t_test(data_div, s ~ imm_ret)
  # res.seas <-  p.sed <- rstatix::pairwise_t_test(data_div, s ~ imm_seas)
  res.aov <- rstatix::kruskal_test(data_div, s ~ imm_time)
  p.sed <- rstatix::wilcox_test(data_div, s ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.05)
  tapply(data_div$s, data_div$imm_time, mean)
  tapply(data_div$s, data_div$imm_time, sd)
  
  library(insight)
  library(parameters)
  library(effectsize)
  library(pwr)
  
  n <- pwr.anova.test(k = 3, f = 0.99, sig.level = 0.05, power = 0.8)
  model <- aov(s ~ imm_time, data = data_div)
  parameters(model, effectsize_type = c("eta","f"))
  
  library(lmPerm)
  model <- aovp(s ~ imm_time, data = data_div, perm = "Prob")
  summary(model)
  

  library(rcompanion)
  
  PT = pairwisePermutationTest(s ~ imm_time, data = data_div,
                               method="bonferroni")
  
  p.sed$p <- round(PT$p.adjust,3)
  
  
  ?pairwisePercentileTest
  #Experimental design is quite ok to run ANOVA
  
  
  
  # Define colors
  colors <- c("6m" = "#CC66CC", "1y" = "#1B9E77", "2y" = "#FF7F00")
  
  # Calculate mean for each immersion time
  means <- data_div %>%
    group_by(imm_time) %>%
    summarise(mean_s = mean(s), .groups = "drop")
  
  ff <- ggplot(data_div, aes(x = fct_relevel(imm_time, "6m", "1y", "2y"), y = s, color = imm_time)) +
    geom_jitter(width = 0.1, size = 3, alpha = 0.7) + # Points with individual colors
    geom_point(data = means, aes(x = imm_time, y = mean_s), color = "grey30", size = 4) + # Mean points in grey
    scale_color_manual(values = colors) + # Assign colors
    labs(title = "",
         x = "Immersion time",
         y = "Average MSP richness in an ARMS") +
    theme_classic() +
    theme(legend.position = "none") + # Remove legend
    stat_pvalue_manual(p.sed, label = "p") +
    annotate(geom="text", x=1, y=30, label = paste0("N = ",sum(data_div$imm_time == "6m")),
             color="black") +
    annotate(geom="text", x=2, y=32, label = paste0("N = ",sum(data_div$imm_time == "1y")),
             color="black") +
    annotate(geom="text", x=3, y=35, label = paste0("N = ",sum(data_div$imm_time == "2y")),
             color="black")
  
  # 
  path_to_box_div <- paste0("outputs/box_div.pdf")
  ggsave(filename =  path_to_box_div , width = 6, height = 6)
  
  #### What are the species common to all stage of colonisation ? ####
  df_pa <- vegan::decostand(df_mean, "pa")
  
  # Trouver les noms des colonnes qui ne présentent pas de 0
  names(df_pa)[colSums(df_pa == 0) == 0]
  
  
  return(path_to_box_div)
      
}  

