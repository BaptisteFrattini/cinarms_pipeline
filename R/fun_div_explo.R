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
  
  tapply(S, meta_mean$imm_time, mean)
  time <- c("6 month", "1 year", "2 years")
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
  res.aov <- rstatix::anova_test(data_div, s ~ imm_time)
  p.sed <- rstatix::pairwise_t_test(data_div, s ~ imm_time, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.05)
  tapply(data_div$s, data_div$imm_time, mean)
  tapply(data_div$s, data_div$imm_time, sd)
  
  library(parameters)
  library(effectsize)
  library(pwr)
  
  n <- pwr.anova.test(k = 3, f = 0.99, sig.level = 0.05, power = 0.8)
  model <- aov(s ~ imm_time, data = data_div)
  parameters(model, effectsize_type = c("eta","f"))
  
  #Experimental design is quite ok to run ANOVA
  
  ff <- ggplot(data_div, aes(x = fct_relevel(imm_time, "6m", "1y", "2y"), y = S)) +
    geom_boxplot(fill =  c("cadetblue2","cadetblue3","cadetblue4")) +
    labs(title = "",
         x = "Immersion time",
         y = "Average species richness in an ARMS") +
    theme_classic() +
    stat_pvalue_manual(p.sed) +
    annotate(geom="text", x=1, y=24, label = paste0("N = ",length(data_div$imm_time[grepl("6m", data_div$imm_time)])),
             color="black") +
    annotate(geom="text", x=2, y=27.75, label = paste0("N = ",length(data_div$imm_time[grepl("1y", data_div$imm_time)])),
             color="black") +
    annotate(geom="text", x=3, y=33, label = paste0("N = ",length(data_div$imm_time[grepl("2y", data_div$imm_time)])),
             color="black")
  
  path_to_box_div <- paste0("outputs/box_div.pdf")
  ggsave(filename =  path_to_box_div , width = 6, height = 6)
  
  #### What are the species common to all stage of colonisation ? ####
  df_pa <- vegan::decostand(df_mean, "pa")
  
  # Trouver les noms des colonnes qui ne présentent pas de 0
  names(df_pa)[colSums(df_pa == 0) == 0]
  
  
  return(path_to_box_div)
      
}  
