

#' Eveness index computation
#'
#' @param metadata_data data
#' @param metadata
#' @return the path to the subseted raw data file
#' @export

fun_pielou <- function(metadata, metadata_data_mean){
 
  # metadata_data_mean = targets::tar_read(mean_metadata_data)
  library(lmPerm)
  library(vegan)
  library(ggplot2)
  library(rstatix)
  library(ggpubr)
  library(forcats)
  # Load data and meta data mean ####
  
  df_mean <- read.csv(metadata_data_mean[!grepl("metadata", metadata_data_mean)], header = TRUE)
  meta_mean <- read.csv(metadata_data_mean[grepl("metadata", metadata_data_mean)], header = TRUE)
  
  ## Pielou's eveness ####

  H <- vegan::diversity(df_mean)
  
  J <- H/log(vegan::specnumber(df_mean))
  
  table = data.frame(imm = meta_mean$imm_time,
                     J = J)
  
  ### plot eveness ####
  
  mod <- aovp(table$J ~ table$imm,
              data = table,
              perm="Exact") #Anova par permutation
  summary(mod)
  # checking residuals
  residuals <- resid(mod)
  plot(table$J,
       residuals,
       xlab="Measurement",
       ylab="Residuals") # not OK
  abline(0,0)
  
  ggplot() +
    geom_qq(aes(sample = rstandard(mod))) +
    geom_abline(color = "red") +
    coord_fixed() # QQplot on residuals --> not OK
  
  
  p.sed <- rstatix::wilcox_test(table, J ~ imm, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.05)
  
  piel <- ggplot(table, aes(x = fct_relevel(imm, "6m", "1y", "2y"), y = J)) +
    geom_boxplot(fill =  c("#CC66CC","#1B9E77","#FF7F00")) +
    labs(title = "",
         x = "set",
         y = "Mean Eveness") +
    scale_x_discrete(labels= c("6m", "1y", "2y")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  path_to_boxplot_piel <- paste0("outputs/piel.pdf")
  ggsave(filename =  path_to_boxplot_piel, plot = piel, width = 12, height = 9)
  
  
  # Load data and meta data all ####
  
  # metadata = targets::tar_read(metadata_data)
  
  df <- read.csv(metadata[!grepl("metadata", metadata)], header = TRUE)
  meta <- read.csv(metadata[grepl("metadata", metadata)], header = TRUE)
  
  ## Pielou's eveness ####
  
  H <- vegan::diversity(df)
  
  J <- H/log(vegan::specnumber(df))
  
  
  table = data.frame(set = meta$arms,
                     J = J)
  
  ### plot eveness ####
  
 
  mod <- aovp(table$J ~ table$set,
              data = table,
              perm="Exact") #Anova par permutation
  summary(mod)
  # checking residuals
  residuals <- resid(mod)
  plot(table$J,
       residuals,
       xlab="Measurement",
       ylab="Residuals") # not OK
  abline(0,0)
  
  ggplot() +
    geom_qq(aes(sample = rstandard(mod))) +
    geom_abline(color = "red") +
    coord_fixed() # QQplot on residuals --> not OK
  

  p.sed <- rstatix::t_test(table, J ~ set, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.05)
  
  piel2 <- ggplot(table, aes(x = fct_relevel(set, "CINA1", "CINA3", "CINA2", "CINA4", "RUNA2"), y = J)) +
    geom_boxplot(fill =  c("#CC66CC","#CC66CC","#1B9E77","#1B9E77","#FF7F00") ) +
    labs(title = "",
         x = "set",
         y = "Mean Eveness") +
    scale_x_discrete(labels= c("CINA1", "CINA3", "CINA2", "CINA4", "RUNA2")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  path_to_boxplot_piel2 <- paste0("outputs/piel2.pdf")
  ggsave(filename =  path_to_boxplot_piel2, plot = piel2, width = 12, height = 9)
  
return(piel2)
  
  
  }
