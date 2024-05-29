

#' Eveness index computation
#'
#' @param metadata_data data
#' @param metadata
#' @return the path to the subseted raw data file
#' @export

fun_pielou <- function(metadata, metadata_data_mean){
 
  # metadata_data_mean = targets::tar_read(mean_metadata_data)
  # metadata = targets::tar_read(metadata_data)
  library(lmPerm)
  library(vegan)
  library(ggplot2)
  library(rstatix)
  library(ggpubr)
  library(forcats)
  library(canaper)
  
  # Pielou index on data mean by ARMS (ARMS/msp matrix) ####
  
  df_mean <- read.csv(metadata_data_mean[!grepl("metadata", metadata_data_mean)], header = TRUE)
  meta_mean <- read.csv(metadata_data_mean[grepl("metadata", metadata_data_mean)], header = TRUE)
  
    ## Compute Pielou's eveness ####

  H <- vegan::diversity(df_mean)
  
  J <- H/log(vegan::specnumber(df_mean))
  
  table = data.frame(arms = meta_mean$arms_name,
                     J = J)
  
      ### plot eveness ####
  
  mod <- aovp(table$J ~ table$arms,
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
  
  
  p.sed <- rstatix::wilcox_test(table, J ~ arms, p.adjust.method = "bonferroni")
  p.sed <- rstatix::add_y_position(test = p.sed, step.increase = 0.05)
  
  piel <- ggplot(table, aes(x = fct_relevel(arms, "CINA1", "CINA3", "CINA2", "CINA4", "RUNA2"), y = J)) +
    geom_boxplot(fill =  c("#CC66CC","#CC66CC","#1B9E77","#1B9E77","#FF7F00") ) +
    labs(title = "",
         x = "set",
         y = "Mean Eveness") +
    scale_x_discrete(labels= c("CINA1", "CINA3", "CINA2", "CINA4", "RUNA2")) +
    theme(legend.position = "none")+
    theme_classic() +
    stat_pvalue_manual(p.sed)
  
  path_to_boxplot_piel <- paste0("outputs/Pielou_mean.pdf")
  ggsave(filename =  path_to_boxplot_piel, plot = piel, width = 12, height = 9)
  

  
    ## Compute Pielou index on null model ####
  df.rounded <- round(df_mean*10000, digits = 0)
  df.rd <- canaper::cpr_rand_comm(df.rounded, "quasiswap_count", 10000)
  df.rd <- df.rd/10000
  rowSums(df_mean)
  rowSums(df.rd)
  colSums(df_mean)
  colSums(df.rd)
  vegan::specnumber(df.rd)
  vegan::specnumber(df_mean)
  ?vegan::commsim()
  
  df.rd <- as.data.frame(vegan::decostand(df.rd, "total"))*100
  
  table.rd = data.frame(set = meta_mean$arms_name,
                        J = J)
  
      ### plot SES for eveness ####
  
  ses <- (table.rd$J - mean(table$J))/sd(table.rd$J)
  
  table_ses <- data.frame(set = meta_mean$arms_name,
                          ses = ses)
  
  
  vv <- ggplot(table_ses, aes(x = fct_relevel(set, "CINA1", "CINA3", "CINA2", "CINA4", "RUNA2"), y = ses)) +
    geom_boxplot(fill =  c("#CC66CC","#CC66CC","#1B9E77","#1B9E77","#FF7F00") ) +
    labs(title = "SES for Pielou index",
         x = "Comparisons",
         y = "") +
    theme(legend.position = "none") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), 
          axis.title.x = element_blank(), 
          axis.title.y = element_blank(),
          axis.text.y = element_text(size = 12)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11), axis.title.x = element_blank(), axis.title.y = element_text(size=12)) +
    geom_hline(yintercept = -1.96, colour = "red") +
    geom_hline(yintercept = 1.96, colour = "red") +
    geom_hline(yintercept = 0, colour = "darkgrey")+
    stat_summary(fun=mean, colour="darkred", geom="point", 
                 shape=18, size=3, show.legend=FALSE)
  
  path_to_boxplot_piel_ses <- paste0("outputs/Pielou_SES_mean.pdf")
  ggsave(filename =  path_to_boxplot_piel_ses, plot = vv, width = 12, height = 9)
  
  
  
  # Pielou index on raw data (plate/msp matrix) ####
  
  df <- read.csv(metadata[!grepl("metadata", metadata)], header = TRUE)
  meta <- read.csv(metadata[grepl("metadata", metadata)], header = TRUE)
  
    ## Compute Pielou's eveness ####
  
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
  
  path_to_boxplot_piel2 <- paste0("outputs/Pielou_full.pdf")
  ggsave(filename =  path_to_boxplot_piel2, plot = piel2, width = 12, height = 9)
  
    ## Compute Pielou index on null model ####
  df.rounded <- round(df*10000, digits = 0)
  df.rd <- canaper::cpr_rand_comm(df.rounded, "quasiswap_count", 10000)
  df.rd <- df.rd/10000
  rowSums(df)
  rowSums(df.rd)
  colSums(df)
  colSums(df.rd)
  vegan::specnumber(df.rd)
  vegan::specnumber(df)
  ?vegan::commsim()
  
  df.rd <- as.data.frame(vegan::decostand(df.rd, "total"))*100

  table.rd = data.frame(set = meta$arms,
                     J = J)
  
      ### plot eveness ####

  piel_random <- ggplot(table.rd, aes(x = fct_relevel(set, "CINA1", "CINA3", "CINA2", "CINA4", "RUNA2"), y = J)) +
    geom_boxplot(fill =  c("#CC66CC","#CC66CC","#1B9E77","#1B9E77","#FF7F00") ) +
    labs(title = "",
         x = "set",
         y = "Mean Eveness") +
    scale_x_discrete(labels= c("CINA1", "CINA3", "CINA2", "CINA4", "RUNA2")) +
    theme(legend.position = "none") +
    theme_classic() 
  
  path_to_boxplot_piel2 <- paste0("outputs/piel2.pdf")
  ggsave(filename =  path_to_boxplot_piel2, plot = piel_random, width = 12, height = 9)
  
      ### plot SES for eveness ####
  
  ses <- (table.rd$J - mean(table$J))/sd(table.rd$J)
  
  table_ses <- data.frame(set = meta$arms,
                          ses = ses)
  
 
  jj <- ggplot(table_ses, aes(x = fct_relevel(set, "CINA1", "CINA3", "CINA2", "CINA4", "RUNA2"), y = ses)) +
    geom_boxplot(fill =  c("#CC66CC","#CC66CC","#1B9E77","#1B9E77","#FF7F00") ) +
    labs(title = "SES for Pielou index",
         x = "Comparisons",
         y = "") +
    theme(legend.position = "none") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), 
          axis.title.x = element_blank(), 
          axis.title.y = element_blank(),
          axis.text.y = element_text(size = 12)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11), axis.title.x = element_blank(), axis.title.y = element_text(size=12)) +
    geom_hline(yintercept = -1.96, colour = "red") +
    geom_hline(yintercept = 1.96, colour = "red") +
    geom_hline(yintercept = 0, colour = "darkgrey")+
    stat_summary(fun=mean, colour="darkred", geom="point", 
                 shape=18, size=3, show.legend=FALSE)
  
  path_to_boxplot_piel_ses <- paste0("outputs/Pielou_SES_full.pdf")
  ggsave(filename =  path_to_boxplot_piel_ses, plot = jj, width = 12, height = 9)
  
  # return ####
return(piel2)
  
  
  }
