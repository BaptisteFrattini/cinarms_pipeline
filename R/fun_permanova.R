#' permanova & simper
#'
#' @param metadata_data_mean the path to the raw data file
#'
#' @return 
#' @export
#' 

fun_perm <- function(metadata_data_mean){
  
  metadata_data_mean = targets::tar_read(mean_metadata_data)
  
  #### Load data and meta data ####
  
  df_mean <- read.csv(metadata_data_mean[!grepl("metadata", metadata_data_mean)], header = TRUE)
  meta_mean <- read.csv(metadata_data_mean[grepl("metadata", metadata_data_mean)], header = TRUE)
  
  #### Compute permanova ####
  
  library(vegan)
  
  # Avec RUNA3 
  
  meta_mean$imm_rec <- paste(meta_mean$immersion_season, meta_mean$recovery_season)
  
  perm <- vegan::adonis2(df_mean ~ imm_time+immersion_season, data = meta_mean, method = "bray", permutations = 99999)

  
  # Sans RUNA3
  df_mean <- df_mean[1:15,]
  meta_mean <- meta_mean[1:15,]
  
  perm <- vegan::adonis2(df_mean ~ imm_time*immersion_season, data = meta_mean, method = "bray", permutations = 99999)
  
  # pairwise_ado <- pairwise.adonis2(df_mean ~ imm_time, data = meta_mean, method = "bray")

  
  sim_imm_tim <- summary(simper(df_mean, meta_mean$imm_time))
  
  sim_imm_tim_1 <- sim_imm_tim[[1]]
  contrib_imm_tim_1 <- sim_imm_tim_1[(sim_imm_tim_1$p < 0.05) | (sim_imm_tim_1$average > 0.1),]
  
  sim_imm_tim_2 <- sim_imm_tim[[2]]
  contrib_imm_tim_2 <- sim_imm_tim_2[(sim_imm_tim_2$p < 0.05) | (sim_imm_tim_2$average > 0.1),]
  
  sim_imm_tim_3 <- sim_imm_tim[[3]]
  contrib_imm_tim_3 <- sim_imm_tim_3[(sim_imm_tim_3$p < 0.05) | (sim_imm_tim_3$average > 0.1),]
  
  levels(as.factor(c(rownames(contrib_imm_tim_1), rownames(contrib_imm_tim_2), rownames(contrib_imm_tim_3))))
  
  
  sim_imm_season <- summary(simper(df_mean, meta_mean$immersion_season))
  sim_imm_recovery <- summary(simper(df_mean, meta_mean$recovery_season))
  
  
  return()  
    
}

  