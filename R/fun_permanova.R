#' permanova
#'
#' @param metadata_data_mean the path to the raw data file
#'
#' @return 
#' @export
#' 


fun_permanova_plot <- function(metadata_data_mean){
  
  # metadata_data_mean = targets::tar_read(mean_metadata_data)
  
  #### Load data and meta data ####
  
  df_mean <- read.csv(metadata_data_mean[!grepl("metadata", metadata_data_mean)], header = TRUE)
  meta_mean <- read.csv(metadata_data_mean[grepl("metadata", metadata_data_mean)], header = TRUE)
  
  #### Compute permanov ####
  
  library(vegan)
  
  perm <- vegan::adonis2(df_mean ~ imm_time*immersion_season+recovery_season, data = meta_mean, method = "bray")
  
  pairwise_ado_ <- pairwise.adonis2(df_mean ~ imm_time, data = meta_mean, method = "bray")
  