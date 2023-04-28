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
                         imm_time = meta_mean$imm_time)
  
  tapply(S, meta_mean$imm_time, mean)
  time <- c("6 month", "1 year", "2 years")
  library(ggpubr)
  library(forcats)
  library(ggplot2)
  ff <- ggplot(data_div, aes(x = fct_relevel(imm_time, "6m", "1y", "2y"), y = S)) +
    geom_boxplot(fill =  c("cadetblue2","cadetblue3","cadetblue4")) +
    labs(title = "",
         x = "Immersion time",
         y = "Average species richness in an ARMS") 
  path_to_box_div <- paste0("outputs/box_div.pdf")
  ggsave(filename =  path_to_box_div , width = 6, height = 6)
  
  #### What are the species common to all stage of colonisation ? ####
  df_pa <- vegan::decostand(df_mean, "pa")
  
  # Trouver les noms des colonnes qui ne présentent pas de 0
  names(df_pa)[colSums(df_pa == 0) == 0]
  
  
  
  return(path_to_box_div)
      
}  
