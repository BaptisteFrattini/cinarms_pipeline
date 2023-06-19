#' Tableau de shannon et pielou
#'
#' @param metadata_data_mean data
#' @return the path to the subseted raw data file
#' @export

fun_tab <- function(metadata_data_mean){
  
  # metadata_data_mean = targets::tar_read(mean_metadata_data)
  df_mean <- read.csv(metadata_data_mean[!grepl("metadata", metadata_data_mean)], header = TRUE)
  meta_mean <- read.csv(metadata_data_mean[grepl("metadata", metadata_data_mean)], header = TRUE)  
  
  # df_mean_pa <- vegan::decostand(df_mean, "pa")
  
  vegan::diversity(df_mean)
  vegan::specnumber(df_mean)
  
  piel <- vegan::diversity(df_mean, index = "shannon") / log(vegan::specnumber(df_mean))
  
  
  div <- data.frame(meta_mean$arms, vegan::specnumber(df_mean), vegan::diversity(df_mean), piel)
  
  write.csv(div, file = "outputs/tabdiv.csv")
  
  return()
  
}