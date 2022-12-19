#' Subset a raw_data table for a sampling campain
#'
#' @param meta_and_data the path to the metadata and data
#' @param dat_thresh_path the path to the file data red

#'
#' @return path to the data file, reduced and mean by site
#' @export
#'

mean_by_arms <- function(meta_and_data) {
   
  #meta_and_data = targets::tar_read("metadata_data") 
   
  meta_path <- meta_and_data[grepl("metadata", meta_and_data)] 
  meta <- read.csv(meta_path)
  arms_name <- meta$arms_name
  
  dat <- read.csv(meta_and_data[!grepl("metadata", meta_and_data)])
  
  tab <- NULL
  U <- NULL
  
  for (i in 1:ncol(dat)) {
    
    U <- tapply(dat[,i], 
                arms_name, 
                mean)
    
    tab <- cbind(tab,
                 U)
  }
  
  N <- colnames(dat)
  colnames(tab) <- N
  tab <- as.data.frame(tab)
  path_to_derived_data <- "data/derived-data"
  dat_mean_name <- "data_mean.csv"
  dat_mean_path <- here::here(path_to_derived_data, dat_mean_name)
  
  write.csv(tab, file = dat_mean_path, row.names = TRUE)
  
  return(dat_mean_path)
}
