#' Mean all species cover of each plates of each arms 
#'
#' @param raw_data the path to the raw data file
#' @param arms_id the ID of the arms to subset for
#'
#' @return the path to the subseted derived data file
#' @export

fun_data_mean <- function(meta_data, arms_id){

  # meta_data = targets::tar_read(metadata_data)
  # arms_id = targets::tar_read(campain_id)
  library(dplyr)
  meta = read.csv(meta_data[grepl("metadata", meta_data)], header = TRUE)
  meta <- meta[,-4]
  data = read.csv(meta_data[!grepl("metadata", meta_data)], header = TRUE)
  data <- data[ , colSums(data) != 0]
  
  df_mean <- data %>% 
    group_by(meta$arms_name) %>% 
    summarise_all(mean, na.rm = TRUE)
  
  df_mean <- as.data.frame(df_mean)
  meta_mean <- data.frame(num = c(1:18),
                          arms = df_mean$`meta$arms_name`)
  
  
  meta_mean$imm_time <- c(rep("6m", 3), rep("1y", 3), rep("6m", 3), rep("1y", 3), rep("2y",6))
  meta_mean$arms_name <- substr(meta_mean$arms, 1, 5)
  hot_season_arms_imm <- c("CINA1", "CINA2", "RUNA2","RUNA3")
  hot_season_arms_rec <- c("RUNA2", "RUNA3", "CINA2", "CINA3")
  meta_mean$immersion_season <- ifelse(meta_mean$arms_name %in% hot_season_arms_imm, "imm_hot", "imm_cool")
  meta_mean$recovery_season <- ifelse(meta_mean$arms_name  %in% hot_season_arms_rec, "rec_hot", "rec_cool")
  
  
  rownames(df_mean) <- df_mean$`meta$arms_name`
  df_mean <- df_mean[,-1]
  
  

  data_out_f_path <- paste0("data/derived-data/data_mean_", arms_id, ".csv") #Nom du fichier de data généré
  meta_out_f_path <- paste0("data/derived-data/metadata_mean_", arms_id, ".csv") #Nom du fichier de 
  #metadata généré dans ce dossier
  
  write.csv(meta_mean, file = meta_out_f_path, row.names = FALSE)
  write.csv(df_mean, file = data_out_f_path, row.names = FALSE)
  
  res <- c(data_out_f_path, meta_out_f_path)
  names(res) <- c("path_data_mean", "path_meta_mean")
  
  return(res)
}
  
  
  
  
  
  
  
  