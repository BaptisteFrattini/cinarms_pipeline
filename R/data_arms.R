#' Subset a raw_data table for a sampling campain
#'
#' @param raw_data the path to the raw data file
#' @param arms_id the ID of the arms to subset for
#'
#' @return the path to the subseted raw data file
#' @export

data_arms <- function(raw_data, arms_id, arms_id_2years){
  # raw_data = targets::tar_read(raw_data)
  # arms_id = targets::tar_read(campain_id)
  # arms_id_2years = targets::tar_read(arms_id_2y)

  dat_path <- here::here(raw_data)
  data <- read.table(dat_path, 
                     header = TRUE, 
                     sep = ";", 
                     dec = ",")
  
  data 
  
  dat <- data[data$prefixe == arms_id, ]
  dat2y <- data[data$arms_name == arms_id_2years, ]
  dat <- data.frame(rbind(dat, dat2y))
  
  
  meta_names <- as.vector(colnames(dat[,c(1:19)]))
  
  
  meta_data <- dat[, meta_names]
 
  
  library(dplyr)
  meta_data$arms = substr(meta_data$arms_name, 1, 5)
  hot_season_arms_imm <- c("CINA1", "CINA2", "RUNA2", "RUNA3")
  hot_season_arms_rec <- c("RUNA2", "CINA2", "CINA3", "RUNA3")
  meta_data$immersion_season <- ifelse(meta_data$arms %in% hot_season_arms_imm, "imm_hot", "imm_cold")
  meta_data$recovery_season <- ifelse(meta_data$arms  %in% hot_season_arms_rec, "rec_hot", "rec_cold")
  
  meta_data$imm_time <- ifelse(meta_data$incubation == 0.5, "6m",
                               ifelse(meta_data$incubation == 1.0, "1y",
                                      ifelse(meta_data$incubation == 2.0, "2y", NA)))
  
  meta_data <- meta_data[,-c(3,4,5,6,7,18)]
  
  dat <- dat[, !(names(dat) %in% meta_names)]
  
  # clean data for zero sum columns
  
  dat <- dat[ , colSums(dat) != 0]
  
  sums <- rowSums(dat)
  
  for(i in 1:nrow(dat)) {
    dat[i,] <- dat[i,]/sums[i]*100
  }
  

  
  out_d_path <- "data/derived-data" #Nom du chemin qui mène à data derived
  out_f_name <- paste0("data_", arms_id, "_",substr(arms_id_2years[1],1,5), ".csv") #Nom du fichier de data généré
  #dans ce dossier
  meta_out_f_name <- paste0("metadata_", arms_id, "_", substr(arms_id_2years[1],1,5), ".csv") #Nom du fichier de 
  #metadata généré dans ce dossier
  out_f_path <- here::here(out_d_path, out_f_name) #
  meta_out_f_path <- here::here(out_d_path, meta_out_f_name)
  write.csv(dat, file = out_f_path, row.names = FALSE)
  write.csv(meta_data, file = meta_out_f_path, row.names = FALSE)
  
  res <- c(out_f_path, meta_out_f_path)
  names(res) <- c("path_data", "path_meta")
  
  return(res)
  
}
