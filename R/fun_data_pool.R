#' pool msp with full data
#'
#' @param metadata_data_mean the path to the raw data file
#'
#' @return 
#' @export
#' 

fun_pool_full <- function(meta_data){
  
  # meta_data = targets::tar_read(metadata_data)
  
  data <- read.csv(meta_data[!grepl("metadata", meta_data)], header = TRUE)
  meta <- read.csv(meta_data[grepl("metadata", meta_data)], header = TRUE)
  
  spo_columns <- grep("spo", names(data), value = TRUE)
  spo_columns <- spo_columns[-14]
  spo_mean <- rowSums(data[,spo_columns])
  
  bryo_columns <- grep("bryo", names(data), value = TRUE)
  bryo_mean <- rowSums(data[,bryo_columns])
  
  
  ascc_columns <- grep("ascc", names(data), value = TRUE)
  ascc_mean <- rowSums(data[,ascc_columns])
  
  ascs_columns <- grep("ascs", names(data), value = TRUE)
  ascs_mean <- rowSums(data[,ascs_columns])
  
  for_columns <- grep("for", names(data), value = TRUE)
  for_mean <- rowSums(data[,for_columns])
  
  algae_columns <- grep("algae", names(data), value = TRUE)
  algae_mean <- rowSums(data[,algae_columns])
  
  worm_columns <- grep("worm", names(data), value = TRUE)
  worm_mean <- rowSums(data[,worm_columns])
  
  prokariot_columns <- c(grep("biofilm", names(data), value = TRUE), grep("Cyanob", names(data), value = TRUE))
  prokariot_mean <- rowSums(data[,prokariot_columns])
  
  rest_columns <- names(data)
  rest_columns <- rest_columns[c(1,3,17,55,56,58)]
  rest_data <- data[,rest_columns]
  
  data_pool <- data.frame(porifera = spo_mean,
                          bryozoa = bryo_mean,
                          ascidiacea_c = ascc_mean,
                          ascidiacea_s = ascs_mean,
                          foraminifera = for_mean,
                          other_algae = algae_mean,
                          annelida = worm_mean,
                          prokariotic_biotas = prokariot_mean)
  
  data_pool <- data.frame(cbind(data_pool,rest_data))
  
  data_pool_path <- paste0("data/derived-data/data_pool.csv") #Nom du fichier de data généré
  
  
  write.csv(data_pool, file = data_pool_path, row.names = FALSE)
  
  return(data_pool_path)
  
  
}