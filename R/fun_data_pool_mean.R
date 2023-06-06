#' pool msp with mean data
#'
#' @param metadata_data_mean the path to the raw data file
#'
#' @return 
#' @export
#' 

fun_pool_mean <- function(metadata_data_mean){
  
  # metadata_data_mean = targets::tar_read(mean_metadata_data)
  
  df_mean <- read.csv(metadata_data_mean[!grepl("metadata", metadata_data_mean)], header = TRUE)
  meta_mean <- read.csv(metadata_data_mean[grepl("metadata", metadata_data_mean)], header = TRUE)
  
  spo_columns <- grep("spo", names(df_mean), value = TRUE)
  spo_columns <- spo_columns[-14]
  spo_mean <- rowSums(df_mean[,spo_columns])
  
  bryo_columns <- grep("bryo", names(df_mean), value = TRUE)
  bryo_mean <- rowSums(df_mean[,bryo_columns])
  
  ascc_columns <- grep("ascc", names(df_mean), value = TRUE)
  ascc_mean <- rowSums(df_mean[,ascc_columns])
  
  ascs_columns <- grep("ascs", names(df_mean), value = TRUE)
  ascs_mean <- rowSums(df_mean[,ascs_columns])
  
  for_columns <- grep("for", names(df_mean), value = TRUE)
  for_mean <- rowSums(df_mean[,for_columns])
  
  algae_columns <- grep("algae", names(df_mean), value = TRUE)
  algae_mean <- rowSums(df_mean[,algae_columns])
  
  worm_columns <- grep("worm", names(df_mean), value = TRUE)
  worm_mean <- rowSums(df_mean[,worm_columns])
  
  prokariot_columns <- c(grep("biofilm", names(df_mean), value = TRUE), grep("Cyanob", names(df_mean), value = TRUE))
  prokariot_mean <- rowSums(df_mean[,prokariot_columns])
  
  rest_columns <- names(df_mean)
  rest_columns <- rest_columns[c(1,3,17,55,56,58)]
  rest_data <- df_mean[,rest_columns]
  
  data_pool <- data.frame(porifera = spo_mean,
                          bryozoa = bryo_mean,
                          ascidiacea_c = ascc_mean,
                          ascidiacea_s = ascs_mean,
                          foraminifera = for_mean,
                          other_algae = algae_mean,
                          annelida = worm_mean,
                          prokariotic_biotas = prokariot_mean)
  
  data_pool <- data.frame(cbind(data_pool,rest_data))
  rownames(data_pool) <- meta_mean$arms
  
  data_pool_mean_path <- paste0("data/derived-data/data_pool_mean.csv") #Nom du fichier de data généré
  
 
  write.csv(data_pool, file = data_pool_mean_path, row.names = TRUE)
  
  return(data_pool_mean_path)

  
  
}