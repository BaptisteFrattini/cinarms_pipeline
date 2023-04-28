#' venn diagram
#'
#' @param metadata_data_mean the path to the raw data file
#'
#' @return 
#' @export
#' 

fun_Venn <- function(metadata_data_mean){
  
  # metadata_data_mean = targets::tar_read(mean_metadata_data)
  library(vegan)
  #### Load data and meta data ####
  
  df_mean <- read.csv(metadata_data_mean[!grepl("metadata", metadata_data_mean)], header = TRUE)
  meta_mean <- read.csv(metadata_data_mean[grepl("metadata", metadata_data_mean)], header = TRUE)
  
  df_tim <- aggregate(df_mean, list(meta_mean$imm_time), mean)
  row.names(df_tim) <- df_tim[,1]
  df_tim <- df_tim[,-1]
  df_tim_pa <- decostand(df_tim, "pa")
  df_tim_pa <- data.frame(t(df_tim_pa))
  df_tim_pa <- data.frame(msp = rownames(df_tim_pa),
                          six_month = df_tim_pa$X6m,
                          one_year = df_tim_pa$X1y,
                          two_years = df_tim_pa$X2y)
  
  
  my_table <- df_tim_pa[order(df_tim_pa$six_month, decreasing = TRUE), ]
  my_table_with_zero <- subset(my_table, rowSums(my_table == 0) > 0)
  
  library(ggvenn)
  six_month <- my_table$msp[my_table$six_month == 1]
  one_year <- my_table$msp[my_table$one_year == 1]
  two_years <- my_table$msp[my_table$two_year == 1]
  
  x <- list(Six_month = six_month,
            One_year = one_year,
            Two_years = two_years)
  Venn_imm_tim_path <- here::here("outputs/Venn_imm_tim.pdf")
  pdf(file = Venn_imm_tim_path, width = 15, height = 15)
  
  ggvenn(x, fill_color = c("limegreen", "maroon1", "navy"))
 
  
  dev.off()
  
  
  return(Venn_imm_tim_path)
   
}
  


