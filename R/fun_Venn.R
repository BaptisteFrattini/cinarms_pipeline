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
  
  #### Venn diag with immersion time ####
  
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
  rownames(my_table_with_zero) <- my_table_with_zero$msp
  my_table_with_zero <- my_table_with_zero[,-1]
  # Créer un vecteur de caractères avec les noms des MSP présentes dans la colonne "six_month"
  six_month_msp <- rownames(my_table_with_zero)[my_table_with_zero["six_month"] == 1]
  
  one_year_msp <- rownames(my_table_with_zero)[my_table_with_zero["one_year"] == 1]
  
  two_years_msp <- rownames(my_table_with_zero)[my_table_with_zero["two_years"] == 1]
  
  # Trouver les MSP présentes dans les deux colonnes (intersection)
  msp_six_one <- intersect(six_month_msp, one_year_msp)
  
  msp_six_two <- intersect(six_month_msp, two_years_msp)
  
  msp_one_two <- intersect(two_years_msp, one_year_msp)

  
  # récupérer les espèces présentes uniquement dans la colonne "six_month"
  only_six_month_msp <- setdiff(six_month_msp, union(one_year_msp, two_years_msp))
  
  only_one_year_msp <- setdiff(one_year_msp, union(six_month_msp, two_years_msp))
  
  only_two_years_msp <- setdiff(two_years_msp, union(six_month_msp, one_year_msp))
  
  
  
  
  
  library(ggvenn)
  six_month <- my_table$msp[my_table$six_month == 1]
  one_year <- my_table$msp[my_table$one_year == 1]
  two_years <- my_table$msp[my_table$two_year == 1]
  
  x <- list(Six_month = six_month,
            One_year = one_year,
            Two_years = two_years)
  Venn_imm_tim_path <- here::here("outputs/Venn_imm_tim.pdf")
  pdf(file = Venn_imm_tim_path, width = 7.5, height = 7.5)
  
  venn <- ggvenn(x, fill_color = c("limegreen", "maroon1", "navy"))
  
  
  venn
  
  dev.off()
  

  #### Venn diag with immersion season ####
  df_mean <- df_mean[1:12,]
  meta_mean <- meta_mean[1:12,]
  df_imsea <- aggregate(df_mean, list(meta_mean$immersion_season), mean)
  row.names(df_imsea) <- df_imsea[,1]
  df_imsea <- df_imsea[,-1]
  df_imsea_pa <- decostand(df_imsea, "pa")
  my_table <- data.frame(t(df_imsea_pa))
  
  my_table_with_zero <- subset(my_table, rowSums(my_table == 0) > 0)
  
  only_cool_msp <- rownames(my_table_with_zero)[my_table_with_zero["imm_cool"] == 1]
  only_hot_msp <- rownames(my_table_with_zero)[my_table_with_zero["imm_hot"] == 1]
  
  cool_msp <- rownames(my_table)[my_table["imm_cool"] == 1]
  hot_msp <- rownames(my_table)[my_table["imm_hot"] == 1]
  
  msp_cool_hot <- intersect(cool_msp, hot_msp)
  
  x <- list(Cool = cool_msp,
            Hot = hot_msp)
  Venn_immersion_season_path <- here::here("outputs/Venn_immersion_season.pdf")
  pdf(file = Venn_immersion_season_path, width = 5, height = 5)
  
  venn <- ggvenn(x, fill_color = c("royalblue1", "red3"))
  
  
  venn
  
  dev.off()
  
  #### Venn diag with recovery season ####
  
  df_recsea <- aggregate(df_mean, list(meta_mean$recovery_season), mean)
  row.names(df_recsea) <- df_recsea[,1]
  df_recsea <- df_recsea[,-1]
  df_recsea_pa <- decostand(df_recsea, "pa")
  my_table <- data.frame(t(df_recsea_pa))
  
  my_table_with_zero <- subset(my_table, rowSums(my_table == 0) > 0)
  
  only_cool_msp <- rownames(my_table_with_zero)[my_table_with_zero["rec_cool"] == 1]
  only_hot_msp <- rownames(my_table_with_zero)[my_table_with_zero["rec_hot"] == 1]
  
  cool_msp <- rownames(my_table)[my_table["rec_cool"] == 1]
  hot_msp <- rownames(my_table)[my_table["rec_hot"] == 1]
  
  msp_cool_hot <- intersect(cool_msp, hot_msp)
  
  x <- list(Cool = cool_msp,
            Hot = hot_msp)
  Venn_recovery_season_path <- here::here("outputs/Venn_recovery_season.pdf")
  pdf(file = Venn_recovery_season_path, width = 5, height = 5)
  
  venn <- ggvenn(x, fill_color = c("royalblue1", "red3"))
  
  
  venn
  
  dev.off()
  return(Venn_imm_tim_path)
   
}
  


