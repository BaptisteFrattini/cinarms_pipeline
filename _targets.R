library(targets)

tar_source()

list(
  tar_target(raw_data, "data/raw-data/Data_sans_UNAV(2).csv", format = "file")
  
  ,tar_target(campain_id, "CINA")
  
  ,tar_target(arms_id_2y, c("RUNA2A", "RUNA2B","RUNA2C"))
  
  ,tar_target(metadata_data, data_arms(raw_data = raw_data, 
                                       arms_id = campain_id,
                                       arms_id_2years = arms_id_2y), format = "file")
  
  ,tar_target(data_pool, fun_pool_full(meta_data = metadata_data))
  
  ,tar_target(mean_metadata_data, fun_data_mean(meta_data = metadata_data, 
                                                arms_id = campain_id), format = "file") 
  
  ,tar_target(data_pool_mean, fun_pool_mean(metadata_data_mean = mean_metadata_data))
  
  ,tar_target(nmds_plot, fun_nmds_plot(metadata_data_mean = mean_metadata_data))
  
  ,tar_target(permanova, fun_perm(metadata_data_mean = mean_metadata_data)) 
  
  ,tar_target(venn_plot, fun_Venn(metadata_data_mean = mean_metadata_data)) 
  
  ,tar_target(boxplot_pool, boxplot_explo(data_full_pool = data_pool, 
                                          meta_data = metadata_data))
  
  ,tar_target(div_explo, diversity_explo(metadata_data_mean = mean_metadata_data)) 
  
  ,tar_target(div_beta_decomp, beta_div_decomp(metadata_data_mean = mean_metadata_data)) 
  
  ,tar_target(div_beta_decomp_full, beta_div_decomp_full(metadata_data_mean = mean_metadata_data))
  
  ,tar_target(div_alpha, fun_alpha_div(metadata_data_mean = mean_metadata_data))
  
  ,tar_target(PCA_fun, fun_PCA(metadata_data_mean = mean_metadata_data,
                               data_mean_pool = data_pool_mean))
  
  )