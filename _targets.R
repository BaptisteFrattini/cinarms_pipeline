library(targets)

targets::tar_source()

list(
  tar_target(raw_data, "data/raw-data/Data_sans_UNAV(2).csv", format = "file")
  
  ,tar_target(campain_id, "CINA")
  
  ,tar_target(arms_id_2y, c("RUNA2A", "RUNA2B","RUNA2C"))
  
  ,tar_target(metadata_data, data_arms(raw_data = raw_data, 
                                       arms_id = campain_id,
                                       arms_id_2years = arms_id_2y,
                                       micro_habitat = micro_hab), format = "file")
  
  ,tar_target(data_pool, fun_pool_full(meta_data = metadata_data))
  
  ,tar_target(mean_metadata_data, fun_data_mean(meta_data = metadata_data, 
                                                arms_id = campain_id), format = "file") 
  
  ,tar_target(data_pool_mean, fun_pool_mean(metadata_data_mean = mean_metadata_data))
  
  ,tar_target(nmds_plot, fun_nmds_plot(metadata_data_mean = mean_metadata_data))
  
  ,tar_target(permanova, fun_perm(metadata_data_mean = mean_metadata_data)) 
  
  ,tar_target(venn_plot, fun_Venn(metadata_data_mean = mean_metadata_data)) 
  
  ,tar_target(boxplot_pool_alt, boxplot_explo_alt(data_full_pool = data_pool, 
                                                  meta_data = metadata_data))
  
  ,tar_target(SAD, fun_SAD(metadata_data_mean = mean_metadata_data))
  
  ,tar_target(div_explo, diversity_explo(metadata_data_mean = mean_metadata_data)) 
  
  ,tar_target(div_beta_decomp, beta_div_decomp(metadata_data_mean = mean_metadata_data)) 
  
  ,tar_target(div_beta_decomp_full, beta_div_decomp_full(metadata_data_mean = mean_metadata_data))
  
  ,tar_target(div_beta_microhab, fun_beta_microhab(meta_data = metadata_data))
  
  ,tar_target(PCoA_microhab, fun_PCoA_microhab(meta_data = metadata_data))
  
  ,tar_target(div_alpha, fun_alpha_div(metadata = metadata_data))
  
  ,tar_target(PCA, fun_PCA(metadata_data_mean = mean_metadata_data,
                               data_mean_pool = data_pool_mean))
  
  ,tar_target(tab_div, fun_tab(metadata_data_mean = mean_metadata_data))
  
  ,tar_target(stack_c_chart, fun_stacked_col_chart(metadata_data_mean = mean_metadata_data,
                                                   data_mean_pool = data_pool_mean))
  
  ,tar_target(null_model, fun_null_model(metadata_data_mean = mean_metadata_data))

  ,tar_target(pielou, fun_pielou(metadata = metadata_data,
                                 metadata_data_mean = mean_metadata_data))  
  
  )