library(targets)

tar_source()

list(
  tar_target(raw_data, "data/raw-data/Data_sans_UNAV(2).csv", format = "file")
  
  ,tar_target(campain_id, "CINA")
  
  ,tar_target(arms_id_2y, c("RUNA2A", "RUNA2B","RUNA2C"))
  
  ,tar_target(metadata_data, data_arms(raw_data = raw_data, 
                                       arms_id = campain_id,
                                       arms_id_2years = arms_id_2y), format = "file")
  
  ,tar_target(mean_metadata_data, fun_data_mean(meta_data = metadata_data, 
                                                arms_id = campain_id), format = "file") 
  
  ,tar_target(nmds_plot, fun_nmds_plot(metadata_data_mean = mean_metadata_data))

  )