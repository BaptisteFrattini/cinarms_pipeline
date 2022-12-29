library(targets)

tar_source()

list(
  tar_target(raw_data, "data/raw-data/Data_sans_UNAV-NR-OROS.csv", format = "file")
  
  ,tar_target(campain_id, "RUNA")
  
  ,tar_target(metadata_data, data_arms(raw_data = raw_data, 
                                       arms_id = campain_id), format = "file")
  
  ,tar_target(arms_mean, mean_by_arms(meta_and_data = metadata_data))
  
  ,tar_target(b_div, decomp_b_div(meta_and_data = metadata_data, 
                                  mean_by_arms = arms_mean,
                                  arms_id = campain_id))
  ,tar_target(ns, north_south(meta_and_data = metadata_data, 
                                  mean_by_arms = arms_mean,
                                  arms_id = campain_id))
  
  )