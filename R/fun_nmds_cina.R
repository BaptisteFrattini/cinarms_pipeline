#' NMDS
#'
#' @param metadata_data_mean the path to the raw data file
#'
#' @return the path to the NMDS plot
#' @export

fun_nmds_plot <- function(metadata_data_mean){
  
  # metadata_data_mean = targets::tar_read(mean_metadata_data)
  
  #### Load data and meta data ####
  
  df_mean <- read.csv(metadata_data_mean[!grepl("metadata", metadata_data_mean)], header = TRUE)
  meta_mean <- read.csv(metadata_data_mean[grepl("metadata", metadata_data_mean)], header = TRUE)

  #### Plot NMDS ####
  library(vegan)
  
  #IMMERSION SEASON
  NMDS_imm_path <- here::here("outputs/NMDS/NMDS_imm.pdf")
  pdf(file = NMDS_imm_path, width = 10, height = 10)
  
  ord <- metaMDS(df_mean, distance = "bray")
  sppscores(ord) <- df_mean
  plot(ord, main="NMDS (dist=bray)")
  imm <- meta_mean$immersion_season
  orditorp(ord,display="species",col="red",air=0.01)
  ordihull(ord, imm, col=1:2, lwd=3)
  ordiellipse(ord, imm, col=1:2, kind = "ehull", lwd=3)
  ordiellipse(ord, imm, col=1:2, draw="polygon")
  ordispider(ord,  imm, col=1:2, label = TRUE)
  dev.off()
  
  #RECOVERY SEASON
  # plot(ord, main="NMDS (dist=bray)")
  NMDS_rec_path <- here::here("outputs/NMDS/NMDS_rec.pdf")
  pdf(file = NMDS_rec_path, width = 10, height = 10)
  
  ord <- metaMDS(df_mean, distance = "bray")
  plot(ord, main="NMDS (dist=bray)")
  sppscores(ord) <- df_mean
  orditorp(ord,display="species",col="red",air=0.01)
  rec <- meta_mean$recovery_season
  ordihull(ord, rec, col=3:4, lwd=3)
  ordiellipse(ord, rec, col=3:4, kind = "ehull", lwd=3)
  ordiellipse(ord, rec, col=3:4, draw="polygon")
  ordispider(ord,  rec, col=3:4, label = TRUE)
  dev.off()
  
  #IMMERSION TIME
  NMDS_tim_path <- here::here("outputs/NMDS/NMDS_tim.pdf")
  pdf(file = NMDS_tim_path, width = 10, height = 10)
  sppscores(ord) <- df_mean
  plot(ord, main="NMDS (dist=bray)")
  tim <- meta_mean$imm_time
  orditorp(ord,display="species",col="red",air=0.01)
  ordihull(ord, tim, col=5:7, lwd=3)
  ordiellipse(ord, tim, col=5:7, kind = "ehull", lwd=3)
  ordiellipse(ord, tim, col=5:7, draw="polygon")
  ordispider(ord,  tim, col=5:7, label = TRUE)
  #text(ord, display = "spec", cex=0.7, col="blue")
  
  dev.off()
  
  #### Jaccard dissimilarity ####
  #IMMERSION SEASON
  NMDS_imm_path_jac <- here::here("outputs/NMDS/NMDS_imm_jac.pdf")
  pdf(file = NMDS_imm_path_jac, width = 10, height = 10)
 
  ord <- metaMDS(df_mean, distance = "jaccard")
  plot(ord, main="NMDS (dist=jacc)")
  imm <- meta_mean$immersion_season
  ordihull(ord, imm, col=1:2, lwd=3)
  ordiellipse(ord, imm, col=1:2, kind = "ehull", lwd=3)
  ordiellipse(ord, imm, col=1:2, draw="polygon")
  ordispider(ord,  imm, col=1:2, label = TRUE)
  dev.off()
  
  #RECOVERY SEASON
  # plot(ord, main="NMDS (dist=bray)")
  NMDS_rec_path_jac <- here::here("outputs/NMDS/NMDS_rec_jac.pdf")
  pdf(file = NMDS_rec_path_jac, width = 10, height = 10)
  
  ord <- metaMDS(df_mean, distance = "jaccard")
  plot(ord, main="NMDS (dist=jacc)")
  rec <- meta_mean$recovery_season
  ordihull(ord, rec, col=3:4, lwd=3)
  ordiellipse(ord, rec, col=3:4, kind = "ehull", lwd=3)
  ordiellipse(ord, rec, col=3:4, draw="polygon")
  ordispider(ord,  rec, col=3:4, label = TRUE)
  dev.off()
  
  #IMMERSION TIME
  NMDS_tim_path_jac <- here::here("outputs/NMDS/NMDS_tim_jac.pdf")
  pdf(file = NMDS_tim_path_jac, width = 10, height = 10)
  
  ord <- metaMDS(df_mean, distance = "jaccard")
  plot(ord, main="NMDS (dist=jacc)")
  tim <- meta_mean$imm_time
  ordihull(ord, tim, col=5:7, lwd=3)
  ordiellipse(ord, tim, col=5:7, kind = "ehull", lwd=3)
  ordiellipse(ord, tim, col=5:7, draw="polygon")
  ordispider(ord,  tim, col=5:7, label = TRUE)
  text(ord, display = "spec", cex=0.7, col="blue")
  dev.off()
  
  return(NMDS_tim_path)
}
