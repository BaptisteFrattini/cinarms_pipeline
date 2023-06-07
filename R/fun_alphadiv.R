#' species area curve to know if tree arms per site is enough
#'
#' @param metadata_data data
#' @return the path to the subseted raw data file
#' @export

fun_alpha_div <- function(metadata){
  
  # metadata = targets::tar_read(metadata_data)
  
  df_mean <- read.csv(metadata[!grepl("metadata", metadata)], header = TRUE)
  meta_mean <- read.csv(metadata[grepl("metadata", metadata)], header = TRUE)
  df_mean_pa <- vegan::decostand(df_mean, "pa")
 
  #### All ARMS pooled ####
  
  s <- vegan::specaccum(df_mean_pa, method = "random", permutations = 999,
                        conditioned =TRUE) 
  pool <- vegan::specpool(df_mean_pa)
  
  div_alpha_name <- paste0("div_alpha_all.pdf")
  div_alpha_path <- here::here("outputs/", div_alpha_name)
  pdf(file =  div_alpha_path, width = 8, height = 6)
  
  plot(s, 
       ci.type="poly",
       col="blue", 
       lwd=2,
       ci.lty=0,
       ci.col="lightblue",
       xlab="Number of plates analysed",
       ylab="# morpho-species detected in all 15 ARMS")
  
  boxplot(s, col="yellow", add=TRUE, pch="+")
  
  
  text(150,
       10,
       paste("S = ", pool$Species, " ; \n Estimation of species richness (Chao) = ", round(pool$chao,2),"±",round(pool$chao.se,2)),
       cex = 0.85)
  
  dev.off()
  
  #### Only ARMS of 6 months ####
  
  df_mean_pa_six <- subset(df_mean_pa, meta_mean$imm_time == "6m")
  
  s <- vegan::specaccum(df_mean_pa_six, method = "random", permutations = 999,
                        conditioned =TRUE) 
  s
  pool <- vegan::specpool(df_mean_pa_six)
  
  div_alpha_name <- paste0("div_alpha_six.pdf")
  div_alpha_path <- here::here("outputs/", div_alpha_name)
  pdf(file =  div_alpha_path, width = 8, height = 6)
  
  plot(s, 
       ci.type="poly",
       col="blue", 
       lwd=2,
       ci.lty=0,
       ci.col="lightblue",
       xlab="Number of plates analysed",
       ylab="# morpho-species detected in 6-month ARMS")
  
  boxplot(s, col="yellow", add=TRUE, pch="+")
  
  text(60,
       10,
       paste("S = ", pool$Species, " ; \n Estimation of species richness (Chao) = ", round(pool$chao,2),"±",round(pool$chao.se,2)),
       cex = 0.85)
  
  dev.off()
  
  #### Only ARMS of 1 year ####
  df_mean_pa_one <- subset(df_mean_pa, meta_mean$imm_time == "1y")
  
  s <- vegan::specaccum(df_mean_pa_one, method = "random", permutations = 999,
                        conditioned =TRUE) 
  s
  pool <- vegan::specpool(df_mean_pa_one)
  
  div_alpha_name <- paste0("div_alpha_one.pdf")
  div_alpha_path <- here::here("outputs/", div_alpha_name)
  pdf(file =  div_alpha_path, width = 8, height = 6)
  
  plot(s, 
       ci.type="poly",
       col="blue", 
       lwd=2,
       ci.lty=0,
       ci.col="lightblue",
       xlab="Number of plates analysed",
       ylab="# morpho-species detected in one-year ARMS")
  
  boxplot(s, col="yellow", add=TRUE, pch="+")
  
  
  text(45,
       10,
       paste("S = ", pool$Species, " ; \n Estimation of species richness (Chao) = ", round(pool$chao,2),"±",round(pool$chao.se,2)),
       cex = 0.85)
  
  dev.off()
  
  #### Only ARMS of 2 years ####
  df_mean_pa_two <- subset(df_mean_pa, meta_mean$imm_time == "2y")
  
  s <- vegan::specaccum(df_mean_pa_two, method = "random", permutations = 999,
                        conditioned =TRUE) 
  s
  pool <- vegan::specpool(df_mean_pa_two)
  
  div_alpha_name <- paste0("div_alpha_two.pdf")
  div_alpha_path <- here::here("outputs/", div_alpha_name)
  pdf(file =  div_alpha_path, width = 8, height = 6)
  
  plot(s, 
       ci.type="poly",
       col="blue", 
       lwd=2,
       ci.lty=0,
       ci.col="lightblue",
       xlab="Number of plates analysed",
       ylab="# morpho-species detected in two-year ARMS")
  
  boxplot(s, col="yellow", add=TRUE, pch="+")
  
  
  text(10,
       6,
       paste("S = ", pool$Species, " ; \n Estimation of species richness (Chao) = ", round(pool$chao,2),"±",round(pool$chao.se,2)),
       cex = 0.85)
  
  dev.off()
  
  #### Only ARMS deployed in hot season ####
  df_mean_pa_hot <- subset(df_mean_pa, meta_mean$immersion_season == "imm_hot")
  
  
  s <- vegan::specaccum(df_mean_pa_hot, method = "random", permutations = 999,
                        conditioned =TRUE) 
  s
  pool <- vegan::specpool(df_mean_pa_hot)
  
  div_alpha_name <- paste0("div_alpha_hot.pdf")
  div_alpha_path <- here::here("outputs/", div_alpha_name)
  pdf(file =  div_alpha_path, width = 8, height = 6)
  
  plot(s, 
       ci.type="poly",
       col="blue", 
       lwd=2,
       ci.lty=0,
       ci.col="lightblue",
       xlab="Number of plates analysed",
       ylab="# morpho-species detected in ARMS deployed in hot season")
  
  boxplot(s, col="yellow", add=TRUE, pch="+")
  
  text(60,
       10,
       paste("S = ", pool$Species, " ; \n Estimation of species richness (Chao) = ", round(pool$chao,2),"±",round(pool$chao.se,2)),
       cex = 0.85)
  
  dev.off()
  
  #### Only ARMS deployed in cool season ####
  df_mean_pa_cool <- subset(df_mean_pa, meta_mean$immersion_season == "imm_cold")
  
  
  s <- vegan::specaccum(df_mean_pa_cool, method = "random", permutations = 999,
                        conditioned =TRUE) 
  pool <- vegan::specpool(df_mean_pa_cool)
  
  div_alpha_name <- paste0("div_alpha_cool.pdf")
  div_alpha_path <- here::here("outputs/", div_alpha_name)
  pdf(file =  div_alpha_path, width = 8, height = 6)
  
  plot(s, 
       ci.type="poly",
       col="blue", 
       lwd=2,
       ci.lty=0,
       ci.col="lightblue",
       xlab="Number of plates analysed",
       ylab="# morpho-species detected in ARMS deployed in cool season")
  
  boxplot(s, col="yellow", add=TRUE, pch="+")
  
  
  text(65,
       10,
       paste("S = ", pool$Species, " ; \n Estimation of species richness (Chao) = ", round(pool$chao,2),"±",round(pool$chao.se,2)),
       cex = 0.85)
  
  dev.off()  
  
  #### return ####
 return(div_alpha_path) 
}
