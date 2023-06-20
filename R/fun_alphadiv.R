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
 
  div_alpha_name <- paste0("div_alpha_total.pdf")
  div_alpha_path <- here::here("outputs/", div_alpha_name)
  pdf(file =  div_alpha_path, width = 12, height = 7)
  
  par(mfrow = c(2, 3))
  #### All ARMS pooled ####
  
  s <- vegan::specaccum(df_mean_pa, method = "random", permutations = 999,
                        conditioned =TRUE) 
  pool <- vegan::specpool(df_mean_pa)
  
  plot(s, 
       ci.type="poly",
       col="blue", 
       lwd=2,
       ci.lty=0,
       ci.col="lightblue",
       xlab="Number of plates analysed",
       ylab="# morpho-species detected in all 15 ARMS",
       ylim=c(1,80))
  
  boxplot(s, col="yellow", add=TRUE, pch="+")
  
  
  text(120,
       10,
       paste("S = ", pool$Species, " ; \n Estimation of species richness (Chao) = ", round(pool$chao,2),"±",round(pool$chao.se,2)),
       cex = 0.85)
  

  
  #### Only ARMS of 6 months ####
  
  df_mean_pa_six <- subset(df_mean_pa, meta_mean$imm_time == "6m")
  
  s <- vegan::specaccum(df_mean_pa_six, method = "random", permutations = 999,
                        conditioned =TRUE) 
  s
  pool <- vegan::specpool(df_mean_pa_six)
  
  
  plot(s, 
       ci.type="poly",
       col="blue", 
       lwd=2,
       ci.lty=0,
       ci.col="lightblue",
       xlab="Number of plates analysed",
       ylab="# morpho-species detected in 6-month ARMS",
       ylim=c(1,80))
  
  boxplot(s, col="yellow", add=TRUE, pch="+")
  
  text(50,
       10,
       paste("S = ", pool$Species, " ; \n Estimation of species richness (Chao) = ", round(pool$chao,2),"±",round(pool$chao.se,2)),
       cex = 0.85)
  
  
  
  #### Only ARMS of 1 year ####
  df_mean_pa_one <- subset(df_mean_pa, meta_mean$imm_time == "1y")
  
  s <- vegan::specaccum(df_mean_pa_one, method = "random", permutations = 999,
                        conditioned =TRUE) 
  s
  pool <- vegan::specpool(df_mean_pa_one)
 
  
  plot(s, 
       ci.type="poly",
       col="blue", 
       lwd=2,
       ci.lty=0,
       ci.col="lightblue",
       xlab="Number of plates analysed",
       ylab="# morpho-species detected in one-year ARMS",
       ylim=c(1,80))
  
  boxplot(s, col="yellow", add=TRUE, pch="+")
  
  
  text(51,
       10,
       paste("S = ", pool$Species, " ; \n Estimation of species richness (Chao) = ", round(pool$chao,2),"±",round(pool$chao.se,2)),
       cex = 0.85)
  
  
  #### Only ARMS of 2 years ####
  df_mean_pa_two <- subset(df_mean_pa, meta_mean$imm_time == "2y")
  
  s <- vegan::specaccum(df_mean_pa_two, method = "random", permutations = 999,
                        conditioned =TRUE) 
  s
  pool <- vegan::specpool(df_mean_pa_two)
  

  
  plot(s, 
       ci.type="poly",
       col="blue", 
       lwd=2,
       ci.lty=0,
       ci.col="lightblue",
       xlab="Number of plates analysed",
       ylab="# morpho-species detected in two-year ARMS",
       ylim=c(1,80))
  
  boxplot(s, col="yellow", add=TRUE, pch="+")
  
  
  text(25,
       6,
       paste("S = ", pool$Species, " ; \n Estimation of species richness (Chao) = ", round(pool$chao,2),"±",round(pool$chao.se,2)),
       cex = 0.85)
  

  
  #### Only ARMS deployed in hot season ####
  df_mean_pa_hot <- subset(df_mean_pa, meta_mean$immersion_season == "imm_hot")
  
  
  s <- vegan::specaccum(df_mean_pa_hot, method = "random", permutations = 999,
                        conditioned =TRUE) 
  s
  pool <- vegan::specpool(df_mean_pa_hot)
  
  
  
  plot(s, 
       ci.type="poly",
       col="blue", 
       lwd=2,
       ci.lty=0,
       ci.col="lightblue",
       xlab="Number of plates analysed",
       ylab="# morpho-species detected in ARMS deployed in hot season",
       ylim=c(1,80))
  
  boxplot(s, col="yellow", add=TRUE, pch="+")
  
  text(70,
       10,
       paste("S = ", pool$Species, " ; \n Estimation of species richness (Chao) = ", round(pool$chao,2),"±",round(pool$chao.se,2)),
       cex = 0.85)
  
 
  
  #### Only ARMS deployed in cool season ####
  df_mean_pa_cool <- subset(df_mean_pa, meta_mean$immersion_season == "imm_cold")
  
  
  s <- vegan::specaccum(df_mean_pa_cool, method = "random", permutations = 999,
                        conditioned =TRUE) 
  pool <- vegan::specpool(df_mean_pa_cool)
  

  
  plot(s, 
       ci.type="poly",
       col="blue", 
       lwd=2,
       ci.lty=0,
       ci.col="lightblue",
       xlab="Number of plates analysed",
       ylab="# morpho-species detected in ARMS deployed in cool season",
       ylim=c(1,80))
  
  boxplot(s, col="yellow", add=TRUE, pch="+")
  
  
  text(49,
       10,
       paste("S = ", pool$Species, " ; \n Estimation of species richness (Chao) = ", round(pool$chao,2),"±",round(pool$chao.se,2)),
       cex = 0.85)
  
  dev.off()  
  
  
  #### return ####
 return(div_alpha_path) 
}
