#' species area curve to know if tree arms per site is enough
#'
#' @param metadata_data data
#' @return the path to the subseted raw data file
#' @export

fun_alpha_div <- function(metadata){
  
  # metadata = targets::tar_read(metadata_data)
  library(ggplot2)
  df_mean <- read.csv(metadata[!grepl("metadata", metadata)], header = TRUE)
  meta_mean <- read.csv(metadata[grepl("metadata", metadata)], header = TRUE)
  df_mean_pa <- vegan::decostand(df_mean, "pa")
 
  div_alpha_name <- paste0("div_alpha_total.pdf")
  div_alpha_path <- here::here("outputs/", div_alpha_name)
  pdf(file =  div_alpha_path, width = 7.65, height = 8.5)
  
  par(mfrow = c(3, 2),
      mar = c(4, 4, 1, 2))
  
  # #### Species abundance distribution ####
  # 
  # species_abundance <- colSums(df_mean)
  # sad_data <- data.frame(Species = names(species_abundance),
  #                        Abundance = species_abundance)
  # sad_data <- sad_data[order(-sad_data$Abundance), ]
  # 
  # sad_data <- sad_data[-c(2,5), ]
  # 
  # yp <- ggplot2::ggplot(sad_data, aes(x = reorder(Species, -Abundance), y = Abundance)) +
  #       geom_bar(stat = "identity") +
  #       labs(title = "Species Abundance Distribution",
  #            x = "Species",
  #            y = "Abundance") +
  #   theme_minimal() +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1))
  # 
  # 
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
       xlab="",
       ylab="# MSPs detected in all 15 ARMS",
       ylim=c(1,80))
  
  # boxplot(s, col="yellow", add=TRUE, pch="+")
  
  
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
       xlab="",
       ylab="# MSPs detected in 6-month ARMS",
       ylim=c(1,80))
  
  # boxplot(s, col="yellow", add=TRUE, pch="+")
  
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
       xlab="",
       ylab="# MSPs detected in one-year ARMS",
       ylim=c(1,80))
  
  # boxplot(s, col="yellow", add=TRUE, pch="+")
  
  
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
       xlab="",
       ylab="# MSPs detected in two-year ARMS",
       ylim=c(1,80))
  
  # boxplot(s, col="yellow", add=TRUE, pch="+")
  
  
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
       xlab="Number of plate faces analysed",
       ylab="# MSPs detected in ARMS deployed in hot season",
       ylim=c(1,80))
  
  # boxplot(s, col="yellow", add=TRUE, pch="+")
  
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
       xlab="Number of plate faces analysed",
       ylab="# MSPs detected in ARMS deployed in cool season",
       ylim=c(1,80))
  
  # boxplot(s, col="yellow", add=TRUE, pch="+")
  
  
  text(49,
       10,
       paste("S = ", pool$Species, " ; \n Estimation of species richness (Chao) = ", round(pool$chao,2),"±",round(pool$chao.se,2)),
       cex = 0.85)
  
  dev.off()  
  
  
  #### return ####
 return(div_alpha_path) 
}
