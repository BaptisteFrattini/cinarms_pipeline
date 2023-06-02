#' Boost regression trees to know if tree arms per site is enough
#'
#' @param metadata_data_mean data
#' @return the path to the subseted raw data file
#' @export

fun_alpha_div <- function(metadata_data_mean){
  
  # metadata_data_mean = targets::tar_read(mean_metadata_data)
  
  df_mean <- read.csv(metadata_data_mean[!grepl("metadata", metadata_data_mean)], header = TRUE)
  meta_mean <- read.csv(metadata_data_mean[grepl("metadata", metadata_data_mean)], header = TRUE)
  
  df_mean_pa <- vegan::decostand(df_mean, "pa")
  
  s <- vegan::specaccum(df_mean_pa, method = "random", permutations = 999,
                        conditioned =TRUE) 
  s
  pool <- vegan::specpool(df_mean_pa)
  
  div_alpha_name <- paste0("div_alpha.pdf")
  div_alpha_path <- here::here("outputs/", div_alpha_name)
  pdf(file =  div_alpha_path, width = 6, height = 6)
  
  plot(s, 
       ci.type="poly",
       col="blue", 
       lwd=2,
       ci.lty=0,
       ci.col="lightblue",
       xlab="Number of ARMS units",
       ylab="Average number of morpho-species detected")
  
  boxplot(s, col="yellow", add=TRUE, pch="+")
  
  
  text(8.5,
       10,
       paste("S = ", pool$Species, " ; \n Estimation of species richness (Chao) = ", round(pool$chao,2),"±",round(pool$chao.se,2)),
       cex = 0.85)
  
  dev.off()  
 return(div_alpha_path) 
}