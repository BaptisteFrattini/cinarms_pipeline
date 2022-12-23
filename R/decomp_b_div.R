#' Subset a raw_data table for a sampling campain
#'
#' @param meta_and_data the path to the metadata and data
#' @param dat_thresh_path the path to the file data red
#'https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12029
#' @return path to the data file, reduced and mean by site
#' @export
#'

decomp_b_div <- function(meta_and_data, mean_by_arms) {
  #mean_by_arms = targets::tar_read("arms_mean") 
  #meta_and_data = targets::tar_read("metadata_data") 
  
  #### data load ####
  data <- read.csv(mean_by_arms)
  rownames(data) <- data[,1]
  data <- data[,-1]
  matrix.hel <- vegan::decostand(data, "hellinger")
  matrix.pa <- vegan::decostand(data, "pa")

  #### Bray curtis & Jaccard computing ####
  betapart::beta.multi.abund(matrix.hel)
  B.pair.abund <- betapart::beta.pair.abund(matrix.hel)
  B.pair.pa <- betapart::beta.pair(matrix.pa, index.family = "jaccard")
  
  #### triangle plot ####
  
  
  decomp.ab <- cbind((1 - B.pair.abund$beta.bray), 
                     B.pair.abund$beta.bray.bal,
                     B.pair.abund$beta.bray.gra)
  colnames(decomp.ab) <- c("1 - BC Similarity", "Balanced variation", "abundance gradiant") 
  mean( B.pair.abund$beta.bray.bal)
  
  triangle_name <- paste0("triangle_", arms_id, ".pdf")
  triangle_path <- here::here("outputs",  triangle_name)
  pdf(file =  triangle_path, width = 11, height = 7)
  par(mfrow = c(1,2))
  
  ade4::triangle.plot(as.data.frame(decomp.ab), 
                      labeltriangle = FALSE, 
                      addmean = TRUE,
                      show.position = TRUE,
                      scale = TRUE,
                      cpoint = 0.1)
  
  text(0.61, 0.2, "1 - Bray-Curtis \n dissimilarity", cex  = 1)
  text(0.1, -0.6, "Balanced variation", cex = 1)
  text(-0.57, 0.2, "Abundance \n gradiant", cex = 1) 
  
  
  
  decomp.pa <-  cbind((1 - B.pair.pa$beta.jac), 
                           B.pair.pa$beta.jne,
                           B.pair.pa$beta.jtu)
  
  colnames(decomp.ab) <- c("1 - Jaccard diss", "Nestedness", "Turnover") 
  
  ade4::triangle.plot(as.data.frame(decomp.pa), 
                      labeltriangle = FALSE, 
                      addmean = TRUE,
                      show.position = TRUE,
                      scale = TRUE,
                      cpoint = 0.1)
  
  
  text(0.58, 0.12, "Turnover", cex  = 1)
  text(0.1, -0.6, "Nestedness", cex = 1)
  text(-0.57, 0.2, "1 - Jaccard \n dissimilarity", cex = 1) 
  
  
  #### Mantel with abundance based dissimilarity decomposition ####
  meta_path <- meta_and_data[grepl("metadata", meta_and_data)] 
  meta <- read.csv(meta_path)

  arms_name <- meta$arms_name
  latitude <- as.numeric(meta$latitude)
  longitude <- as.numeric(meta$longitude)
  tab <- as.data.frame(cbind(arms_name, latitude, longitude))
  lat.pool <- tapply(latitude, arms_name, mean)
  long.pool <- tapply(longitude, arms_name, mean)
  tab <- as.data.frame(cbind(long.pool,lat.pool,levels(arms_name)))
  
  
  a <- geosphere::distGeo(as.numeric(tab[1,]),as.numeric(tab[5,])) 
  
  matrix.dist <- geosphere::distm(tab)
  row.names(matrix.dist) <- row.names(tab)
  colnames(matrix.dist) <- colnames(matrix.dist)
  matrix.dist <- matrix.dist/1000
  matrix.dist <- as.dist(matrix.dist)
  matrix.dist <- round(matrix.dist, 3)
  
  matrix.bal <- B.pair.abund$beta.bray.bal
  matrix.gra <- B.pair.abund$beta.bray.gra
  matrix.bray <- B.pair.abund$beta.bray
  matrix.turn <- B.pair.pa$beta.jtu
  matrix.nest <- B.pair.pa$beta.jne
  matrix.jacc <- B.pair.pa$beta.jac
  
  
  
  
  correl
  
  #### plot mantel ####
  
  triangle_name <- paste0("triangle_", arms_id, ".pdf")
  triangle_path <- here::here("outputs",  triangle_name)
  pdf(file =  triangle_path, width = 10, height = 7)

  library(ggplot2)  
  #correl = vegan::mantel(matrix.dist, matrix.jacc, method = "pearson")
  ##### jacc #####
  aa = as.vector(matrix.jacc)
  tt = as.vector(matrix.dist)
  #new data frame with vectorized distance matrices
  mat = data.frame(aa,tt)
  
  mm1 = ggplot(mat, aes(y = aa, x = tt)) + 
    geom_point(size = 3, alpha = 0.5, color = "black") + 
    labs(x = "Geographic distance (km)",
         y = "Jaccard dissimilarity") +
    geom_smooth(method = "gam", 
                colour = "red", 
                alpha = 0.2, 
                fill = "red") +
    theme( axis.text.x = element_text(face = "bold",
                                      colour = "black",
                                      size = 12), 
           axis.text.y = element_text(face = "bold",
                                      size = 11, 
                                      colour = "black"), 
           axis.title = element_text(face = "bold", 
                                     size = 14, 
                                     colour = "black"), 
           panel.background = element_blank(), 
           panel.border = element_rect(fill = NA,
                                       colour = "black")) 
  
  ##### nest #####
  aa = as.vector(matrix.nest)
  tt = as.vector(matrix.dist)
  #new data frame with vectorized distance matrices
  mat = data.frame(aa,tt)
  
  mm2 = ggplot(mat, aes(y = aa, x = tt)) + 
    geom_point(size = 3, alpha = 0.5, color = "black") + 
    labs(x = "Geographic distance (km)",
         y = "Nestedness component") +
    geom_smooth(method = "gam", 
                colour = "red", 
                alpha = 0.2, 
                fill = "red") +
    theme( axis.text.x = element_text(face = "bold",
                                      colour = "black",
                                      size = 12), 
           axis.text.y = element_text(face = "bold",
                                      size = 11, 
                                      colour = "black"), 
           axis.title = element_text(face = "bold", 
                                     size = 14, 
                                     colour = "black"), 
           panel.background = element_blank(), 
           panel.border = element_rect(fill = NA,
                                       colour = "black")) 

  ##### turn #####
  aa = as.vector(matrix.turn)
  tt = as.vector(matrix.dist)
  #new data frame with vectorized distance matrices
  mat = data.frame(aa,tt)
  
  mm3 = ggplot(mat, aes(y = aa, x = tt)) + 
    geom_point(size = 3, alpha = 0.5, color = "black") + 
    labs(x = "Geographic distance (km)",
         y = "Turnover component") +
    geom_smooth(method = "gam", 
                colour = "red", 
                alpha = 0.2, 
                fill = "red") +
    theme( axis.text.x = element_text(face = "bold",
                                      colour = "black",
                                      size = 12), 
           axis.text.y = element_text(face = "bold",
                                      size = 11, 
                                      colour = "black"), 
           axis.title = element_text(face = "bold", 
                                     size = 14, 
                                     colour = "black"), 
           panel.background = element_blank(), 
           panel.border = element_rect(fill = NA,
                                       colour = "black")) 
  
  ##### bray #####
  aa = as.vector(matrix.bray)
  tt = as.vector(matrix.dist)
  #new data frame with vectorized distance matrices
  mat = data.frame(aa,tt)
  
  mm4 = ggplot(mat, aes(y = aa, x = tt)) + 
    geom_point(size = 3, alpha = 0.5, color = "black") + 
    labs(x = "Geographic distance (km)",
         y = "Bray-Curtis dissimilarity") +
    geom_smooth(method = "gam", 
                colour = "red", 
                alpha = 0.2, 
                fill = "red") +
    theme( axis.text.x = element_text(face = "bold",
                                      colour = "black",
                                      size = 12), 
           axis.text.y = element_text(face = "bold",
                                      size = 11, 
                                      colour = "black"), 
           axis.title = element_text(face = "bold", 
                                     size = 14, 
                                     colour = "black"), 
           panel.background = element_blank(), 
           panel.border = element_rect(fill = NA,
                                       colour = "black")) 
  ##### gra #####
  aa = as.vector(matrix.gra)
  tt = as.vector(matrix.dist)
  #new data frame with vectorized distance matrices
  mat = data.frame(aa,tt)
  
  mm5 = ggplot(mat, aes(y = aa, x = tt)) + 
    geom_point(size = 3, alpha = 0.5, color = "black") + 
    labs(x = "Geographic distance (km)",
         y = "Abundance gradiant component") +
    geom_smooth(method = "gam", 
                colour = "red", 
                alpha = 0.2, 
                fill = "red") +
    theme( axis.text.x = element_text(face = "bold",
                                      colour = "black",
                                      size = 12), 
           axis.text.y = element_text(face = "bold",
                                      size = 11, 
                                      colour = "black"), 
           axis.title = element_text(face = "bold", 
                                     size = 14, 
                                     colour = "black"), 
           panel.background = element_blank(), 
           panel.border = element_rect(fill = NA,
                                       colour = "black")) 
  
  ##### bal #####
  aa = as.vector(matrix.bal)
  tt = as.vector(matrix.dist)
  #new data frame with vectorized distance matrices
  mat = data.frame(aa,tt)
  
  mm6 = ggplot(mat, aes(y = aa, x = tt)) + 
    geom_point(size = 3, alpha = 0.5, color = "black") + 
    labs(x = "Geographic distance (km)",
         y = "Balanced variation component") +
    geom_smooth(method = "gam", 
                colour = "red", 
                alpha = 0.2, 
                fill = "red") +
    theme( axis.text.x = element_text(face = "bold",
                                      colour = "black",
                                      size = 12), 
           axis.text.y = element_text(face = "bold",
                                      size = 11, 
                                      colour = "black"), 
           axis.title = element_text(face = "bold", 
                                     size = 14, 
                                     colour = "black"), 
           panel.background = element_blank(), 
           panel.border = element_rect(fill = NA,
                                       colour = "black")) 
  
  ##### plot #####
  cowplot::plot_grid(mm1, mm2, mm3, mm4, mm5, mm6, ncol = 3, nrow = 2)
  #### heatmap ####
  
  source("R/coldiss.R")
  library(gclus)
  library(vegan)
  coldiss(as.dist(decostand(matrix.turn,"standardize")), byrank = TRUE, diag = TRUE)
  
  }
