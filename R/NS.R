#' Triangle plot and mantel test/graphs
#'
#' @param meta_and_data the path to the metadata and data
#' @param dat_thresh_path the path to the file data red
#'https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12029
#' @return path to the data file, reduced and mean by site
#' @export
#'
#'


north_south <- function (meta_and_data, mean_by_arms, arms_id) {
  
  #mean_by_arms = targets::tar_read("arms_mean") 
  #meta_and_data = targets::tar_read("metadata_data") 
  #arms_id = "RUNA"
  
  
  #### data load and transform ####
  data <- read.csv(mean_by_arms)
  rownames(data) <- data[,1]
  data <- data[,-1]
  matrix.hel <- vegan::decostand(data, "hellinger")
  
  vector <- c(substr(row.names(matrix.hel),1,5))
  
  tab <- aggregate(matrix.hel, by = list(vector), mean)
  
  rownames(tab) <- tab$Group.1
  matrix.h.ag <- tab[,-1]
  
  #### Computing Bdiv comp (abun based) ####
  B.pair.abund <- betapart::beta.pair.abund(matrix.h.ag)
  
  diag.bal <- as.matrix(B.pair.abund$beta.bray.bal)
  diag.bal <- as.data.frame(diag.bal)
  diag.bal <- diag.bal[-1,]
  diag.bal <- as.matrix(diag.bal)
  diag.bal <- diag(diag.bal)
  
  diag.gra <- as.matrix(B.pair.abund$beta.bray.gra)
  diag.gra <- as.data.frame(diag.gra)
  diag.gra <- diag.gra[-1,]
  diag.gra <- as.matrix(diag.gra)
  diag.gra <- diag(diag.gra)
  
  diag.bray <- as.matrix(B.pair.abund$beta.bray)
  diag.bray <- as.data.frame(diag.bray)
  diag.bray <- diag.bray[-1,]
  diag.bray <- as.matrix(diag.bray)
  diag.bray <- diag(diag.bray)

  #### Computing Bdiv comp (pa based) ####
  matrix.pa.ag <- vegan::decostand(matrix.h.ag, "pa")
  B.pair.pa <- betapart::beta.pair(matrix.pa.ag, "jaccard")
  
  diag.turn <- as.matrix(B.pair.pa$beta.jtu)
  diag.turn <- as.data.frame(diag.turn)
  diag.turn <- diag.turn[-1,]
  diag.turn <- as.matrix(diag.turn)
  diag.turn <- diag(diag.turn)
  
  diag.nest <- as.matrix(B.pair.pa$beta.jne)
  diag.nest <- as.data.frame(diag.nest)
  diag.nest <- diag.nest[-1,]
  diag.nest <- as.matrix(diag.nest)
  diag.nest <- diag(diag.nest)
  
  diag.jacc <- as.matrix(B.pair.pa$beta.jac)
  diag.jacc <- as.data.frame(diag.jacc)
  diag.jacc <-diag.jacc[-1,]
  diag.jacc <- as.matrix(diag.jacc)
  diag.jacc <- diag(diag.jacc)
 
  
  north_south_name <- paste0("north_south_", arms_id, ".pdf")
  north_south_path <- here::here("outputs",  north_south_name)
  pdf(file =  north_south_path, width = 11, height = 7)
  par(mfrow = c(1,2))
  #### Plot abundance based ####
  label.pairs = c("RUNA1-2", 
                  "RUNA2-3", 
                  "RUNA3-4", 
                  "RUNA4-5", 
                  "RUNA5-6", 
                  "RUNA6-7", 
                  "RUNA7-8",
                  "RUNA8-9")
  
 
  plot(diag.bray, xaxt = "n", xlab = "", ylim = c(0,0.4) )
  
  points(diag.bray,
          pch = 16,
          cex = 2,
          col = "black",
          bg = "black")
  lines(diag.bray, col = "black")
  
  points(diag.bal,
    pch = 22,
    cex = 2,
    col = "white",
    bg = "blue")
  lines(diag.bal, col = "blue") 

  points(diag.gra,
    pch = 24,
    cex = 2,
    col = "white",
    bg = "red")
  lines(diag.gra, col = "red") 
  
  legend("top", c("Bray-Curtis Disimilarity", "Balanced variation", "Abundance gradiant"), pch = c(16, 15, 17), col = c("black", "blue", "red"), cex = 0.7)
  
  axis(side = 1, 1:8, labels = label.pairs, las = 2, cex.axis = 0.9)
  
  #### Graph pa based ####
  plot(diag.jacc, xaxt = "n", xlab = "", ylim = c(0,0.6) )
  
  points(diag.jacc,
         pch = 16,
         cex = 2,
         col = "black",
         bg = "black")
  lines(diag.jacc, col = "black")
  
  points(diag.nest,
         pch = 22,
         cex = 2,
         col = "white",
         bg = "blue")
  lines(diag.nest, col = "blue") 
  
  points(diag.turn,
         pch = 24,
         cex = 2,
         col = "white",
         bg = "red")
  lines(diag.turn, col = "red") 
  
  legend("top", c("Jaccard Disimilarity", "Nestedness", "Turnover"), pch = c(16, 15, 17), col = c("black", "blue", "red"), cex = 0.7)
  
  axis(side = 1, 1:8, labels = label.pairs, las = 2, cex.axis = 0.9)
  dev.off()
  #### coldiss ####
  heatmap_name <- paste0("heatmap_", arms_id, ".pdf")
  heatmap_path <- here::here("outputs",  heatmap_name)
  pdf(file = heatmap_path, width = 10, height = 10)

  source("R/coldiss.R")
  library(gclus)
  library(vegan)
  library(viridis)
 
  coldiss(as.dist(scale(as.matrix(B.pair.abund$beta.bray))), byrank = TRUE, diag = FALSE, nc = 20)
  mtext("Bray-Curtis index", side = 3)
  #coldiss(as.dist(scale(as.matrix(B.pair.abund$beta.bray.bal))), byrank = TRUE, diag = FALSE, nc = 20)
  #mtext("Balanced variation", side = 3)
  #coldiss(as.dist(scale(as.matrix(B.pair.abund$beta.bray.gra))), byrank = TRUE, diag = FALSE, nc = 20)
  #mtext("Abundance gradiant", side = 3)
  
  
  #coldiss(as.dist(scale(as.matrix(B.pair.pa$beta.jac))), byrank = TRUE, diag = FALSE, nc = 20)
  #mtext("Jaccard index", side = 3)
  #coldiss(as.dist(scale(as.matrix(B.pair.pa$beta.jtu))), byrank = TRUE, diag = TRUE, nc = 20)
  #coldiss(as.dist(scale(as.matrix(B.pair.pa$beta.jne))), byrank = TRUE, diag = FALSE, nc = 20)
  
  
  dev.off()
  return(c(north_south_path,heatmap_path))
}

