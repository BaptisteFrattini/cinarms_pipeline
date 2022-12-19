#' Subset a raw_data table for a sampling campain
#'
#' @param meta_and_data the path to the metadata and data
#' @param dat_thresh_path the path to the file data red
#'https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12029
#' @return path to the data file, reduced and mean by site
#' @export
#'

decomp_b_div <- function(mean_by_arms) {
  #mean_by_arms = targets::tar_read("arms_mean") 
  data <- read.csv(mean_by_arms)
  rownames(data) <- data[,1]
  data <- data[,-1]
  matrix.hel <- vegan::decostand(data, "hellinger")
  #Jaccard
  decomp.bas.J <- adespatial::beta.div.comp(matrix.hel, coef = "BJ", quant = FALSE) 
 
  decomp.bas.J.3 <- cbind((1 - decomp.bas.J$D), decomp.bas.J$repl, decomp.bas.J$rich)

  colnames(decomp.bas.J.3) <- c("Similarity", "Repl", "RichDiff") 
  
  ade4::triangle.plot(as.data.frame(decomp.bas.J.3[, c(3, 1, 2)]), 
                                 labeltriangle = FALSE, 
                                 addmean = TRUE,
                                 show.position = FALSE)
  
  text(-0.45, 0.5, "Nestedness", cex = 1.5)
  text(0.4, 0.5, "Turnover", cex = 1.5)
  text(0, -0.6, "Jaccard similarity", cex = 1.5) 
  
  stats::heatmap(as.matrix(decomp.bas.J$rich))
  
  betaset <- betapart::beta.pair(matrix.pa, index.family = "jaccard")

  summary(bd)
  plot(bd, label.cex = 0.5)
  anova(bd)
  permutest(bd)
  bd
  
  
  
  
  
  
  
}
