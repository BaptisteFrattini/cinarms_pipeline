#' beta diversity decomposition
#'
#' @param metadata_data_mean data
#' @return the path to the subseted raw data file
#' @export

beta_div_decomp <- function(metadata_data_mean){
  
  # metadata_data_mean = targets::tar_read(mean_metadata_data)
  library(ggpubr)
  library(forcats)
  library(reshape2)
  library(ggplot2)
  df_mean <- read.csv(metadata_data_mean[!grepl("metadata", metadata_data_mean)], header = TRUE)
  meta_mean <- read.csv(metadata_data_mean[grepl("metadata", metadata_data_mean)], header = TRUE)
  
  #   # Return the null matrix
  # 
  # 
  # 
  # n=999
  # tab.jacc <- matrix(nrow = n, ncol = 1)
  # 
  # tab.nest <- matrix(nrow = n, ncol = 1) 
  # 
  # tab.turn <- matrix(nrow = n, ncol = 1)  
  # 
  # #### boucle ####
  # 
  # for (i in 1:999) {
  # 
  # df_mean <- read.csv(metadata_data_mean[!grepl("metadata", metadata_data_mean)], header = TRUE)
  # df_mean <- vegan::decostand(df_mean, "pa")
  # null_models <- vegan::nullmodel(df_mean, method = "swap")
  # sim <- simulate(null_models, nsim = 1, seed = NULL,
  #          burnin = 0, thin = 1)
  # 
  # #resampling
  # df_mean <- permute_columns_preserve_richness(df_mean)
  # matrix.pa <- df_mean
  # rownames(matrix.pa) <- meta_mean$arms
  # 
  # B.pair.pa <- betapart::beta.pair(matrix.pa, index.family = "jaccard")
  # mat.turn <- B.pair.pa$beta.jtu
  # mat.nest <- B.pair.pa$beta.jne
  # mat.jacc <- B.pair.pa$beta.jac
  # 
  # ####  inter/intra par set ####
  # #turn
  # df.turn <- melt(as.matrix(mat.turn), varnames = c("row", "col"))
  # df.turn <- subset(df.turn, row != col)
  # 
  # # Créez un tableau vide de longueur n
  # tab.turn[i,] <- mean(df.turn$value)
  # 
  # #nest
  # df.nest <- melt(as.matrix(mat.nest), varnames = c("row", "col"))
  # df.nest <- subset(df.nest, row != col)
  # 
  # tab.nest[i,] <- mean(df.nest$value)
  # 
  # #jac
  # df.jacc <- melt(as.matrix(mat.jacc), varnames = c("row", "col"))
  # df.jacc <- subset(df.jacc, row != col)
  # 
  # tab.jacc[i,] <- mean(df.jacc$value)
  # 
  # tab.jacc.null <- as.data.frame(tab.jacc)
  # tab.nest.null <- as.data.frame(tab.nest)
  # tab.turn.null <- as.data.frame(tab.turn)
  # 
  # }
  # 
  # mean.jacc.null <- mean(tab.jacc.null[,1])
  # mean.nest.null <- mean(tab.nest.null[,1])
  # mean.turn.null <- mean(tab.turn.null[,1])
  # 
  # #### compute obsrved values ####
  # df_mean <- read.csv(metadata_data_mean[!grepl("metadata", metadata_data_mean)], header = TRUE)
  # meta_mean <- read.csv(metadata_data_mean[grepl("metadata", metadata_data_mean)], header = TRUE)
  # 
  # matrix.pa <- vegan::decostand(df_mean, "pa")
  # rownames(matrix.pa) <- meta_mean$arms
  # colnames(matrix.pa) <- meta_mean$arms
  # 
  # B.pair.pa <- betapart::beta.pair(matrix.pa, index.family = "jaccard")
  # mat.turn <- B.pair.pa$beta.jtu
  # mat.nest <- B.pair.pa$beta.jne
  # mat.jacc <- B.pair.pa$beta.jac
  # 
  # #turn
  # df.turn <- melt(as.matrix(mat.turn), varnames = c("row", "col"))
  # df.turn <- subset(df.turn, row != col)
  # 
  # df.turn$row <- substr(df.turn$row, 1, 5)
  # df.turn$col <- substr(df.turn$col, 1, 5)
  # df.turn$same_value <- ifelse(df.turn$row == df.turn$col, "Yes", "No")
  # 
  # #nest
  # df.nest <- melt(as.matrix(mat.nest), varnames = c("row", "col"))
  # df.nest <- subset(df.nest, row != col)
  # 
  # df.nest$row <- substr(df.nest$row, 1, 5)
  # df.nest$col <- substr(df.nest$col, 1, 5)
  # df.nest$same_value <- ifelse(df.nest$row == df.nest$col, "Yes", "No")
  # 
  # #jac
  # df.jacc <- melt(as.matrix(mat.jacc), varnames = c("row", "col"))
  # df.jacc <- subset(df.jacc, row != col)
  # 
  # df.jacc$row <- substr(df.jacc$row, 1, 5)
  # df.jacc$col <- substr(df.jacc$col, 1, 5)
  # df.jacc$same_value <- ifelse(df.jacc$row == df.jacc$col, "Yes", "No")
  # 
  # #### compare null and observed using t.test ####
  # turn <- tapply(df.turn$value, df.turn$same_value, mean)
  # nest <- tapply(df.nest$value, df.nest$same_value, mean)
  # jacc <- tapply(df.jacc$value, df.jacc$same_value, mean)
  # 
  # 
  # turn.test.yes <- t.test(x = tab.turn.null$Yes, mu = turn[2], alternative = "two.sided")
  # nest.test.yes <- t.test(x = tab.nest.null$Yes, mu = nest[2], alternative = "two.sided")
  # jacc.test.yes <- t.test(x = tab.jacc.null$Yes, mu = jacc[2], alternative = "two.sided")
  # 
  # turn.test.no <- t.test(x = tab.turn.null$No, mu = turn[1], alternative = "two.sided")
  # nest.test.no <- t.test(x = tab.nest.null$No, mu = nest[1], alternative = "two.sided")
  # jacc.test.no <- t.test(x = tab.jacc.null$No, mu = jacc[1], alternative = "two.sided")
  # 
  # 
  # 
  return(meta_mean)
}