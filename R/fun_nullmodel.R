#' beta diversity decomposition
#'
#' @param metadata_data_mean data
#' @return the path to the subseted raw data file
#' @export

fun_null_model <- function(metadata_data_mean){
  
  # metadata_data_mean = targets::tar_read(mean_metadata_data)
  library(ggpubr)
  library(forcats)
  library(reshape2)
  library(ggplot2)
  library(EcoSimR)
  
  
  df_mean <- read.csv(metadata_data_mean[!grepl("metadata", metadata_data_mean)], header = TRUE)
  meta_mean <- read.csv(metadata_data_mean[grepl("metadata", metadata_data_mean)], header = TRUE)
  
  df_mean <- vegan::decostand(df_mean, "pa")
  
  #### Compute beta-diversity index on null model ####
  
  # Define the number of iterations for the null model
  num_iterations <- 1000  # You can adjust this number
  
  # Run the null model using the swapping algorithm
  null_model <- sim9(df_mean, nReps = num_iterations, metric = "c_score", algo = "sim9")
  
    ### parameters justification ###
    # C-Score :
    # EcoSim generates 1000 random matrices as the default. The default 
    # co-occurrence index is Stone and Robert's (1990) C-score. The default 
    # randomization algorithm maintains fixed sums for row and column constraints. 
    # Thus, each matrix generated has the same row and column totals as the
    # original matrix (Connor and Simberloff 1979). As described in the tutorial, 
    # this algorithm has good Type I properties (low chance of falsely rejecting 
    # the null hypothesis when it is true), but also has good power for detecting 
    # non-random patterns in noisy data sets.
    # sim9 : 
    # Fixed rows-fixed columns This simulation maintains fixed row and column 
    # sums. Thus, no degenerate matrices are produced. Although an earlier 
    # version of this model by Connor and Simberloff (1979) was widely 
    # criticized, the version implmented in EcoSim has a good Type I error rate, 
    # and is powerful at detecting patterns in noisy data sets, particularly 
    # when used with the C-score. This model cannot be used with the V-ratio, 
    # which is determined exclusively by row and column sums. RECOMMENDED.
  
    
   
  # Return the null matrix
  null_model_data <- null_model$Randomized.Data
  
  B.pair.pa <- betapart::beta.pair(null_model_data, index.family = "jaccard")
  
  mat.turn <- B.pair.pa$beta.jtu
  mat.nest <- B.pair.pa$beta.jne
  mat.jacc <- B.pair.pa$beta.jac
  
  #turn
  df.turn <- melt(as.matrix(mat.turn), varnames = c("row", "col"))
  df.turn <- subset(df.turn, row != col)
  #nest
  df.nest <- melt(as.matrix(mat.nest), varnames = c("row", "col"))
  df.nest <- subset(df.nest, row != col)
  #jac
  df.jacc <- melt(as.matrix(mat.jacc), varnames = c("row", "col"))
  df.jacc <- subset(df.jacc, row != col)
  
  tab.jacc.null <- as.data.frame(df.jacc)
  tab.nest.null <- as.data.frame(df.nest)
  tab.turn.null <- as.data.frame(df.turn)
  
  #### Compute beta-diversity index on observed data ####
  df_mean <- read.csv(metadata_data_mean[!grepl("metadata", metadata_data_mean)], header = TRUE)
  meta_mean <- read.csv(metadata_data_mean[grepl("metadata", metadata_data_mean)], header = TRUE)

  matrix.pa <- vegan::decostand(df_mean, "pa")
  rownames(matrix.pa) <- meta_mean$arms
  colnames(matrix.pa) <- meta_mean$arms

  obs.B.pair.pa <- betapart::beta.pair(matrix.pa, index.family = "jaccard")
  obs.mat.turn <- obs.B.pair.pa$beta.jtu
  obs.mat.nest <- obs.B.pair.pa$beta.jne
  obs.mat.jacc <- obs.B.pair.pa$beta.jac

  #turn
  obs.df.turn <- melt(as.matrix(obs.mat.turn), varnames = c("row", "col"))
  obs.df.turn <- subset(obs.df.turn, row != col)

  obs.df.turn$row <- substr(obs.df.turn$row, 1, 5)
  obs.df.turn$col <- substr(obs.df.turn$col, 1, 5)
  obs.df.turn$same_value <- ifelse(obs.df.turn$row == obs.df.turn$col, "Yes", "No")

  
  #nest
  obs.df.nest <- melt(as.matrix(obs.mat.nest), varnames = c("row", "col"))
  obs.df.nest <- subset(obs.df.nest, row != col)
  
  obs.df.nest$row <- substr(obs.df.nest$row, 1, 5)
  obs.df.nest$col <- substr(obs.df.nest$col, 1, 5)
  obs.df.nest$same_value <- ifelse(obs.df.nest$row == obs.df.nest$col, "Yes", "No")
  
  #jacc
  obs.df.jacc <- melt(as.matrix(obs.mat.jacc), varnames = c("row", "col"))
  obs.df.jacc <- subset(obs.df.jacc, row != col)
  
  obs.df.jacc$row <- substr(obs.df.jacc$row, 1, 5)
  obs.df.jacc$col <- substr(obs.df.jacc$col, 1, 5)
  obs.df.jacc$same_value <- ifelse(obs.df.jacc$row == obs.df.jacc$col, "Yes", "No")

  #### compare null and observed using t.test ####
  obs.turn <- tapply(obs.df.turn$value, obs.df.turn$same_value, mean)
  obs.nest <- tapply(obs.df.nest$value, obs.df.nest$same_value, mean)
  obs.jacc <- tapply(obs.df.jacc$value, obs.df.jacc$same_value, mean)

  obs.df.jacc.yes <- subset(obs.df.jacc, obs.df.jacc$same_value=="Yes")
  obs.df.jacc.no <- subset(obs.df.jacc, obs.df.jacc$same_value=="No")
  
  obs.df.turn.yes <- subset(obs.df.turn, obs.df.turn$same_value=="Yes")
  obs.df.turn.no <- subset(obs.df.turn, obs.df.turn$same_value=="No")
  
  obs.df.nest.yes <- subset(obs.df.nest, obs.df.nest$same_value=="Yes")
  obs.df.nest.no <- subset(obs.df.nest, obs.df.nest$same_value=="No")
  
  # Condition d'application t.test/wilcox.test
  ggplot(tab.jacc.null, aes(x=value)) + 
    geom_density()
  ggpubr::ggqqplot(tab.jacc.null$value) 
  shapiro.test(tab.jacc.null$value) 
  
  ggplot(tab.nest.null, aes(x=value)) + 
    geom_density()
  ggpubr::ggqqplot(tab.nest.null$value) 
  shapiro.test(tab.nest.null$value) 
  
  ggplot(tab.turn.null, aes(x=value)) + 
    geom_density()
  ggpubr::ggqqplot(tab.turn.null$value) 
  shapiro.test(tab.turn.null$value) 
  
  #wilcox.test
  # Comparison between a mean value and a set of value (conformation test) 

  turn.test.yes <- wilcox.test(x = tab.turn.null$value, mu = obs.turn[2], alternative = "two.sided")
  nest.test.yes <- wilcox.test(x = tab.nest.null$value, mu = obs.nest[2], alternative = "two.sided")
  jacc.test.yes <- wilcox.test(x = tab.jacc.null$value, mu = obs.jacc[2], alternative = "two.sided")
  
  turn.test.no <- wilcox.test(x = tab.turn.null$value, mu = obs.turn[1], alternative = "two.sided")
  nest.test.no <- wilcox.test(x = tab.nest.null$value, mu = obs.nest[1], alternative = "two.sided")
  jacc.test.no <- wilcox.test(x = tab.jacc.null$value, mu = obs.jacc[1], alternative = "two.sided")
  # Comparison between two set of value (mean comparison)
  
  install.packages("effsize")
  library(effsize)
  
  cohen_d <- cohen.d(tab.turn.null$value, obs.df.turn.no$value)
  
  turn.test.yes <- wilcox.test(x = tab.turn.null$value, mu = obs.turn[2], alternative = "two.sided")
  nest.test.yes <- wilcox.test(x = tab.nest.null$value, mu = obs.nest[2], alternative = "two.sided")
  jacc.test.yes <- wilcox.test(x = tab.jacc.null$value, mu = obs.jacc[2], alternative = "two.sided")
  
  
  
  
  #### boucle ####
  # tab.jacc <- matrix(nrow = n, ncol = 1)
  # 
  # tab.nest <- matrix(nrow = n, ncol = 1)
  # 
  # tab.turn <- matrix(nrow = n, ncol = 1)
  # for (i in 1:999) {
  #   
  #   df_mean <- read.csv(metadata_data_mean[!grepl("metadata", metadata_data_mean)], header = TRUE)
  #   df_mean <- vegan::decostand(df_mean, "pa")
  #   
  #   # Run the null model using the swapping algorithm
  #   null_model <- sim9(df_mean, nReps = 1000, metric = "c_score", algo = "sim9")
  #   
  #   # Return the null matrix
  #   null_model_data <- null_model$Randomized.Data
  #   
  #   B.pair.pa <- betapart::beta.pair(null_model_data, index.family = "jaccard")
  #   
  #   mat.turn <- B.pair.pa$beta.jtu
  #   mat.nest <- B.pair.pa$beta.jne
  #   mat.jacc <- B.pair.pa$beta.jac
  #   
  #   ####  inter/intra par set ####
  #   #turn
  #   df.turn <- melt(as.matrix(mat.turn), varnames = c("row", "col"))
  #   df.turn <- subset(df.turn, row != col)
  #   
  #   # Créez un tableau vide de longueur n
  #   tab.turn[i,] <- mean(df.turn$value)
  #   
  #   #nest
  #   df.nest <- melt(as.matrix(mat.nest), varnames = c("row", "col"))
  #   df.nest <- subset(df.nest, row != col)
  #   
  #   tab.nest[i,] <- mean(df.nest$value)
  #   
  #   #jac
  #   df.jacc <- melt(as.matrix(mat.jacc), varnames = c("row", "col"))
  #   df.jacc <- subset(df.jacc, row != col)
  #   
  #   tab.jacc[i,] <- mean(df.jacc$value)
  #   
  #   tab.jacc.null <- as.data.frame(tab.jacc)
  #   tab.nest.null <- as.data.frame(tab.nest)
  #   tab.turn.null <- as.data.frame(tab.turn)
  #   
  # }
  # 
  # mean.jacc.null <- mean(tab.jacc.null[,1])
  # mean.nest.null <- mean(tab.nest.null[,1])
  # mean.turn.null <- mean(tab.turn.null[,1])
  


  return(meta_mean)
  return (cohen_d)
}