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
  library(effsize)
  library(dplyr)
  # install.packages("effsize")
  
  df_mean <- read.csv(metadata_data_mean[!grepl("metadata", metadata_data_mean)], header = TRUE)
  meta_mean <- read.csv(metadata_data_mean[grepl("metadata", metadata_data_mean)], header = TRUE)
  
  df_mean <- vegan::decostand(df_mean, "pa")
  
  #### ab based
  

  # out <- vegan::nullmodel(df_mean, "quasiswap")
  # 
  # b <- stats::simulate(out, nsim = 1000, seed = NULL,
  #          burnin = 0, thin = 1)
  # summary(b)
  
  # df_mean <- round(df_mean,0)
  
  # df_mean <- vegan::decostand(df_mean, "total")
  
  # rowSums(df_mean)
  # ?vegan::decostand
  # x <- vegan::permatswap(df_mean, method = "swsh_samp", fixedmar="both", shuffle = "both",
  #            strata = NULL, mtype = "count", times = 99, 
  #            burnin = 0, thin = 1)
  # summary(x)
  # x$perm
  #### Compute beta-diversity index on null model ####
  
  # Define the number of iterations for the null model
  num_iterations <- 1000  # You can adjust this number
  
  # Run the null model using the swapping algorithm
  df_mean <- t(df_mean)
  # null_model <- sim9(df_mean, nReps = num_iterations, metric = "c_score", algo = "sim9")
  null_model <- sim9(df_mean, nReps = num_iterations, metric = "c_score", algo = "sim2")
  ?sim9
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
  
  plot(null_model,type="burn_in")
  plot(null_model,type="hist")
  plot(null_model,type="cooc")
  
  summary(null_model)
  # 
  # plot(null_model_2,type="burn_in")
  # plot(null_model_2,type="hist")
  # plot(null_model_2,type="cooc")
  #   
  # summary(null_model_2)
  
  # Global SES using C Score = 6 (>>1,96) --> In community assembly, it would significate 
  # that niche processes are playing a role. Since we have different parameters 
  # of immersion time and season, "over dispersion" is attributed to these parameters.
  
  # Return the null matrix
  null_model_data <- t(null_model$Randomized.Data)
  rowSums(null_model_data)
  colSums(null_model_data)
  rownames(null_model_data) <- meta_mean$arms
  
  B.pair.pa <- betapart::beta.pair(null_model_data, index.family = "jaccard")
  
  mat.turn <- B.pair.pa$beta.jtu
  mat.nest <- B.pair.pa$beta.jne
  mat.jacc <- B.pair.pa$beta.jac
  
  #turn
  df.turn <- melt(as.matrix(mat.turn), varnames = c("row", "col"))
  df.turn <- subset(df.turn, row != col)
  df.turn <- df.turn %>%
    mutate(
      row = as.character(row),  # Convert 'row' to character
      col = as.character(col),  # Convert 'col' to character
      combined = paste(pmin(row, col), pmax(row, col), sep = "")
    )
  df.turn <- df.turn[!duplicated(df.turn$combined), ]
  nrow(df.turn)
  
  #nest
  df.nest <- melt(as.matrix(mat.nest), varnames = c("row", "col"))
  df.nest <- subset(df.nest, row != col)
  df.nest <- df.nest %>%
    mutate(
      row = as.character(row),  # Convert 'row' to character
      col = as.character(col),  # Convert 'col' to character
      combined = paste(pmin(row, col), pmax(row, col), sep = "")
    )
  df.nest <- df.nest[!duplicated(df.nest$combined), ]
  nrow(df.nest)
  
  #jac
  df.jacc <- melt(as.matrix(mat.jacc), varnames = c("row", "col"))
  df.jacc <- subset(df.jacc, row != col)
  df.jacc <- df.jacc %>%
    mutate(
      row = as.character(row),  # Convert 'row' to character
      col = as.character(col),  # Convert 'col' to character
      combined = paste(pmin(row, col), pmax(row, col), sep = "")
    )
  df.jacc <- df.jacc[!duplicated(df.jacc$combined), ]
  nrow(df.jacc)
  
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
  obs.df.turn <- obs.df.turn %>%
    mutate(
      row = as.character(row),  # Convert 'row' to character
      col = as.character(col),  # Convert 'col' to character
      combined = paste(pmin(row, col), pmax(row, col), sep = "")
    )
  obs.df.turn <- obs.df.turn[!duplicated(obs.df.turn$combined), ]
  nrow(obs.df.turn)
  
  
  obs.df.turn$row <- substr(obs.df.turn$row, 1, 5)
  obs.df.turn$col <- substr(obs.df.turn$col, 1, 5)
  obs.df.turn$same_value <- ifelse(obs.df.turn$row == obs.df.turn$col, "Yes", "No")

  
  #nest
  obs.df.nest <- melt(as.matrix(obs.mat.nest), varnames = c("row", "col"))
  obs.df.nest <- subset(obs.df.nest, row != col)
  obs.df.nest <- obs.df.nest %>%
    mutate(
      row = as.character(row),  # Convert 'row' to character
      col = as.character(col),  # Convert 'col' to character
      combined = paste(pmin(row, col), pmax(row, col), sep = "")
    )
  obs.df.nest <- obs.df.nest[!duplicated(obs.df.nest$combined), ]
  nrow(obs.df.nest)
  
  obs.df.nest$row <- substr(obs.df.nest$row, 1, 5)
  obs.df.nest$col <- substr(obs.df.nest$col, 1, 5)
  obs.df.nest$same_value <- ifelse(obs.df.nest$row == obs.df.nest$col, "Yes", "No")
  
  #jacc
  obs.df.jacc <- melt(as.matrix(obs.mat.jacc), varnames = c("row", "col"))
  obs.df.jacc <- subset(obs.df.jacc, row != col)
  obs.df.jacc <- obs.df.jacc %>%
    mutate(
      row = as.character(row),  # Convert 'row' to character
      col = as.character(col),  # Convert 'col' to character
      combined = paste(pmin(row, col), pmax(row, col), sep = "")
    )
  obs.df.jacc <- obs.df.jacc[!duplicated(obs.df.jacc$combined), ]
  nrow(obs.df.jacc)
  
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
  
  #### Compute the null deviation data ####
  
  #### JACCARD
  
  nrow(obs.df.jacc)
  nrow(tab.jacc.null)
  null.dev.jacc <- (obs.df.jacc$value - mean(tab.jacc.null$value))/sd(tab.jacc.null$value)
  
  tab.null.dev.jacc <- data.frame(obs.df.jacc$row, obs.df.jacc$col, obs.df.jacc$same_value, null.dev.jacc)
  
  tapply(tab.null.dev.jacc$null.dev, tab.null.dev.jacc$obs.df.jacc.same_value, mean)
  tapply(tab.null.dev.jacc$null.dev, tab.null.dev.jacc$obs.df.jacc.same_value, sd)
  
  # As proposed by Oscar B. Vitorino Júnior et al, 2016 or Tom R. Bishop et al, (2015)
  
  ggplot(tab.null.dev.jacc, aes(x=null.dev.jacc)) +
  geom_density()
  
  #SES values significantly different from random expectations
  subset_tab.null.dev.jacc <- tab.null.dev.jacc[tab.null.dev.jacc$null.dev > 1.96 | tab.null.dev.jacc$null.dev < -1.96, ] 
  
  tapply(subset_tab.null.dev.jacc$null.dev, subset_tab.null.dev.jacc$obs.df.jacc.same_value, mean)
  tapply(subset_tab.null.dev.jacc$null.dev, subset_tab.null.dev.jacc$obs.df.jacc.same_value, sd)
  
  #### TURNOVER
  
  nrow(obs.df.turn)
  nrow(tab.turn.null)
  
  null.dev.turn <- (obs.df.turn$value - mean(tab.turn.null$value))/sd(tab.turn.null$value)
  
  tab.null.dev.turn <- data.frame(obs.df.turn$row, obs.df.turn$col, obs.df.turn$same_value, null.dev.turn)
  
  tapply(tab.null.dev.turn$null.dev, tab.null.dev.turn$obs.df.turn.same_value, mean)
  tapply(tab.null.dev.turn$null.dev, tab.null.dev.turn$obs.df.turn.same_value, sd)
  
  # As proposed by Oscar B. Vitorino Júnior et al, 2016 or Tom R. Bishop et al, (2015)
  
  ggplot(tab.null.dev.turn, aes(x=null.dev.turn)) + 
    geom_density()
  
  #SES values significantly different from random expectations
  subset_tab.null.dev.turn <- tab.null.dev.turn[tab.null.dev.turn$null.dev > 1.96 | tab.null.dev.turn$null.dev < -1.96, ] 
  
  tapply(subset_tab.null.dev.turn$null.dev, subset_tab.null.dev.turn$obs.df.turn.same_value, mean)
  tapply(subset_tab.null.dev.turn$null.dev, subset_tab.null.dev.turn$obs.df.turn.same_value, sd)
  
  #### NESTEDNESS
  
  nrow(obs.df.nest)
  nrow(tab.nest.null)
  
  null.dev.nest <- (obs.df.nest$value - mean(tab.nest.null$value))/sd(tab.nest.null$value)
  
  tab.null.dev.nest <- data.frame(obs.df.nest$row, obs.df.nest$col, obs.df.nest$same_value, null.dev.nest)
  
  tapply(tab.null.dev.nest$null.dev, tab.null.dev.nest$obs.df.nest.same_value, mean)
  tapply(tab.null.dev.nest$null.dev, tab.null.dev.nest$obs.df.nest.same_value, sd)
  
  # As proposed by Oscar B. Vitorino Júnior et al, 2016 or Tom R. Bishop et al, (2015)
  
  ggplot(tab.null.dev.nest, aes(x=null.dev.nest)) + 
    geom_density()
  
  #SES values significantly different from random expectations
  subset_tab.null.dev.nest <- tab.null.dev.nest[tab.null.dev.nest$null.dev > 1.96 | tab.null.dev.nest$null.dev < -1.96, ] 
  
  tapply(subset_tab.null.dev.nest$null.dev, subset_tab.null.dev.nest$obs.df.nest.same_value, mean)
  tapply(subset_tab.null.dev.nest$null.dev, subset_tab.null.dev.nest$obs.df.nest.same_value, sd)
  
  
  #### Represent results ####
  
  #### intra set
  subset_tab.null.dev.jacc
  subset_tab.null.dev.turn
  subset_tab.null.dev.nest
  
  #### JACCARD
  
  subset.tab.null.dev.jacc <- subset(tab.null.dev.jacc, tab.null.dev.jacc$obs.df.jacc.same_value == "Yes")
  intra = c("between ARMS of \n the CINA1 set", "between ARMS of \n the CINA3 set","between ARMS of \n the CINA2 set","between ARMS of \n the CINA4 set","between ARMS of \n the RUNA2 set")
  kk <- ggplot(subset.tab.null.dev.jacc, aes(x = fct_relevel(obs.df.jacc.col, "CINA1", "CINA3", "CINA2", "CINA4", "RUNA2"), y = null.dev.jacc)) +
    geom_boxplot(fill =  c("#CC66CC","#CC66CC","#1B9E77","#1B9E77","#FF7F00") ) +
    labs(title = "Jaccard dissimilarity",
         x = "Comparisons",
         y = "Standardized Effect Size (SES)") +
    theme(legend.position = "none") +
    scale_x_discrete(labels=intra) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.title.x = element_blank(), axis.title.y = element_text(size=12)) +
    geom_hline(yintercept = -1.96, colour = "red")+
    geom_hline(yintercept = 1.96, colour = "red") +
    geom_hline(yintercept = 0, colour = "darkgrey") +
    ylim(min = -4, max = 4)
  
  #### TURNOVER
  
  subset.tab.null.dev.turn <- subset(tab.null.dev.turn, tab.null.dev.turn$obs.df.turn.same_value == "Yes")
  intra = c("between ARMS of \n the CINA1 set", "between ARMS of \n the CINA3 set","between ARMS of \n the CINA2 set","between ARMS of \n the CINA4 set","between ARMS of \n the RUNA2 set")
  jj <- ggplot(subset.tab.null.dev.turn, aes(x = fct_relevel(obs.df.turn.col, "CINA1", "CINA3", "CINA2", "CINA4", "RUNA2"), y = null.dev.turn)) +
    geom_boxplot(fill =  c("#CC66CC","#CC66CC","#1B9E77","#1B9E77","#FF7F00") ) +
    labs(title = "Turnover component",
         x = "Comparisons",
         y = "Standardized Effect Size (SES)") +
    theme(legend.position = "none") +
    scale_x_discrete(labels=intra) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.title.x = element_blank(), axis.title.y = element_text(size=12)) +
    geom_hline(yintercept = -1.96, colour = "red") +
    geom_hline(yintercept = 1.96, colour = "red") +
    geom_hline(yintercept = 0, colour = "darkgrey")+
    ylim(min = -4, max = 4)
  
  #### NESTEDNESS
  
  subset.tab.null.dev.nest <- subset(tab.null.dev.nest, tab.null.dev.nest$obs.df.nest.same_value == "Yes")
  intra = c("between ARMS of \n the CINA1 set", "between ARMS of \n the CINA3 set","between ARMS of \n the CINA2 set","between ARMS of \n the CINA4 set","between ARMS of \n the RUNA2 set")
  hh <- ggplot(subset.tab.null.dev.nest, aes(x = fct_relevel(obs.df.nest.col, "CINA1", "CINA3", "CINA2", "CINA4", "RUNA2"), y = null.dev.nest)) +
    geom_boxplot(fill =  c("#CC66CC","#CC66CC","#1B9E77","#1B9E77","#FF7F00") ) +
    labs(title = "Nestedness component",
         x = "Comparisons",
         y = "Standardized Effect Size (SES)") +
    theme(legend.position = "none") +
    scale_x_discrete(labels=intra) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.title.x = element_blank(), axis.title.y = element_text(size=12)) +
    geom_hline(yintercept = -1.96, colour = "red") +
    geom_hline(yintercept = 1.96, colour = "red") +
    geom_hline(yintercept = 0, colour = "darkgrey")+
    ylim(min = -4, max = 4)
  
  #### plot
  
  fin <- cowplot::plot_grid(kk,jj,hh,
                            ncol = 1,
                            nrow = 3)
  
  path_to_boxplot <- paste0("outputs/null_model/boxplot_null_model.pdf")
  ggsave(filename =  path_to_boxplot, plot = fin, width = 8, height = 14)
  
  #### intra-inter
  
  #### JACCARD
  intrainter = c("between ARMS of \n the same set", "between ARMS of \n different sets")
  a <- ggplot(tab.null.dev.jacc, aes(x = fct_relevel(obs.df.jacc.same_value, "Yes", "No"), y = null.dev.jacc)) +
    geom_boxplot(fill =  c("lightblue","lightblue") ) +
    labs(title = "Jaccard dissimilarity",
         x = "Comparisons",
         y = "Standardized Effect Size (SES)") +
    theme(legend.position = "none") +
    scale_x_discrete(labels=intrainter) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.title.x = element_blank(), axis.title.y = element_text(size=12)) +
    geom_hline(yintercept = -1.96, colour = "red") +
    geom_hline(yintercept = 1.96, colour = "red") +
    geom_hline(yintercept = 0, colour = "darkgrey")+
    ylim(min = -4, max = 4)
  a
  #### TURNOVER
  intrainter = c("between ARMS of \n the same set", "between ARMS of \n different sets")
  b <- ggplot(tab.null.dev.turn, aes(x = fct_relevel(obs.df.turn.same_value, "Yes", "No"), y = null.dev.turn)) +
    geom_boxplot(fill =  c("lightblue","lightblue") ) +
    labs(title = "Turnover component",
         x = "Comparisons",
         y = "Standardized Effect Size (SES)") +
    theme(legend.position = "none") +
    scale_x_discrete(labels=intrainter) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.title.x = element_blank(), axis.title.y = element_text(size=12)) +
    geom_hline(yintercept = -1.96, colour = "red") +
    geom_hline(yintercept = 1.96, colour = "red") +
    geom_hline(yintercept = 0, colour = "darkgrey")+
    ylim(min = -4, max = 4)
  b
  #### NESTEDNESS
  intrainter = c("between ARMS of \n the same set", "between ARMS of \n different sets")
  c <- ggplot(tab.null.dev.nest, aes(x = fct_relevel(obs.df.nest.same_value, "Yes", "No"), y = null.dev.nest)) +
    geom_boxplot(fill =  c("lightblue","lightblue") ) +
    labs(title = "Nestedness component",
         x = "Comparisons",
         y = "Standardized Effect Size (SES)") +
    theme(legend.position = "none") +
    scale_x_discrete(labels=intrainter) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.title.x = element_blank(), axis.title.y = element_text(size=12)) +
    geom_hline(yintercept = -1.96, colour = "red") +
    geom_hline(yintercept = 1.96, colour = "red") +
    geom_hline(yintercept = 0, colour = "darkgrey")+
    ylim(min = -4, max = 4)
  c
  
  fin <- cowplot::plot_grid(kk,jj,hh,
                            ncol = 3,
                            nrow = 1)
  
  path_to_boxplot <- paste0("outputs/null_model/boxplot_null_model_full.pdf")
  ggsave(filename =  path_to_boxplot, plot = fin, width = 16, height = 7)
  #### bilateral test ####
  
  turn.test.yes <- wilcox.test(x = tab.turn.null$value, mu = obs.turn[2], alternative = "two.sided")
  nest.test.yes <- wilcox.test(x = tab.nest.null$value, mu = obs.nest[2], alternative = "two.sided")
  jacc.test.yes <- wilcox.test(x = tab.jacc.null$value, mu = obs.jacc[2], alternative = "two.sided")
  
  
  return (path_to_boxplot)
}