#' PCOA on turnover component
#'
#' @param metadata_data_mean data
#' @return the path to the subseted raw data file
#' @export

fun_PCA <- function(metadata_data_mean, data_mean_pool){
  
  # metadata_data_mean = targets::tar_read(mean_metadata_data)
  # data_mean_pool = targets::tar_read(data_pool_mean)
  
  #### Load data and meta data ========
  library(ggpubr)
  library(forcats)
  library(reshape2)
  library(ggplot2)
  df_mean <- read.csv(metadata_data_mean[!grepl("metadata", metadata_data_mean)], header = TRUE)
  meta_mean <- read.csv(metadata_data_mean[grepl("metadata", metadata_data_mean)], header = TRUE)
  rownames(df_mean) = meta_mean$arms
  colnames(meta_mean)[5] <- "imm_season"
  colnames(meta_mean)[6] <- "ret_season" 
  
  matrix.pa <- vegan::decostand(df_mean, "pa")
  rownames(matrix.pa) <- meta_mean$arms
  colnames(matrix.pa) <- meta_mean$arms
  
  B.pair.pa <- betapart::beta.pair(matrix.pa, index.family = "jaccard")
  mat.turn <- B.pair.pa$beta.jtu
  mat.nest <- B.pair.pa$beta.jne
  mat.jacc <- B.pair.pa$beta.jac
  mat.bray <- vegan::vegdist(df_mean, method = "bray")
  
  #### With all species ----------
  #### PCOA ####
  #Imm season an imm time
  library(plotly)
  library(ggfortify)
  
  #turn
  pca.res.turn <- prcomp(mat.turn,
                    center = TRUE,
                    scale. = TRUE)
  
  plot_pca_turn <- ggplot2::autoplot(pca.res.turn,
                                data = meta_mean,
                                colour = "imm_time",
                                shape = "imm_season") + 
                   labs(title = "Turnover component")

          
  plot_pca_turn

  #jacc
  pca.res.jacc <- prcomp(mat.jacc,
                         center = TRUE,
                         scale. = TRUE)
  
  plot_pca_jacc <- ggplot2::autoplot(pca.res.jacc,
                                     data = meta_mean,
                                     colour = "imm_time",
                                     shape = "imm_season") + 
                   labs(title = "Jaccard distance")
  plot_pca_jacc
  
  #bray
  pca.res.bray <- prcomp(mat.bray,
                         center = TRUE,
                         scale. = TRUE)
  
  plot_pca_bray <- ggplot2::autoplot(pca.res.bray,
                                     data = meta_mean,
                                     colour = "imm_time",
                                     shape = "imm_season") + 
                   labs(title = "Bray-Curtis distance")
  plot_pca_bray
  
  # retrieval and imm time
  
  pca.res.turn <- prcomp(mat.turn,
                          center = TRUE,
                          scale. = TRUE)
  
  plot_pca_turn.2 <- ggplot2::autoplot(pca.res.turn,
                                     data = meta_mean,
                                     colour = "imm_time",
                                     shape = "ret_season") + 
    labs(title = "Turnover component")
  
  plot_pca_turn.2
  
  pca.res.jacc <- prcomp(mat.jacc,
                         center = TRUE,
                         scale. = TRUE)
  
  plot_pca_jacc.2 <- ggplot2::autoplot(pca.res.jacc,
                                     data = meta_mean,
                                     colour = "imm_time",
                                     shape = "ret_season") + 
    labs(title = "Jaccard distance")
  plot_pca_jacc.2
  
  pca.res.bray <- prcomp(mat.bray,
                         center = TRUE,
                         scale. = TRUE)
  
  plot_pca_bray.2 <- ggplot2::autoplot(pca.res.bray,
                                     data = meta_mean,
                                     colour = "imm_time",
                                     shape = "ret_season") + 
    labs(title = "Bray-Curtis distance")
  plot_pca_bray.2
  
  fin <- cowplot::plot_grid(plot_pca_turn,
                            plot_pca_jacc,
                            plot_pca_bray,
                            plot_pca_turn.2,
                            plot_pca_jacc.2,
                            plot_pca_bray.2,
                            ncol = 3,
                            nrow = 2)
  
  path_to_PCA <- paste0("outputs/beta/PCA.pdf")
  ggsave(filename =  path_to_PCA, plot = fin, width = 19, height = 14)
  
  #with fviz and selection of var 
  sim_imm_tim <- summary(vegan::simper(df_mean, meta_mean$imm_time))
 
  sim_imm_tim_1 <- sim_imm_tim[[1]]
  contrib_imm_tim_1 <- sim_imm_tim_1[(sim_imm_tim_1$p < 0.05) & (sim_imm_tim_1$average > 0.0005),]
  
  sim_imm_tim_2 <- sim_imm_tim[[2]]
  contrib_imm_tim_2 <- sim_imm_tim_2[(sim_imm_tim_2$p < 0.05) & (sim_imm_tim_2$average > 0.0005),]
  
  sim_imm_tim_3 <- sim_imm_tim[[3]]
  contrib_imm_tim_3 <- sim_imm_tim_3[(sim_imm_tim_3$p < 0.05) & (sim_imm_tim_3$average > 0.0005),]
  
  sel <- levels(as.factor(c(rownames(contrib_imm_tim_1), rownames(contrib_imm_tim_2), rownames(contrib_imm_tim_3))))
  
  library(factoextra)
  pca.res.full <- prcomp(df_mean,
                         center = TRUE,
                         scale. = TRUE)
  
  w <- fviz_pca_biplot(pca.res.full,
                       col.ind = meta_mean$imm_time,
                       addEllipses = TRUE,
                       ellipse.type = "convex",
                       select.var = list(contrib = 20),
                       col.var = "purple",
                       repel = TRUE,
                       pointsize = 3,
                       labelsize = 5)
  
  meta_mean$set <- paste0(meta_mean$imm_time,"_",meta_mean$imm_season)
  
  y <- fviz_pca_biplot(pca.res.full,
                       col.ind = meta_mean$set,
                       addEllipses = TRUE,
                       ellipse.type = "convex",
                       select.var = list(contrib = 20),
                       col.var = "purple",
                       repel = TRUE,
                       pointsize = 2,
                       labelsize = 5) + scale_color_manual(values=c("coral","coral","dodgerblue","forestgreen","forestgreen"))
  
  path_to_PCA_select_set <- paste0("outputs/PCA_select_set.pdf")
  ggsave(filename =  path_to_PCA_select_set, plot = y, width = 12, height = 10)
  
  path_to_PCA_select <- paste0("outputs/PCA_select.pdf")
  ggsave(filename =  path_to_PCA_select, plot = w, width = 12, height = 10)
  

  #### With species pool ----------
  #### Load data and meta data ========
  data_pool <- read.csv(data_mean_pool, header = TRUE, row.names = "X")
  mat.bray <- vegan::vegdist(data_pool, method = "bray")
  
  pca.res.bray <- prcomp(data_pool,
                         center = TRUE,
                         scale. = TRUE)
  ?prcomp
  
  a <- fviz_pca_biplot(pca.res.bray,
                       col.ind = meta_mean$imm_time,
                       addEllipses = TRUE,
                       ellipse.type = "convex",
                       col.var = "purple")
  
  
  b <- fviz_pca_biplot(pca.res.bray,
                  col.ind = meta_mean$imm_season,
                  addEllipses = TRUE,
                  ellipse.type = "convex",
                  col.var = "purple")
  
  
  c <- fviz_pca_biplot(pca.res.bray, 
                  col.ind = meta_mean$ret_season,
                  addEllipses = TRUE,
                  ellipse.type = "convex",
                  col.var = "purple")
  
  meta_mean$comb <- paste(meta_mean$imm_time, meta_mean$imm_season)
  
  
  d <- fviz_pca_biplot(pca.res.bray,
                  col.ind = meta_mean$comb,
                  addEllipses = TRUE,
                  ellipse.type = "convex",
                  col.var = "purple")
  
  
  # plot_pca_bray <- ggplot2::autoplot(pca.res.bray,
  #                                    data = meta_mean,
  #                                    colour = "imm_time",
  #                                    shape = "imm_season",
  #                                    label = TRUE,
  #                                    loadings.label.vjust = 14) + 
  #                  labs(title = "Bray-Curtis distance")
  # plot_pca_bray
 
  fin <- cowplot::plot_grid(a,
                            b,
                            c,
                            d,
                            ncol = 2,
                            nrow = 2)
  
  path_to_PCA_pool <- paste0("outputs/PCA_pool.pdf")
  ggsave(filename =  path_to_PCA_pool, plot = fin, width = 15, height = 12)
  
  
  return(path_to_PCA)
}
  
  