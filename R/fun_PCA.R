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
  
  
  
  #### With species pool ----------
  #### Load data and meta data ========
  data_pool <- read.csv(data_mean_pool, header = TRUE, row.names = "X")
  mat.bray <- vegan::vegdist(data_pool, method = "bray")
  
  pca.res.bray <- prcomp(data_pool,
                         center = TRUE,
                         scale. = TRUE)
  ?prcomp
  library(factoextra)
  a <- fviz_pca_biplot(pca.res.bray,
                       col.ind = meta_mean$imm_time,
                       addEllipses = TRUE,
                       ellipse.type = "confidence",
                       col.var = "purple")
  
  
  b <- fviz_pca_biplot(pca.res.bray,
                  col.ind = meta_mean$imm_season,
                  addEllipses = TRUE,
                  ellipse.type = "confidence",
                  col.var = "purple")
  
  
  c <- fviz_pca_biplot(pca.res.bray, 
                  col.ind = meta_mean$ret_season,
                  addEllipses = TRUE,
                  ellipse.type = "confidence",
                  col.var = "purple")
  
  meta_mean$comb <- paste(meta_mean$imm_time, meta_mean$imm_season)
  
  
  d <- fviz_pca_biplot(pca.res.bray,
                  col.ind = meta_mean$comb,
                  addEllipses = TRUE,
                  ellipse.type = "confidence",
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
  
  path_to_PCA_pool <- paste0("outputs/beta/PCA_pool.pdf")
  ggsave(filename =  path_to_PCA_pool, plot = fin, width = 10, height = 8)
  
  
  return(path_to_PCA)
}
  
  