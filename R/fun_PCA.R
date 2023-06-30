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
  #### PCA ####
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
                       select.var = list(contrib = 35),
                       col.var = "purple",
                       repel = TRUE,
                       pointsize = 3,
                       labelsize = 5)
  
  meta_mean$set <- paste0(meta_mean$imm_time,"_",meta_mean$imm_season)
  fviz_screeplot(pca.res.full)
  y <- fviz_pca_biplot(pca.res.full,
                       col.ind = meta_mean$set,
                       addEllipses = TRUE,
                       select.var = list(contrib = 20),
                       ellipse.type = "convex",
                       col.var = "black",
                       repel = TRUE,
                       pointsize = 2,
                       labelsize = 5,
                       palette = c("firebrick3","firebrick3","dodgerblue3","forestgreen","forestgreen") ) 
  # + scale_color_manual(values=c("firebrick3","firebrick3","dodgerblue3","forestgreen","forestgreen"))
  
  #same but with bray curtis distance
  bc <- vegan::vegdist(df_mean, method = "bray")

  pca.res.bc <- prcomp(bc,
                       center = TRUE,
                       scale. = TRUE)
  ybray <- fviz_pca_ind(pca.res.bc,
                        col.ind = meta_mean$set,
                        addEllipses = TRUE,
                        ellipse.type = "convex",
                        col.var = "black",
                        repel = TRUE,
                        pointsize = 2,
                        labelsize = 5,
                        palette = c("firebrick3","firebrick3","dodgerblue3","forestgreen","forestgreen"))
  
  
  
  meta_mean$set2 <- paste0(meta_mean$imm_time,"_",meta_mean$ret_season)
  
  x <- fviz_pca_biplot(pca.res.full,
                       col.ind = meta_mean$set2,
                       addEllipses = TRUE,
                       ellipse.type = "convex",
                       select.var = list(contrib =35),
                       col.var = "black",
                       repel = TRUE,
                       pointsize = 2,
                       labelsize = 5,
                       palette = c("firebrick3","firebrick3","dodgerblue3","forestgreen","forestgreen") ) 
  # + scale_color_manual(values=c("firebrick3","firebrick3","dodgerblue3","forestgreen","forestgreen"))
  
  #blank
  
  v <- fviz_pca_biplot(pca.res.full,
                       col.ind = meta_mean$set,
                       addEllipses = TRUE,
                       ellipse.type = "convex",
                       select.var = list(contrib = 35),
                       geom.var = "arrow",
                       col.var = "black",
                       repel = TRUE,
                       pointsize = 2,
                       labelsize = 5,
                       palette = c("firebrick3","firebrick3","dodgerblue3","forestgreen","forestgreen") ) 

  i <- fviz_pca_biplot(pca.res.full,
                       col.ind = meta_mean$set2,
                       addEllipses = TRUE,
                       ellipse.type = "convex",
                       select.var = list(contrib = 35),
                       geom.var = "arrow",
                       col.var = "black",
                       repel = TRUE,
                       pointsize = 2,
                       labelsize = 5,
                       palette = c("firebrick3","firebrick3","dodgerblue3","forestgreen","forestgreen") ) 
  
  path_to_PCA_select_set <- paste0("outputs/PCA/PCA_select_set.pdf")
  ggsave(filename =  path_to_PCA_select_set, plot = y, width = 12, height = 10)
  
  path_to_PCA_select_set2 <- paste0("outputs/PCA/PCA_select_set2.pdf")
  ggsave(filename =  path_to_PCA_select_set2, plot = x, width = 12, height = 10)
  
  path_to_PCA_select <- paste0("outputs/PCA/PCA_select.pdf")
  ggsave(filename =  path_to_PCA_select, plot = w, width = 12, height = 10)
  
  path_to_PCA_select_set_blank <- paste0("outputs/PCA/PCA_select_set_blank.pdf")
  ggsave(filename = path_to_PCA_select_set_blank , plot = v, width = 12, height = 10)
  
  path_to_PCA_select_set2_blank <- paste0("outputs/PCA/PCA_select_set2_blank.pdf")
  ggsave(filename =  path_to_PCA_select_set2_blank, plot = i, width = 12, height = 10)
  
  path_to_PCA_select_set_bray <- paste0("outputs/PCA/PCA_select_set_bray.pdf")
  ggsave(filename =  path_to_PCA_select_set_bray, plot = ybray, width = 12, height = 10)
  
  
  
  #### PCoA ####
  
  df_bray <- vegan::vegdist(df_mean, method = "bray")
  pcoa_res <- ape::pcoa(df_bray)
  sites.scores <- pcoa_res$vectors[,c(1,2)]
  
  library(BiodiversityR)
  pcoa_res2 <- cmdscale(df_bray)
  species.scores <- BiodiversityR::add.spec.scores(pcoa_res2, df_mean, method="cor.scores", multi=1, Rscale=F, scaling="1")
  species.scores <- species.scores$cproj
  
  sim_imm_tim <- summary(vegan::simper(df_mean, meta_mean$imm_time))
  
  ##1
  sim_imm_tim_1 <- sim_imm_tim[[1]]
  contrib_imm_tim_1 <- sim_imm_tim_1[(sim_imm_tim_1$p < 0.05) & (sim_imm_tim_1$average > 0.0005),]
  
  sim_imm_tim_2 <- sim_imm_tim[[2]]
  contrib_imm_tim_2 <- sim_imm_tim_2[(sim_imm_tim_2$p < 0.05) & (sim_imm_tim_2$average > 0.0005),]
  
  sim_imm_tim_3 <- sim_imm_tim[[3]]
  contrib_imm_tim_3 <- sim_imm_tim_3[(sim_imm_tim_3$p < 0.05) & (sim_imm_tim_3$average > 0.0005),]
  
  
  ##2
  sim_imm_tim_1 <- sim_imm_tim[[1]]
  contrib_imm_tim_1 <- sim_imm_tim_1[(sim_imm_tim_1$p < 0.05),]
  contrib_imm_tim_1bis <- sim_imm_tim_1[1:10,])
  sim_imm_tim_2 <- sim_imm_tim[[2]]
  contrib_imm_tim_2 <- sim_imm_tim_2[(sim_imm_tim_2$p < 0.05) | (sim_imm_tim_2[1:10,]),]
  
  sim_imm_tim_3 <- sim_imm_tim[[3]]
  contrib_imm_tim_3 <- sim_imm_tim_3[(sim_imm_tim_3$p < 0.05) | (sim_imm_tim_3[1:10,]),]
  
  
  
  sel <- levels(as.factor(c(rownames(contrib_imm_tim_1), rownames(contrib_imm_tim_2), rownames(contrib_imm_tim_3))))
  
  species.scores <- species.scores[sel,]
  species.scores <- species.scores*0.2
  
  colnames(sites.scores) <- c("Dim1", "Dim2")
  

  biplot <- ggplot() +
    geom_point(data = sites.scores, aes(x = Dim1, y = Dim2, color = meta_mean$imm_time, shape = meta_mean$imm_season, size = 3)) +
    geom_segment(data = species.scores, aes(x = 0, y = 0, xend = Dim1, yend = Dim2), arrow = arrow(length = unit(0.03, "npc"))) +
    ggrepel::geom_text_repel(data = species.scores, aes(x = Dim1, y = Dim2, label = rownames(species.scores)), box.padding = 0.5, max.overlaps = Inf) +
    ggrepel::geom_text_repel(data = sites.scores, aes(x = Dim1, y = Dim2, label = rownames(sites.scores), color = meta_mean$imm_time, fontface = "bold"), vjust = -1.5) +
    labs(x = "Dimension 1", y = "Dimension 2") +
    theme_minimal()
  
  # Display the biplot
  print(biplot)
  
  
  pcoa.res <- ape::pcoa(df_bray)
  ?stats::princomp

  autoplot(pcoa.res)
  
  #### With species pool ----------
  #### Load data and meta data ========
  data_pool <- read.csv(data_mean_pool, header = TRUE, row.names = "X")
  mat.bray <- vegan::vegdist(data_pool, method = "bray")
  
  pca.res.bray <- prcomp(data_pool,
                         center = TRUE,
                         scale. = TRUE)

  
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
  
  d2 <- fviz_pca_biplot(pca.res.bray,
                  col.ind = meta_mean$comb,
                  addEllipses = TRUE,
                  ellipse.type = "convex",
                  # geom.var = "arrow",
                  col.var = "black",
                  repel = TRUE,
                  pointsize = 2,
                  labelsize = 5,
                  palette = c("firebrick3","firebrick3","dodgerblue3","forestgreen","forestgreen") ) 
  
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
  
  path_to_PCA_pool <- paste0("outputs/PCA/PCA_poool.pdf")
  ggsave(filename =  path_to_PCA_pool, plot = fin, width = 12, height = 10)
  
  
  path_to_PCA_set_pool <- paste0("outputs/PCA/PCA_set_pool.pdf")
  ggsave(filename =  path_to_PCA_set_pool, plot = d2, width = 12, height = 10)
  

  
  return(path_to_PCA)
}
  
  