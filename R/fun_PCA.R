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
  library(plotly)
  library(ggfortify)
  library(factoextra)
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
  # ape::biplot.pcoa(pcoa_res, df_mean, dir.axis2 = -1)
  
  sites.scores <- pcoa_res$vectors[,c(1,2)]

  
  library(BiodiversityR)
  pcoa_res2 <- cmdscale(df_bray)
  species.scores <- BiodiversityR::add.spec.scores(pcoa_res2 , df_mean, method="cor.scores", multi=1, Rscale=F, scaling="1")
  species.scores <- species.scores$cproj
  nrow(species.scores)
  
  sim_imm_tim <- summary(vegan::simper(df_mean, meta_mean$imm_time, permutations = 9999))
##1
  # sim_imm_tim_1 <- sim_imm_tim[[1]]
  # contrib_imm_tim_1 <- sim_imm_tim_1[(sim_imm_tim_1$p < 0.0011) | (sim_imm_tim_1$average > 0.01),]
  # 
  # sim_imm_tim_2 <- sim_imm_tim[[2]]
  # contrib_imm_tim_2 <- sim_imm_tim_2[(sim_imm_tim_2$p < 0.0011) | (sim_imm_tim_2$average > 0.01),]
  # 
  # sim_imm_tim_3 <- sim_imm_tim[[3]]
  # contrib_imm_tim_3 <- sim_imm_tim_3[(sim_imm_tim_3$p < 0.0011) | (sim_imm_tim_3$average > 0.01),]
  # 
  
  ##2
  sim_imm_tim_1 <- sim_imm_tim[[1]]
  contrib_imm_tim_1 <- sim_imm_tim_1[(sim_imm_tim_1$cumsum < 0.95),]
  
  sim_imm_tim_2 <- sim_imm_tim[[2]]
  contrib_imm_tim_2 <- sim_imm_tim_2[(sim_imm_tim_2$cumsum < 0.95),]
  
  sim_imm_tim_3 <- sim_imm_tim[[3]]
  contrib_imm_tim_3 <- sim_imm_tim_3[(sim_imm_tim_3$cumsum < 0.95),]

  
  
  sel <- levels(as.factor(c(rownames(contrib_imm_tim_1), rownames(contrib_imm_tim_2), rownames(contrib_imm_tim_3))))
  
  species.scores <- species.scores[sel,]
  species.scores <- species.scores*0.25
  #### labels ####
  # install.packages("ggtext")
  # library(ggtext)
  # label2 = c(
  #   "<i>Ascc</i> <i>Aplidium</i> sp",
  #   "Bare plate",
  #   "Bivalvia",
  #   "<i>Ascc</i> <i>Botryllus_gregalis</i>",
  #   "<i>Ascc</i> <i>Botryllus_tuberatus</i>",
  #   "Brown_erect_macroalgae",
  #   "Calcareous_worm_tube",
  #   "Crustose_coralline_algae",
  #   "Cirripedia",
  #   "Cyanobacteria",
  #   "<i>Ascc</i> <i>Didemnidae</i> sp2.",
  #   "<i>Ascc</i> <i>Didemnidae</i> sp3.",
  #   "<i>Bry</i> <i>Disporella_novaehollandiae</i>",
  #   "<i>Ascc</i> <i>Eusynstyela_hartmeyeri</i>",
  #   "Green_biofilm",
  #   "Hydrozoa",
  #   "<i>For</i> <i>Miniacina</i> sp1. (adult)",
  #   "<i>For</i> <i>Miniacina</i> sp1. (juv)",
  #   "<i>For</i> <i>Miniacina</i> sp2.",
  #   "Por_msp11",
  #   "Por_msp12",
  #   "Ascc_msp15",
  #   "Por_msp15",
  #   "Por_msp3",
  #   "Por_msp4",
  #   "Ascs_msp5",
  #   "Bryo_msp5",
  #   "<i>Ascs</i> <i>Polycarpa</i> sp.",
  #   "Red_biofilm",
  #   "Red_erect_macroalgae",
  #   "Sediments",
  #   "Soft_worm_tube",
  #   "<i>Bry</i> <i>Tubulipora</i> sp.",
  #   "<i>Bry</i> <i>Watersipora_subtorquata</i>"
  # )
  # 
 #  labels1 <- as.character(c(substitute(paste("Ascc_",italic("Aplidium"),"_sp")),
 #               "Bare plate",
 #               "Bivalvia",
 #               substitute(paste("Ascc_",italic("Botryllus_gregalis"))),
 #               substitute(paste("Ascc_",italic("Botryllus_tuberatus"))),
 #               "Brown_erect_macroalgae",
 #               "Calcareous_worm_tube",
 #               "Crustose_coralline_algae",
 #               "Cirripedia",
 #               "Cyanobacteria",
 #               substitute(paste("Ascc_",italic("Didemnidae"),"_sp2.")),
 #               substitute(paste("Ascc_",italic("Didemnidae"),"_sp3.")),
 #               substitute(paste("Bry_",italic("Disporella_novaehollandiae"))),
 #               substitute(paste("Ascc_",italic("Eusynstyela_hartmeyeri"))),
 #               "Green_biofilm",
 #               "Hydrozoa",
 #               substitute(paste("For_",italic("Miniacina"),"_sp1. (adult)")),
 #               substitute(paste("For_",italic("Miniacina"),"_sp1. (juv)")),
 #               substitute(paste("For_",italic("Miniacina"),"_sp2.")),
 #               "Por_msp11",
 #               "Por_msp12",
 #               "Ascc_msp15",
 #               "Por_msp15",
 #               "Por_msp3",                       
 #               "Por_msp4",
 #               "Ascs_msp5",
 #               "Bryo_msp5",
 #               substitute(paste("Ascs_",italic("Polycarpa"),"_sp.")),
 #               "Red_biofilm",
 #               "Red_erect_macroalgae",
 #               "Sediments",
 #               "Soft_worm_tube",
 #               substitute(paste("Bry_",italic("Tubulipora"),"_sp.")),
 #               substitute(paste("Bry_",italic("Watersipora_subtorquata")))))
 #  
 # labels1 <- as.character(c(
 #    expression(Ascc[italic("Aplidium")]*" sp"),
 #    "Bare plate",
 #    "Bivalvia",
 #    expression(Ascc[italic("Botryllus_gregalis")]),
 #    expression(Ascc[italic("Botryllus_tuberatus")]),
 #    "Brown_erect_macroalgae",
 #    "Calcareous_worm_tube",
 #    "Crustose_coralline_algae",
 #    "Cirripedia",
 #    "Cyanobacteria",
 #    expression(Ascc[italic("Didemnidae")]*" sp2."),
 #    expression(Ascc[italic("Didemnidae")]*" sp3."),
 #    expression(Bry[italic("Disporella_novaehollandiae")]),
 #    expression(Ascc[italic("Eusynstyela_hartmeyeri")]),
 #    "Green_biofilm",
 #    "Hydrozoa",
 #    expression(For[italic("Miniacina")]*" sp1. (adult)"),
 #    expression(For[italic("Miniacina")]*" sp1. (juv)"),
 #    expression(For[italic("Miniacina")]*" sp2."),
 #    "Por_msp11",
 #    "Por_msp12",
 #    "Ascc_msp15",
 #    "Por_msp15",
 #    "Por_msp3",
 #    "Por_msp4",
 #    "Ascs_msp5",
 #    "Bryo_msp5",
 #    expression(Ascs[italic("Polycarpa")]*" sp."),
 #    "Red_biofilm",
 #    "Red_erect_macroalgae",
 #    "Sediments",
 #    "Soft_worm_tube",
 #    expression(Bry[italic("Tubulipora")]*" sp."),
 #    expression(Bry[italic("Watersipora_subtorquata")])))

  #### suite ####
  colnames(sites.scores) <- c("Dim1", "Dim2")

  biplot1 <- ggplot() +
    geom_point(data = sites.scores, aes(x = Dim1, y = Dim2, color = meta_mean$imm_time, shape = meta_mean$imm_season, size = 3)) +
    scale_color_manual(values = c("6m" = "#CC66CC", "1y" = "#1B9E77", "2y" = "#FF7F00")) +  # Adjust colors as needed
    geom_segment(data = species.scores, aes(x = 0, y = 0, xend = Dim1, yend = Dim2), arrow = arrow(length = unit(0.03, "npc"))) +
    ggrepel::geom_text_repel(data = species.scores, aes(x = Dim1, y = Dim2, label = rownames(species.scores)), box.padding = 0.5, max.overlaps = Inf) +
    ggrepel::geom_text_repel(data = sites.scores, aes(x = Dim1, y = Dim2, label = rownames(sites.scores), color = meta_mean$imm_time, fontface = "bold"), vjust = -1.5) +
    labs(x = "Dimension 1", y = "Dimension 2") +
    theme_minimal()
  
  sites.scores <- cbind(sites.scores, meta_mean)
  
  biplot1_test <- ggplot() +
    # Polygons per arm (transparent)
    geom_polygon(
      data = sites.scores,
      aes(x = Dim1, y = Dim2, group = arms_name, fill = imm_time),
      alpha = 0.2, color = "grey40", linewidth = 0.5
    ) +
    # Points
    geom_point(
      data = sites.scores,
      aes(x = Dim1, y = Dim2, color = imm_time, shape = imm_season),
      size = 3
    ) +
    # Site labels
    ggrepel::geom_text_repel(
      data = sites.scores,
      aes(x = Dim1, y = Dim2, label = rownames(sites.scores), color = imm_time),
      fontface = "bold",
      vjust = -1.5
    ) +
    scale_color_manual(values = c("6m" = "#CC66CC", "1y" = "#1B9E77", "2y" = "#FF7F00")) +
    scale_fill_manual(values = c("6m" = "#CC66CC", "1y" = "#1B9E77", "2y" = "#FF7F00")) +
    labs(x = "Dimension 1", y = "Dimension 2") +
    theme_minimal()
  
  path_to_PCoA_full <- paste0("outputs/PCA/PCoA_full_blank.pdf")
  ggsave(filename =  path_to_PCoA_full, plot = biplot1, width = 12, height = 10)
  
  
  biplot1_ret <- ggplot() +
    geom_point(data = sites.scores, aes(x = Dim1, y = Dim2, color = meta_mean$imm_time, shape = meta_mean$ret_season, size = 3)) +
    scale_color_manual(values = c("6m" = "#CC66CC", "1y" = "#1B9E77", "2y" = "#FF7F00")) +  # Adjust colors as needed
    geom_segment(data = species.scores, aes(x = 0, y = 0, xend = Dim1, yend = Dim2), arrow = arrow(length = unit(0.03, "npc"))) +
    # ggrepel::geom_text_repel(data = species.scores, aes(x = Dim1, y = Dim2, label = rownames(species.scores)), box.padding = 0.5, max.overlaps = Inf) +
    ggrepel::geom_text_repel(data = sites.scores, aes(x = Dim1, y = Dim2, label = rownames(sites.scores), color = meta_mean$imm_time, fontface = "bold"), vjust = -1.5) +
    labs(x = "Dimension 1", y = "Dimension 2") +
    theme_minimal()
  
  path_to_PCoA_full <- paste0("outputs/PCA/PCoA_full_blank_ret.pdf")
  ggsave(filename =  path_to_PCoA_full, plot = biplot1_ret, width = 12, height = 10)
  
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
                  palette = c("#CC66CC","#CC66CC","#1B9E77","#FF7F00","#FF7F00") ) 
  
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
  
  #### pcoa ####
  
  df_bray <- vegan::vegdist(data_pool, method = "bray")
  pcoa_res <- ape::pcoa(df_bray)
  sites.scores <- pcoa_res$vectors[,c(1,2)]
  
  pcoa_res2 <- cmdscale(df_bray)
  species.scores <- BiodiversityR::add.spec.scores(pcoa_res2 , data_pool, method="cor.scores", multi=1, Rscale=F, scaling="1")
  species.scores <- species.scores$cproj

  species.scores <- species.scores*0.25

   
  colnames(sites.scores) <- c("Dim1", "Dim2")
  
  label3 <- c("Porifera", "Bryozoa", "Ascidiacea (colonial)", "Ascidiacea (solitary)", "Foraminifera",      
  "Other algae", "Annelida", "Prokariotic biotas", "Bivalvia", "Hydrozoa",          
   "Cirripedia", "Bare plate", "Sediments", "CCA")
  
  biplot2 <- ggplot() +
    geom_point(data = sites.scores, aes(x = Dim1, y = Dim2, color = meta_mean$imm_time, shape = meta_mean$imm_season, size = 3)) +
    scale_color_manual(values = c("6m" = "#CC66CC", "1y" = "#1B9E77", "2y" = "#FF7F00")) +  # Adjust colors as needed
    geom_segment(data = species.scores, aes(x = 0, y = 0, xend = Dim1, yend = Dim2), arrow = arrow(length = unit(0.03, "npc"))) +
    ggrepel::geom_text_repel(data = species.scores, aes(x = Dim1, y = Dim2, label = label3), box.padding = 0.5, max.overlaps = Inf) +
    ggrepel::geom_text_repel(data = sites.scores, aes(x = Dim1, y = Dim2, label = rownames(sites.scores), color = meta_mean$imm_time, fontface = "bold"), vjust = -1.5) +
    labs(x = "Dimension 1", y = "Dimension 2") +
    theme_minimal()
  
  sites.scores <- cbind(sites.scores, meta_mean)
  
  biplot2_test <- ggplot() +
    # Polygons per arm (transparent)
    geom_polygon(
      data = sites.scores,
      aes(x = Dim1, y = Dim2, group = arms_name, fill = imm_time),
      alpha = 0.2, color = "grey40", linewidth = 0.5
    ) +
    # Points
    geom_point(
      data = sites.scores,
      aes(x = Dim1, y = Dim2, color = imm_time, shape = imm_season),
      size = 3
    ) +
    geom_segment(
      data = species.scores,
      aes(x = 0, y = 0, xend = Dim1, yend = Dim2),
      arrow = arrow(length = unit(0.03, "npc")),
      color = "black"
    ) +
    # Species labels (no arrows)
    ggrepel::geom_text_repel(
      data = species.scores,
      aes(x = Dim1, y = Dim2, label = rownames(species.scores)),
      box.padding = 0.5,
      max.overlaps = Inf
    ) +
    # Site labels
    ggrepel::geom_text_repel(
      data = sites.scores,
      aes(x = Dim1, y = Dim2, label = rownames(sites.scores), color = imm_time),
      fontface = "bold",
      vjust = -1.5
    ) +
    scale_color_manual(values = c("6m" = "#CC66CC", "1y" = "#1B9E77", "2y" = "#FF7F00")) +
    scale_fill_manual(values = c("6m" = "#CC66CC", "1y" = "#1B9E77", "2y" = "#FF7F00")) +
    labs(x = "Dimension 1", y = "Dimension 2") +
    theme_minimal()
  
  
  path_to_PCoA_pool <- paste0("outputs/PCA/PCoA_pool_blank.pdf")
  ggsave(filename =  path_to_PCoA_pool, plot = biplot2, width = 12, height = 10)
  
  
  biplot2_ret <- ggplot() +
    geom_point(data = sites.scores, aes(x = Dim1, y = Dim2, color = meta_mean$imm_time, shape = meta_mean$ret_season, size = 3)) +
    scale_color_manual(values = c("6m" = "#CC66CC", "1y" = "#1B9E77", "2y" = "#FF7F00")) +  # Adjust colors as needed
    geom_segment(data = species.scores, aes(x = 0, y = 0, xend = Dim1, yend = Dim2), arrow = arrow(length = unit(0.03, "npc"))) +
    # ggrepel::geom_text_repel(data = species.scores, aes(x = Dim1, y = Dim2, label = label3), box.padding = 0.5, max.overlaps = Inf) +
    ggrepel::geom_text_repel(data = sites.scores, aes(x = Dim1, y = Dim2, label = rownames(sites.scores), color = meta_mean$imm_time, fontface = "bold"), vjust = -1.5) +
    labs(x = "Dimension 1", y = "Dimension 2") +
    theme_minimal()
  
  path_to_PCoA_pool <- paste0("outputs/PCA/PCoA_pool_blank_ret.pdf")
  ggsave(filename =  path_to_PCoA_pool, plot = biplot2_ret, width = 12, height = 10)
  
  return(path_to_PCA)
}
  
  