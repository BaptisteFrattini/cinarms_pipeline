#' PCoA microhab
#'
#' @param metadata_data_mean data
#' @return the path to the subseted raw data file
#' @export

fun_PCoA_microhab <- function(meta_data){
  
  # Import data and packages ####
  
  # packages
  
  remotes::install_deps(upgrade = "never")
  devtools::load_all()
  
  # data
  
  # meta_data = targets::tar_read(metadata_data)
  library(dplyr)
  library(ggplot2)
  meta = read.csv(meta_data[grepl("metadata", meta_data)], header = TRUE)
  data = read.csv(meta_data[!grepl("metadata", meta_data)], header = TRUE)
  data <- data[ , colSums(data) != 0]
  meta <- meta %>%
    mutate(o_c = recode(o_c, "F" = "C"))
  
  rownames(data) <- paste0(meta$arms_name,meta$Orientation,meta$o_c,meta$plate_number)
  
  meta$microhab <- paste0(meta$Orientation, "_", meta$o_c)
  
  
  data_subset <- data %>%
    filter(meta$microhab == "D_C" | meta$microhab == "U_O")

  meta_subset <- meta %>%
    filter(meta$microhab == "D_C" | meta$microhab == "U_O")
  
  # PCOA a 6 mois ####
  data_subset_6m <- data_subset %>%
    filter(meta_subset$imm_time == "6m")
  
  meta_subset_6m <- meta_subset %>%
    filter(meta_subset$imm_time == "6m")
  
  df_bray <- vegan::vegdist(data_subset_6m, method = "bray")
  pcoa_res_score <- ape::pcoa(df_bray)
  score <- pcoa_res_score$values$Rel_corr_eig[c(1,2)]
  pcoa_res <- cmdscale(df_bray)
  sites.scores <- pcoa_res[,c(1,2)]
  
  species.scores <- BiodiversityR::add.spec.scores(pcoa_res, data_subset_6m, method="cor.scores", multi=1, Rscale=F, scaling="1")
  species.scores <- species.scores$cpro
  nrow(species.scores)
  
  sim_imm_tim <- summary(vegan::simper(data_subset_6m, meta_subset_6m$microhab, permutations = 999))
  sim_imm_tim_1 <- sim_imm_tim[[1]]
  contrib_imm_tim_1 <- sim_imm_tim_1[(sim_imm_tim_1$cumsum < 0.95),]

  
  sel <- levels(as.factor(rownames(contrib_imm_tim_1)))
  
  species.scores <- species.scores[sel,]
  species.scores <- species.scores*0.50
  
  colnames(sites.scores) <- c("Dim1", "Dim2")
  
  biplot_6m <- ggplot() +
    geom_point(data = sites.scores, aes(x = Dim1, y = Dim2, color = meta_subset_6m$microhab, shape = meta_subset_6m$microhab), size = 3) +
    scale_color_manual(values = c("D_C" = "#ab49ab", "U_O" = "#d9a6d9")) +  # Adjust colors as needed
    geom_segment(data = species.scores, aes(x = 0, y = 0, xend = Dim1, yend = Dim2), arrow = arrow(length = unit(0.03, "npc"))) +
    ggrepel::geom_text_repel(data = species.scores, aes(x = Dim1, y = Dim2, label = rownames(species.scores)), box.padding = 0.5, max.overlaps = Inf) +
    # ggrepel::geom_text_repel(data = sites.scores, aes(x = Dim1, y = Dim2, label = rownames(sites.scores), color = meta_subset$imm_time, fontface = "bold"), vjust = -1.5) +
    labs(x = paste0("Dimension 1 (", 100*round(score[1], 3),"%)"), y = paste0("Dimension 2 (", 100*round(score[2], 3),"%)")) +
    theme_minimal()
  
  path_to_PCoA_6m <- paste0("outputs/PCA/PCoA_DC_UO_6m.png")
  ggsave(filename =  path_to_PCoA_6m, plot = biplot_6m, width = 8, height = 5)
  
  biplot_6m_b <- ggplot() +
    geom_point(data = sites.scores, aes(x = Dim1, y = Dim2, color = meta_subset_6m$microhab, shape = meta_subset_6m$microhab), size = 3) +
    scale_color_manual(values = c("D_C" = "#ab49ab", "U_O" = "#d9a6d9")) +  # Adjust colors as needed
    geom_segment(data = species.scores, aes(x = 0, y = 0, xend = Dim1, yend = Dim2), arrow = arrow(length = unit(0.03, "npc"))) +
    # ggrepel::geom_text_repel(data = species.scores, aes(x = Dim1, y = Dim2, label = rownames(species.scores)), box.padding = 0.5, max.overlaps = Inf) +
    # ggrepel::geom_text_repel(data = sites.scores, aes(x = Dim1, y = Dim2, label = rownames(sites.scores), color = meta_subset$imm_time, fontface = "bold"), vjust = -1.5) +
    labs(x = paste0("Dimension 1 (", 100*round(score[1], 3),"%)"), y = paste0("Dimension 2 (", 100*round(score[2], 3),"%)")) +
    theme_minimal()
  
  path_to_PCoA_6m_blank <- paste0("outputs/PCA/PCoA_DC_UO_6m_blank.png")
  ggsave(filename =  path_to_PCoA_6m_blank, plot = biplot_6m_b, width = 8, height = 5)
  
  # PCOA a 1y ####
  
  data_subset_1y <- data_subset %>%
    filter(meta_subset$imm_time == "1y")
  
  meta_subset_1y <- meta_subset %>%
    filter(meta_subset$imm_time == "1y")
  
  df_bray <- vegan::vegdist(data_subset_1y, method = "bray")
  pcoa_res_score <- ape::pcoa(df_bray)
  score <- pcoa_res_score$values$Rel_corr_eig[c(1,2)]
  pcoa_res <- cmdscale(df_bray)
  sites.scores <- pcoa_res[,c(1,2)]
  
  species.scores <- BiodiversityR::add.spec.scores(pcoa_res, data_subset_1y, method="cor.scores", multi=1, Rscale=F, scaling="1")
  species.scores <- species.scores$cpro
  nrow(species.scores)
  
  sim_imm_tim <- summary(vegan::simper(data_subset_1y, meta_subset_1y$microhab, permutations = 999))
  sim_imm_tim_1 <- sim_imm_tim[[1]]
  contrib_imm_tim_1 <- sim_imm_tim_1[(sim_imm_tim_1$cumsum < 0.95),]
  
  
  sel <- levels(as.factor(rownames(contrib_imm_tim_1)))
  
  species.scores <- species.scores[sel,]
  species.scores <- species.scores*0.50
  
  colnames(sites.scores) <- c("Dim1", "Dim2")
  
  biplot_1y <- ggplot() +
    geom_point(data = sites.scores, aes(x = Dim1, y = Dim2, color = meta_subset_1y$microhab, shape = meta_subset_1y$microhab), size = 3) +
    scale_color_manual(values = c("D_C" = "#20745b", "U_O" = "#32b58e")) +  # Adjust colors as needed
    geom_segment(data = species.scores, aes(x = 0, y = 0, xend = Dim1, yend = Dim2), arrow = arrow(length = unit(0.03, "npc"))) +
    ggrepel::geom_text_repel(data = species.scores, aes(x = Dim1, y = Dim2, label = rownames(species.scores)), box.padding = 0.5, max.overlaps = Inf) +
    # ggrepel::geom_text_repel(data = sites.scores, aes(x = Dim1, y = Dim2, label = rownames(sites.scores), color = meta_subset$imm_time, fontface = "bold"), vjust = -1.5) +
    labs(x = paste0("Dimension 1 (", 100*round(score[1], 3),"%)"), y = paste0("Dimension 2 (", 100*round(score[2], 3),"%)")) +
    theme_minimal()
  
  path_to_PCoA_1y <- paste0("outputs/PCA/PCoA_DC_UO_1y.png")
  ggsave(filename =  path_to_PCoA_1y, plot = biplot_1y, width = 8, height = 5)
  
  biplot_1y_b <- ggplot() +
    geom_point(data = sites.scores, aes(x = Dim1, y = Dim2, color = meta_subset_1y$microhab, shape = meta_subset_1y$microhab), size = 3) +
    scale_color_manual(values = c("D_C" = "#20745b", "U_O" = "#32b58e")) +  # Adjust colors as needed
    geom_segment(data = species.scores, aes(x = 0, y = 0, xend = Dim1, yend = Dim2), arrow = arrow(length = unit(0.03, "npc"))) +
    # ggrepel::geom_text_repel(data = species.scores, aes(x = Dim1, y = Dim2, label = rownames(species.scores)), box.padding = 0.5, max.overlaps = Inf) +
    # ggrepel::geom_text_repel(data = sites.scores, aes(x = Dim1, y = Dim2, label = rownames(sites.scores), color = meta_subset$imm_time, fontface = "bold"), vjust = -1.5) +
    labs(x = paste0("Dimension 1 (", 100*round(score[1], 3),"%)"), y = paste0("Dimension 2 (", 100*round(score[2], 3),"%)")) +
    theme_minimal()
  
  path_to_PCoA_1y_blank <- paste0("outputs/PCA/PCoA_DC_UO_1y_blank.png")
  ggsave(filename =  path_to_PCoA_1y_blank, plot = biplot_1y_b, width = 8, height = 5)
  
  # PCOA a 2y ####
  
  data_subset_2y <- data_subset %>%
    filter(meta_subset$imm_time == "2y")
  
  meta_subset_2y <- meta_subset %>%
    filter(meta_subset$imm_time == "2y")
  
  df_bray <- vegan::vegdist(data_subset_2y, method = "bray")
  pcoa_res_score <- ape::pcoa(df_bray)
  score <- pcoa_res_score$values$Rel_corr_eig[c(1,2)]
  pcoa_res <- cmdscale(df_bray)
  sites.scores <- pcoa_res[,c(1,2)]
  
  species.scores <- BiodiversityR::add.spec.scores(pcoa_res, data_subset_2y, method="cor.scores", multi=1, Rscale=F, scaling="1")
  species.scores <- species.scores$cpro
  nrow(species.scores)
  
  sim_imm_tim <- summary(vegan::simper(data_subset_2y, meta_subset_2y$microhab, permutations = 999))
  sim_imm_tim_1 <- sim_imm_tim[[1]]
  contrib_imm_tim_1 <- sim_imm_tim_1[(sim_imm_tim_1$cumsum < 0.95),]
  
  
  sel <- levels(as.factor(rownames(contrib_imm_tim_1)))
  
  species.scores <- species.scores[sel,]
  species.scores <- species.scores*0.50
  
  colnames(sites.scores) <- c("Dim1", "Dim2")
  
  biplot_2y <- ggplot() +
    geom_point(data = sites.scores, aes(x = Dim1, y = Dim2, color = meta_subset_2y$microhab, shape = meta_subset_2y$microhab), size = 3) +
    scale_color_manual(values = c("D_C" = "#df6d14", "U_O" = "darkorange")) +  # Adjust colors as needed
    geom_segment(data = species.scores, aes(x = 0, y = 0, xend = Dim1, yend = Dim2), arrow = arrow(length = unit(0.03, "npc"))) +
    ggrepel::geom_text_repel(data = species.scores, aes(x = Dim1, y = Dim2, label = rownames(species.scores)), box.padding = 0.5, max.overlaps = Inf) +
    # ggrepel::geom_text_repel(data = sites.scores, aes(x = Dim1, y = Dim2, label = rownames(sites.scores), color = meta_subset$imm_time, fontface = "bold"), vjust = -1.5) +
    labs(x = paste0("Dimension 1 (", 100*round(score[1], 3),"%)"), y = paste0("Dimension 2 (", 100*round(score[2], 3),"%)")) +
    theme_minimal()
  
  path_to_PCoA_2y <- paste0("outputs/PCA/PCoA_DC_UO_2y.png")
  ggsave(filename =  path_to_PCoA_2y, plot = biplot_2y, width = 8, height = 5)
  
  biplot_2y_b <- ggplot() +
    geom_point(data = sites.scores, aes(x = Dim1, y = Dim2, color = meta_subset_2y$microhab, shape = meta_subset_2y$microhab), size = 3) +
    scale_color_manual(values = c("D_C" = "#df6d14", "U_O" = "darkorange")) +  # Adjust colors as needed
    geom_segment(data = species.scores, aes(x = 0, y = 0, xend = Dim1, yend = Dim2), arrow = arrow(length = unit(0.03, "npc"))) +
    # ggrepel::geom_text_repel(data = species.scores, aes(x = Dim1, y = Dim2, label = rownames(species.scores)), box.padding = 0.5, max.overlaps = Inf) +
    # ggrepel::geom_text_repel(data = sites.scores, aes(x = Dim1, y = Dim2, label = rownames(sites.scores), color = meta_subset$imm_time, fontface = "bold"), vjust = -1.5) +
    labs(x = paste0("Dimension 1 (", 100*round(score[1], 3),"%)"), y = paste0("Dimension 2 (", 100*round(score[2], 3),"%)")) +
    theme_minimal()
  
  path_to_PCoA_2y_blank <- paste0("outputs/PCA/PCoA_DC_UO_2y_blank.png")
  ggsave(filename =  path_to_PCoA_2y_blank, plot = biplot_2y_b, width = 8, height = 5)
 
  return(NULL) 
  
  
}

