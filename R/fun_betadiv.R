#' beta diversity decomposition
#'
#' @param metadata_data_mean data
#' @return the path to the subseted raw data file
#' @export

beta_div_decomp<- function(metadata_data_mean){
  
  # metadata_data_mean = targets::tar_read(mean_metadata_data)
  
  #### Load data and meta data ####
  
  df_mean <- read.csv(metadata_data_mean[!grepl("metadata", metadata_data_mean)], header = TRUE)
  meta_mean <- read.csv(metadata_data_mean[grepl("metadata", metadata_data_mean)], header = TRUE)
  
  matrix.pa <- vegan::decostand(df_mean, "pa")
  rownames(matrix.pa) <- meta_mean$arms
  colnames(matrix.pa) <- meta_mean$arms
  
  B.pair.pa <- betapart::beta.pair(matrix.pa, index.family = "jaccard")
  mat.turn <- B.pair.pa$beta.jtu
  mat.nest <- B.pair.pa$beta.jne
  mat.jacc <- 1-B.pair.pa$beta.jac
  
  library(reshape2)
  
  #### inter/intra ####
  
  mat.turn
  df.turn <- melt(as.matrix(mat.turn), varnames = c("row", "col"))
  df.turn <- subset(df.turn, row != col)
  
  df.turn$row <- substr(df.turn$row, 1, 5)
  df.turn$col <- substr(df.turn$col, 1, 5)
  df.turn$same_value <- ifelse(df.turn$row == df.turn$col, "Yes", "No")
  
  df.nest <- melt(as.matrix(mat.nest), varnames = c("row", "col"))
  df.nest <- subset(df.nest, row != col)
  
  df.nest$row <- substr(df.nest$row, 1, 5)
  df.nest$col <- substr(df.nest$col, 1, 5)
  df.nest$same_value <- ifelse(df.nest$row == df.nest$col, "Yes", "No")
  
  
  df.jacc <- melt(as.matrix(mat.jacc), varnames = c("row", "col"))
  df.jacc <- subset(df.jacc, row != col)
  
  df.jacc$row <- substr(df.jacc$row, 1, 5)
  df.jacc$col <- substr(df.jacc$col, 1, 5)
  df.jacc$same_value <- ifelse(df.jacc$row == df.jacc$col, "Yes", "No")
  
  
  decomp.pa <- as.data.frame(cbind(df.nest$same_value, df.turn$value, df.nest$value  ,df.jacc$value))
  
  colnames(decomp.pa) <- c("intrasite","Turnover", "Nestedness", "Jaccard_diss")
  
  decomp.pa$Jaccard_diss <- abs(as.numeric(decomp.pa$Jaccard_diss))
  decomp.pa$Nestedness <- abs(as.numeric(decomp.pa$Nestedness))
  decomp.pa$Turnover <- abs(as.numeric(decomp.pa$Turnover))
  library(ggplot2)
  library(ggtern)
   
  #   trip <- ggtern(decomp.pa, aes(x = Nestedness, y =Jaccard_diss, z = Turnover, color = intrasite)) +
  #     # Ajouter une couche de points
  #     geom_point(size = 0.5) +
  #     # Ajouter la densité des points
  #     stat_density_tern(alpha = 0.5) +
  #     # Ajouter un titre
  #     labs(title = "Graphique ternaire de similarité écologique") +
  #     # Ajouter des étiquettes d'axes
  #     xlab("Nestedness") +
  #     ylab("Jaccard similarity") +
  #     zlab("Turnover") +
  #     tern_limits(T=0.98, L=0.6, R=0.6) +
  #     # Ajouter une couche de couleur pour les comparaisons intra-site et inter-site
  #     scale_color_manual(values = c("#00A08A", "#FD6467"), labels = c("Inter-site", "Intra-site"))
  # trip
  # 
  mean.tur.intra <- tapply(decomp.pa$Turnover, decomp.pa$intrasite, mean)
  sd.tur.intra <- tapply(decomp.pa$Turnover, decomp.pa$intrasite, sd)
  res.t.test <- t.test(decomp.pa$Turnover ~ decomp.pa$intrasite)
  stat <- res.t.test$statistic
  p <- res.t.test$p.value
  v <- ggplot(decomp.pa, aes(x=intrasite, y=Turnover)) +
    geom_boxplot() +
    xlab(" ") +
    ylab("Turnover component index") +
    scale_x_discrete(labels=c("Among sites comparisons","Within sites comparisons")) +
    annotate(geom = "text",  x = 1.8, y = 0.6, label = paste0("t = 7.19 ; p < 0.05"),
             color = "black", size = 5, hjust = 0, vjust = 1) +
    theme(axis.text=element_text(size=14))

   path_to_turncomp <- paste0("outputs/beta/turnover_comp", arms_id,".pdf")
   ggsave(filename = path_to_turncomp, plot = v, width = 6, height = 6)
   
   #### 6m - 1a - 2a ####
   
   mat.turn
   df.turn <- melt(as.matrix(mat.turn), varnames = c("row", "col"))
   df.turn <- subset(df.turn, row != col)
   
   df.turn$row <- substr(df.turn$row, 1, 5)
   df.turn$col <- substr(df.turn$col, 1, 5)
   df.turn$same_value <- ifelse(df.turn$row == df.turn$col, "Yes", "No")
  
   
}





