#' Stacked column chart
#'
#' @param metadata_data_mean metadata
#' @param data_mean_pool data
#' @return the path to the plot generated
#' @export

fun_stacked_col_chart <- function(metadata_data_mean, data_mean_pool){
  
  # metadata_data_mean = targets::tar_read(mean_metadata_data)
  # data_mean_pool = targets::tar_read(data_pool_mean)

  meta_mean <- read.csv(metadata_data_mean[grepl("metadata", metadata_data_mean)], header = TRUE)
  data_pool <- read.csv(data_mean_pool, header = TRUE, row.names = "X")
  
  #### Stacked column chart ####
  
  data_pool$Ascidiacea <- rowSums(data.frame(data_pool$ascidiacea_s, data_pool$ascidiacea_c))
  Sediment <- data_pool$sediment
  Bare_plate <- data_pool$bare_plate
  
  Sediment = data_pool$sediment
  Bare_plate = data_pool$bare_plate
  
  data_pool <- data_pool[,-c(3,4,12,13)]
  
  
  sums <- colSums(data_pool)
  
  # Ordonner les indices des colonnes en fonction de leur somme
  ordered_indices <- order(sums, decreasing = FALSE)
  
  # Réarranger les colonnes de votre tableau en utilisant l'ordre obtenu
  data_pool <- data_pool[, ordered_indices]
  
  data_pool <- data.frame(cbind(Sediment, Bare_plate, data_pool))
  
  
  spe.T <- colnames(data_pool)
  spe.T <- make.names(spe.T)
  spe.T[c(2,4,6,9,10,11,12,13)] <- c("Bare plate","Other algae", "Porifera", "Bryozoa", "Prokariotic biotas", "Annelida", "Foraminifera", "Crustose coralline algae")
  data_pool$ordre <- c(1:15)
  
  data_pool <- data_pool[order(-data_pool$ordre),]
  data_pool <- data_pool[,-ncol(data_pool)]
  
  data_pool_t <- t(data_pool)
  row.names(data_pool_t) <- spe.T 
  data_pool <- t(data_pool_t)
 
  data_pool <- (data_pool/rowSums(data_pool))*100
  rowSums(data_pool)
  
  
  ?prop.table
  
  pcm <- reshape2::melt(data_pool, id = rownames(as.matrix(data_pool)))
  
  ncol(pcm)
  ####colours####
  col<-NULL
  col[1] <-  "#FAEFD1" 
  col[2] <-  "#899DA4"
  col[3] <-  "#F98400"
  col[4] <-  "#5BBCD6"
  col[5] <-  "#F2AD00"
  col[6] <-  "#FF0000"
  col[7] <-  "#FD6467"
  col[8] <-  "darkgreen"
  col[9] <-  "#FDD262"
  col[10] <- "#00A08A"
  col[11] <- "#446455"
  col[12] <- "#F1BB7B"
  col[13] <- "#C93312"
 
  
  library(ggplot2)
  mx <- ggplot(pcm, aes(x =Var1, fill =Var2, y = value)) + 
    geom_bar(stat = "identity", colour = "black") + 
    theme(axis.text.x = element_text(angle = 90, size = 14, colour = "black", vjust = 0.5, hjust = 1, face= "bold"), 
          axis.title.y = element_text(size = 16, face = "bold"), legend.title = element_text(size = 16, face = "bold"), 
          legend.text = element_text(size = 12, face = "bold", colour = "black"), 
          axis.text.y = element_text(colour = "black", size = 12, face = "bold")) + 
    scale_y_continuous(expand = c(0,0)) + 
    labs(x = "", y = "Average covers of ARMS plate faces (%)", fill = "Categories :") + 
    scale_fill_manual(values = col)
  mx = mx  + coord_flip()
  
  
  mx = mx + guides(fill = guide_legend(reverse = TRUE))

  path_to_stack_c_chart <-  here::here("outputs/stack_c_chart.pdf")
  ggsave(path_to_stack_c_chart, plot = mx, width = 8, height = 8.5)
  
  return(path_to_stack_c_chart)
  
  
}


