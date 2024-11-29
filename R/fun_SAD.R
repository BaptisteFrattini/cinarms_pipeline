#' Species abundance distribution 
#'
#' @param metadata_data_mean data
#' @return the path to the ...
#' @export

fun_SAD <- function(metadata_data_mean){
 
  # metadata_data_mean = targets::tar_read(mean_metadata_data)
  df_mean <- read.csv(metadata_data_mean[!grepl("metadata", metadata_data_mean)], header = TRUE)
  
  library(vegan)
  library(sads)
  
  # Calculer l'abondance totale par espèce
  species_abundance <- as.numeric(colSums(df_mean))
  
  # Créer des classes d'abondance (octaves)
  octaves <- octav(species_abundance, preston = TRUE)
  ?octav
  # Afficher les classes d'abondance
  print(octaves)
  
  # Tracer la courbe de distribution d'abondance des espèces
  plot(octaves, main="Species Abundance Distribution (SAD)",
       xlab="Classes d'abondance (octave)",
       ylab="Nombre d'espèces",
       col="blue", border="black")
  
  mod <- fitsad(species_abundance,'lnorm')
  summary(mod)
  AIC(mod)
  plot(mod)
  
  #### return ####
  return(NULL) 
}
