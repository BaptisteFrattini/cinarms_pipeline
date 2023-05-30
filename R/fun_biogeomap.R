#' testing biogeonat
#'
#' @param metadata_data_mean data
#' @return the path to the subseted raw data file
#' @export

biogeo_test <- function(metadata_data_mean){
 
  # metadata_data_mean = targets::tar_read(mean_metadata_data)
  
  #### Load data and meta data ####
  
  df_mean <- read.csv(metadata_data_mean[!grepl("metadata", metadata_data_mean)], header = TRUE)
  meta_mean <- read.csv(metadata_data_mean[grepl("metadata", metadata_data_mean)], header = TRUE)
  
  #### step 1 ####
  data <- reshape2::melt(df_mean, id = rownames(as.matrix(df_mean)))
  nrow(data)
  arms <- rep(meta_mean$arms, 60)
  data$arms <- arms
  data
  data = data.frame(arms = as.factor(data$arms),
                    msp = as.factor(data$variable),
                    abun = data$value)
  data  
  
  #### step 2 ####
  library(biogeonetworks)
  writePajek(data, 
             site.field = "arms", # Name of your site column
             species.field = "msp", # Name of your species column
             filename = "data/derived-data/cina_biogeonet.net", # Name of the output file
             abundance.field = "abun") # (FACULTATIVE) Name of abundance column
  
  #### step 3 ####
  #{r mapequation, eval = TRUE, cache = TRUE, echo = TRUE, results = "hide"}
  system("infomap --undirected --tree --map -N 100 ./data/derived-data/cina_biogeonet.net ./outputs/biogeonetworkstest")
  
  #### step 4 ####
  cina.clusters <- readInfomapTree("outputs/biogeonetworkstest/cina_biogeonet.tree",
                                   network.summary = TRUE, # Prints a summary of the clustering results
                                   replace.leaf.names = TRUE, # Changes node numbers to actual names for terminal nodes (i.e. site & species names)
                                   db = data, # Specify original database here to know which nodes are sites and which nodes are species
                                   site.field = "arms",  # Site column in db
                                   species.field = "msp") # Species column in db 
  
  head(cina.clusters)
  
  plyr::count(cina.clusters$lvl1)
  
  cina.sites <- getSiteTable(data, # Your bipartite data.frame of STEP 1
                             site.field = "arms", # Name of your site column
                             network = cina.clusters) # Output table from Map Equation
  writeGDF(db = data, # Database of step 1
           site.field = "arms", # Name of site column in your database
           species.field = "msp", # Name of species column in your database
           abundance.field = "abun",
           network = cina.clusters, 
           filename = "outputs/biogeonetworkstest/cina_biogeo.gdf") # Name of the color field in the network
  
  ?writeGDF
  return()
  }  
    