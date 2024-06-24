#' CINARMS: A Research Compendium
#' 
#' @description 
#' A data analysis pipeline for CINARMS
#' 
#' @author Baptiste Frattini \email{baptiste.frattini22@gmail.com}
#' 
#' @date 2023/04/12



## Install Dependencies (not available on CRAN) ----
# install.packages("~/R-projects/tar.gz files/adespatial_0.3-21.tar.gz", repos = NULL, type = "source")


# dependences management

renv::init()
renv::install()
renv::status()
renv::snapshot()


# make the pipeline
targets::tar_visnetwork()
targets::tar_make()

targets::tar_visnetwork()



