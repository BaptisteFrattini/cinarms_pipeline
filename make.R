#' CINARMS: A Research Compendium
#' 
#' @description 
#' A data analysis pipeline for CINARMS
#' 
#' @author Baptiste Frattini \email{baptiste.frattini22@gmail.com}
#' 
#' @date 2023/04/12



## Install Dependencies (listed in DESCRIPTION) ----



# dependences management

renv::init()
renv::install()
renv::status()
renv::snapshot()


# make the pipeline
targets::tar_visnetwork()
targets::tar_make()
targets::tar_visnetwork()



## Global Variables ----

# You can list global variables here (or in a separate R script)


## Run Project ----

# List all R scripts in a sequential order and using the following form:
# source(here::here("analyses", "script_X.R"))
