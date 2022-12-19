#' beta_diversity_arms: A Research Compendium
#' 
#' @description 
#' Computing beta diversity composition
#' 
#' @author Baptiste Frattini \email{baptiste.frattini22@gmail.com}
#' 
#' @date 2022/12/16



## Install Dependencies (listed in DESCRIPTION) ----


## Load Project Addins (R Functions and Packages) ----

devtools::load_all(here::here())
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
