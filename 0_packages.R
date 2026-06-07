###############################################################################
############################### 0 : PACKAGES ###################################
###############################################################################

# Packages used in the project
required_packages <- c(
  "deSolve",
  "ggplot2",
  "reshape2",
  "dplyr",
  "patchwork",
  "gtable",
  "gridExtra",
  "pbapply",
  "tidyr",
  "sensitivity"
)

# Install only missing packages
missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages) > 0) {
  install.packages(missing_packages, repos = "https://cran.rstudio.com")
}

# Load packages
library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)
library(patchwork)
library(gtable)
library(gridExtra)
library(grid)
library(parallel)
library(pbapply)
library(tidyr)
library(sensitivity)
