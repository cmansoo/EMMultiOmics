# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

hello <- function() {
  print("Hello, world!")
}

# setup
install.packages("roxygen2")
install.packages("devtools")
library(devtools)
library(roxygen2)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("cBioPortalData")
