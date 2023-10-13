##**************************************************#
## R script for functions in SPARRAfairness      ####
##**************************************************#
##
## James Liley, Ioanna Thoma
## March 2023
##

## Use command
##  R CMD build --resave-data SPARRAfairness
## to build package

##**************************************************#
## Packages and scripts                          ####
##**************************************************#

#' @import stats
#' @import graphics
#' @import grDevices
#' @import matrixStats
#' @import ranger
#' @import ggplot2
#' @import mvtnorm
#' @import cvAUC
#' @import ggrepel
#' @import patchwork
#' @import scales

require("stats")        # Base stats package
require("graphics")     # Base graphics package
require("grDevices")    # Base graphics devices package
require("matrixStats")  # Matrix row- and column-wise operations
require("ranger")       # Random forests

##**************************************************#
## Data documentation                            ####
##**************************************************#

##' @title All data for fairness measures
##' @description This object contains all data from analysis of fairness measures in SPARRA v3 and v4.
##' @docType data
##' @keywords data fairness sparra
"all_data"

##' @title Decomposition matrix
##' @description Matrix giving frequency of admission types for various groups at various score thresholds. Row names are of the form vX_Y_qZ, where X is version (3 or 4), Y is cohort (e.g., all, over 65, island postcode) and Z is quantile (1-20) of score. Column names are cause of admission or cause of death.
##' @docType data
##' @keywords data fairness sparra
"decomposition_matrix"

