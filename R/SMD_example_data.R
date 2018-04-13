#' Example data set for the package SurrogateMinimalDepth
#'
#' A dataset containing 100 individuals (rows), 1 dependent variable (first column, y), and 200 additional variables (columns 2 to 201):
#' For the simulation of the 200 additional variables nine variables X1,…,X9 called basic variables were simulated. X1 to X6 are causal
#' with constant effect sizes of 1 and X7 to X9 are noncausal with effect size of 0. The outcome y is a linear combination of the causal
#' predictor variables and a normally distributed error term. All basic variables were sampled from a normal distribution
#' (N(0,1)) just like the noise (N(0,0.2)). For each of the six basic variables X1, X2, X3, X7, X8, and X9 ten variables
#' with predefined correlations of 0.9 for X1 and X7, 0.6 for X2 and X8, and 0.3 for X3 and X9 were obtained by \link[WGCNA]{simulateModule} function of
#' the R package WGCNA. The ten variables of each basis variable are labeled: Cp_basicvariable_number. Additional non-correlated and
#' independent predictor variables (cgn) were simulated using the standard normal distribution to reach a total number of 200 variables.
#'
#'
#' @docType data
#' @keywords datasets
#' @name SMD_example_data
#' @usage data(SMD_example_data)
#' @format A matrix with 100 rows and 1001 columns
NULL
