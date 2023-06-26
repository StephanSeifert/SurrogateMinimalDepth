#' Variable selection with Surrogate Minimal Depth (SMD) (MAIN FUNCTION)
#'
#' This function executes SMD applying \link[ranger]{ranger} for random forests generation and a modified version of \link[rpart]{rpart} to find surrogate variables.
#'
#' @param x data.frame of predictor variables with variables in
#'   columns and samples in rows (Note: missing values are not allowed)
#' @param y vector with values of phenotype variable (Note: will be converted to factor if
#'   classification mode is used). For survival forests this is the time variable.
#' @param num.trees number of trees. Default is 500.
#' @param mtry number of variables to possibly split at in each node. Default is no. of variables^(3/4) ("^3/4") as recommended by (Ishwaran 2011). Also possible is "sqrt" and "0.5" to use the square root or half of the no. of variables.
#' @param type mode of prediction ("regression", "classification" or "survival"). Default is regression.
#' @param min.node.size minimal node size. Default is 1.
#' @param num.threads number of threads used for parallel execution. Default is number of CPUs available.
#' @param s predefined number of surrogate splits (it may happen that the actual number of surrogate splits differs in individual nodes). Default is 1 \% of no. of variables.
#' @param status status variable, only applicable to survival data. Use 1 for event and 0 for censoring.
#' @param save.ranger set TRUE if ranger object should be saved. Default is that ranger object is not saved (FALSE).
#' @param create.forest set FALSE if you want to analyze an existing forest. Default is TRUE.
#' @param forest the random forest that should be analyzed if create.forest is set to FALSE. (x and y still have to be given to obtain variable names)
#' @param save.memory Use memory saving (but slower) splitting mode. No effect for survival and GWAS data. Warning: This option slows down the tree growing, use only if you encounter memory problems. (This parameter is transfered to ranger)
#' @param case.weights Weights for sampling of training observations. Observations with larger weights will be selected with higher probability in the bootstrap (or subsampled) samples for the trees.
#'
#' @return list with the following components:
#' \itemize{
#' \item info: list with results from surrmindep function:
#' \itemize{
#' \item depth: mean surrogate minimal depth for each variable.
#' \item selected: variables has been selected (1) or not (0).
#' \item threshold: the threshold that is used for the selection.
#' }
#' \item var: vector of selected variables.
#'
#'\item s: list with the results of count.surrogate function:
#'\itemize{
#' \item s.a: total average number of surrogate variables.
#' \item s.l: average number of surrogate variables in the respective layers.
#'}
#' \item forest: a list containing:
#' #'\itemize{
#' \item trees: list of trees that was created by getTreeranger, addLayer, and addSurrogates functions and that was used for surrogate minimal depth variable importance.
#' \item allvariables: all variable names of the predictor variables that are present in x.
#' }
#'\item ranger: ranger object.
#'
#' }
#' @examples
#' # read data
#' data("SMD_example_data")
#'
#' \donttest{
#' # select variables (usually more trees are needed)
#' set.seed(42)
#' res = var.select.smd(x = SMD_example_data[,2:ncol(SMD_example_data)], y = SMD_example_data[,1],s = 10, num.trees = 10)
#' res$var
#' }
#'@references
##' \itemize{
##'   \item Seifert, S. et al. (2019) Surrogate minimal depth as an importance measure for variables in random forests. Bioinformatics, 35, 3663–3671. \url{https://academic.oup.com/bioinformatics/article/35/19/3663/5368013}
##'   \item Ishwaran, H. et al. (2011) Random survival forests for high-dimensional data. Stat Anal Data Min, 4, 115–132. \url{https://onlinelibrary.wiley.com/doi/abs/10.1002/sam.10103}
##'   \item Ishwaran, H. et al. (2010) High-Dimensional Variable Selection for Survival Data. J. Am. Stat. Assoc., 105, 205–217. \url{http://www.ccs.miami.edu/~hishwaran/papers/IKGML.JASA.2010.pdf}
##'   }
#' @export

var.select.smd = function(x = NULL, y = NULL, num.trees = 500, type = "regression", s = NULL, mtry = NULL, min.node.size = 1,
                          num.threads = NULL, status = NULL, save.ranger = FALSE, create.forest = TRUE, forest = NULL,
                          save.memory = FALSE, case.weights = NULL) {
  if(!is.data.frame(x)){
    stop("x has to be a data frame")
  }
  if (create.forest) {
  ## check data
  if (length(y) != nrow(x)) {
    stop("length of y and number of rows in x are different")
  }

  if (any(is.na(x))) {
    stop("missing values are not allowed")
  }

  variables = colnames(x)    # extract variables names
  nvar = length(variables)   # count number of variables

  ## set global parameters
  if (is.null(mtry)) {
    mtry = floor((nvar)^(3/4))
  }
  if (mtry == "sqrt") {
    mtry = floor(sqrt(nvar))
  }
  if (mtry == "0.5") {
    mtry = floor(0.5*nvar)
  }
  if (mtry == "^3/4") {
    mtry = floor((nvar)^(3/4))
  }



  if (is.null(s)) {
    s = ceiling(nvar*0.01)
  }

  if (s > (nvar - 1)) {
    s = nvar - 1
    warning("s was set to the maximum number that is reasonable (variables-1) ")
  }

  if (type == "classification") {
    y = as.factor(y)
    if (length(levels(y)) > 15) {
      stop("Too much classes defined, classification might be the wrong choice")
    }
  }
  if (type == "regression" && class(y) == "factor") {
    stop("use factor variable for y only for classification! ")
  }
  data = data.frame(y, x)
  if (type == "survival") {
    if (is.null(status)) {
      stop("a status variable named status has to be given for survival analysis")
    }
    data$status = status
    RF = ranger::ranger(data = data, dependent.variable.name = "y",num.trees = num.trees,mtry = mtry,min.node.size = min.node.size,
                        keep.inbag = TRUE, num.threads = num.threads, status.variable.name = "status", save.memory = save.memory,
                        case.weights = case.weights, respect.unordered.factors = "partition")
  }
  if (type == "classification" | type == "regression") {
  RF = ranger::ranger(data = data,dependent.variable.name = "y",num.trees = num.trees,mtry = mtry,min.node.size = min.node.size,
                      keep.inbag = TRUE, num.threads = num.threads, case.weights = case.weights, respect.unordered.factors = "partition")
  }
  trees = getTreeranger(RF = RF,num.trees = num.trees)
  trees.lay = addLayer(trees)
  rm(trees)
  ###AddSurrogates###
  trees.surr = addSurrogates(RF = RF,trees = trees.lay,s = s,Xdata = x, num.threads = num.threads)
  rm(trees.lay)
  # count surrogates
  s = count.surrogates(trees.surr)
  surrminimaldepth.s = surrmindep(forest = list(trees = trees.surr, allvariables = variables), s.l = s$s.l)
  }
  if (!create.forest) {
  if (is.null(forest)) {
    stop("set create.forest to TRUE or analyze an existing random forest specified by parameter forest")
  }
    allvariables = forest[["allvariables"]]
    trees = forest[["trees"]]
    s = count.surrogates(trees)
    surrminimaldepth.s = surrmindep(forest, s.l = s$s.l)
    trees.surr = forest[["trees"]]
    variables = forest[["variables"]]
  }
  if (save.ranger) {
  results = list(info = surrminimaldepth.s,var = names(surrminimaldepth.s$selected[surrminimaldepth.s$selected == 1]),s = s,
                 forest = list(trees = trees.surr, allvariables = variables), ranger = RF)
  }
  else {
    results = list(info = surrminimaldepth.s,var = names(surrminimaldepth.s$selected[surrminimaldepth.s$selected == 1]),s = s,
                   forest = list(trees = trees.surr, allvariables = variables))
  }
  return(results)
}
