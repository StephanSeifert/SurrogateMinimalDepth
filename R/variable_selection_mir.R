#' Variable selection with mutual impurity reduction (MIR)
#'
#' This function executes MIR applying \link[ranger]{ranger} for random forests generation and actual impurity reduction and a modified version of \link[rpart]{rpart} to find surrogate variables.
#'
#' @param x matrix or data.frame of predictor variables with variables in
#'   columns and samples in rows (Note: missing values are not allowed)
#' @param y vector with values of phenotype variable (Note: will be converted to factor if
#'   classification mode is used). For survival forests this is the time variable.
#' @param ntree number of trees. Default is 500.
#' @param mtry number of variables to possibly split at in each node. Default is no. of variables^(3/4) ("^3/4") as recommended by (Ishwaran 2011). Also possible is "sqrt" and "0.5" to use the square root or half of the no. of variables.
#' @param type mode of prediction ("regression" or "classification"). Default is regression.
#' @param min.node.size minimal node size. Default is 1.
#' @param num.threads number of threads used for parallel execution. Default is number of CPUs available.
#' @param s predefined number of surrogate splits (it may happen that the actual number of surrogate splits differs in individual nodes). Default is 1 \% of no. of variables.
#' @param p.t p.value threshold for variable selection. Default is 0.01
#' @param num.p number of permuted variables used to determine p-value for variable selection. Default is 500.
#' @param status status variable, only applicable to survival data. Use 1 for event and 0 for censoring.
#' @param save.ranger set TRUE if ranger object should be saved. Default is that ranger object is not saved (FALSE).
#' @param create.forest set FALSE if you want to analyze an existing forest. Default is TRUE.
#' @param forest the random forest that should be analyzed if create.forest is set to FALSE. (x and y still have to be given to obtain variable names)
#' @param save.memory Use memory saving (but slower) splitting mode. No effect for survival and GWAS data. Warning: This option slows down the tree growing, use only if you encounter memory problems. (This parameter is transfered to ranger)

#' @return list with the following components:
#' \itemize{
#' \item info: list with results containing:
#' \itemize{
#' \item MIR: the calculated variable importance for each variable based on mutual impurity reduction.
#' \item pvalue: the obtained p-values for each variable.
#' \item selected: variables has been selected (1) or not (0).
#' \item parameters: a list that contains the parameters s, type, mtry and p.t that were used.
#' }
#' \item var: vector of selected variables.
#'
#' \item forest: a list containing:
#' \itemize{
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
#' res = var.select.mir(x = SMD_example_data[,2:ncol(SMD_example_data)], y = SMD_example_data[,1],s = 10, ntree = 10)
#' res$var
#' }
#'@references
##' \itemize{
##'   \item Nembrini, S. et al. (2018) The revival of the Gini importance? Bioinformatics, 34, 3711–3718. \url{https://academic.oup.com/bioinformatics/article/34/21/3711/4994791}
##'   \item Seifert, S. et al. (2019) Surrogate minimal depth as an importance measure for variables in random forests. Bioinformatics, 35, 3663–3671. \url{https://academic.oup.com/bioinformatics/article/35/19/3663/5368013}
##'   }
#' @export

var.select.mir = function(x = NULL, y = NULL, ntree = 500, type = "regression", s = NULL, mtry = NULL, min.node.size = 1,
                          num.threads = NULL, status = NULL, save.ranger = FALSE, create.forest = TRUE, forest = NULL,
                          save.memory = FALSE, num.p = 500, p.t = 0.01) {
  if (create.forest) {
    ## check data
    if (length(y) != nrow(x)) {
      stop("length of y and number of rows in x are different")
    }

    if (any(is.na(x))) {
      stop("missing values are not allowed")
    }

    variables = colnames(x)# extract variables names
    nvar = length(variables)   # count number of variables

    ## set global parameters
    if (is.null(mtry)) {
      mtry = floor((nvar + num.p)^(3/4))
    }
    if (mtry == "sqrt") {
      mtry = floor(sqrt(nvar + num.p))
    }
    if (mtry == "0.5") {
      mtry = floor(0.5*(nvar + num.p))
    }
    if (mtry == "^3/4") {
      mtry = floor((nvar + num.p)^(3/4))
    }


    if (is.null(s)) {
      s = ceiling(nvar*0.01)
    }

    if (s > (nvar - 2)) {
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
    # create shadow variables to simulate unimportance variables for the determination of the p-value

  if (num.p > ncol(x)) {
    message("More permuted variables than original variables are needed. They are sampled with replacement")
    var.perm = sample(c(1:ncol(x)),num.p, replace = TRUE)
    x_perm = sapply(var.perm,permute.variable,x=x)
    colnames(x_perm) = paste(variables[var.perm],"_perm", sep = "")
  } else {
    var.perm = sample(c(1:ncol(x)),num.p)
    x_perm = sapply(var.perm,permute.variable,x=x)
    colnames(x_perm) = paste(variables[var.perm],"_perm", sep = "")
  }
    data = data.frame(y, x, x_perm)
    if (type == "survival") {
      if (is.null(status)) {
        stop("a status variables has to be given for survival analysis")
      }
      data$status = status
      RF = ranger::ranger(data = data,dependent.variable.name = "y",num.trees = ntree,mtry = mtry,min.node.size = min.node.size,
                          keep.inbag = TRUE, num.threads = num.threads, dependent.variable.name = "status", save.memory = save.memory,
                          importance ="impurity_corrected")
    }
    if (type == "classification" | type == "regression") {
      RF = ranger::ranger(data = data,dependent.variable.name = "y",num.trees = ntree,mtry = mtry,min.node.size = min.node.size,
                          keep.inbag = TRUE, num.threads = num.threads, importance ="impurity_corrected")
    }
    trees = getTreeranger(RF = RF,ntree = ntree)
    trees.lay = addLayer(trees)
    rm(trees)
    ###AddSurrogates###
    trees.surr = addSurrogates(RF = RF,trees = trees.lay,s = s,Xdata = data[,-1], num.threads = num.threads)
    rm(trees.lay)
    forest = list(trees = trees.surr, allvariables = colnames(data[,-1]))

  }
  if (!create.forest) {
    if (is.null(forest)) {
      stop("set create.forest to TRUE or analyze an existing random forest specified by parameter forest")
    }
    trees.surr = forest[["trees"]]
    variables = forest[["variables"]]
  }

  # determine variable relations
  rel = var.relations(forest = forest,variables = colnames(data[,-1]), candidates = colnames(data[,-1]))
  adj.agree = rel$surr.res
  adj.agree[which(is.na(adj.agree))] = 1
  mir = rowSums(t(t(adj.agree) * RF$variable.importance))
  vimp = mir[(length(variables) + 1):length(mir)]

  # compute p-values using numSmaller function from ranger
  pval <- 1 - ranger:::numSmaller(mir[1:length(variables)], vimp) / length(vimp)
  names(pval) = variables
  selected = as.numeric(pval < p.t)
  names(selected) = names(pval)

  info = list(MIR = mir[1:length(variables)],
              pvalue = pval,
              selected =selected,
              parameters = list(s = s, type = type, mtry = mtry, p.t = p.t))


  if (save.ranger) {
    results = list(info = info,var = names(info$selected[info$selected == 1]),
                   forest = forest, ranger = RF)
  }
  else {
    results = list(info = info,var = names(info$selected[info$selected == 1]),
                   forest = forest)
  }
  return(results)
}

#' permute.variable
#'
#' This is an internal function
#'
#' @keywords internal
permute.variable=function(i=1,x){
  var.perm = sample(x[,i],nrow(x))
  return(var.perm)
}
