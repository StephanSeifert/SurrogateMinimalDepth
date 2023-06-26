#' Investigate variable relations of a specific variable with mean adjusted agreement
#'
#'This function uses the mean adjusted agreement to select variables that are related to a defined variable using a threshold T.
#'The parameter t is used to calculate T: t=1 means that every variable with higher probability than "by chance" is identified
#'as "important". t=2 means the probability has to be twice, etc.
#'Based on the threshold a vector is created containing the related variables.
#'
#' @param variables variable names (string) for which related variables should be searched for (has to be contained in allvariables)
#' @param candidates vector of variable names (strings) that are candidates to be related to the variables (has to be contained in allvariables)
#' @param t variable to calculate threshold. Default is 5.
#' @param select.rel set False if only relations should be calculated and no related variables should be selected.
#' @param num.threads number of threads used for determination of relations. Default is number of CPUs available.
#' @inheritParams var.select.smd
#'
#' @return a list containing:
#' \itemize{
#' \item variables: the variables to which relations are investigated.
#' \item surr.res: a matrix with mean adjusted agreement values with variables in rows and candidates in columns.
#' \item threshold: the threshold used to select related variables.
#' \item var: a list with one vector for each variable containing related variables.
#' \item ranger: ranger object.
#' }
#' @examples
#' # read data
#' data("SMD_example_data")
#' x = SMD_example_data[,2:ncol(SMD_example_data)]
#' y = SMD_example_data[,1]
#' \donttest{
#' # calculate variable relations
#' set.seed(42)
#' res = var.relations(x = x, y = y, s = 10, num.trees = 100, variables = c("X1","X7"), candidates = colnames(x)[1:100], t = 5)
#' res$var
#' }
#'
#' @export

var.relations = function(x = NULL, y = NULL, num.trees = 500, type = "regression", s = NULL, mtry = NULL, min.node.size = 1,
                         num.threads = NULL, status = NULL, save.ranger = FALSE, create.forest = TRUE, forest = NULL,
                         save.memory = FALSE, case.weights = NULL,
                         variables, candidates, t = 5, select.rel = TRUE) {
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

    allvariables = colnames(x)# extract variables names
    nvar = length(allvariables)   # count number of variables
    ## set global parameters
    if (is.null(mtry)) {
      mtry = floor((nvar)^(3/4))
    }
    if (mtry == "sqrt") {
      mtry = floor(sqrt(nvar))
    }
    if (mtry == "0.5") {
      mtry = floor(0.5*(nvar))
    }
    if (mtry == "^3/4") {
      mtry = floor((nvar)^(3/4))
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

    data = data.frame(y, x)

    if (type == "survival") {
      if (is.null(status)) {
        stop("a status variable named status has to be given for survival analysis")
      }
      data$status = status
      RF = ranger::ranger(data = data,dependent.variable.name = "y",num.trees = num.trees,mtry = mtry,min.node.size = min.node.size,
                          keep.inbag = TRUE, num.threads = num.threads, status.variable.name = "status", save.memory = save.memory,
                          case.weights = case.weights)
    }
    if (type == "classification" | type == "regression") {
      RF = ranger::ranger(data = data,dependent.variable.name = "y",num.trees = num.trees,mtry = mtry,min.node.size = min.node.size,
                          keep.inbag = TRUE, num.threads = num.threads, case.weights = case.weights)

    }

    trees = getTreeranger(RF = RF,num.trees = num.trees)
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
  }
  trees = forest[["trees"]]
  allvariables = forest[["allvariables"]]

  if (all(candidates %in% allvariables)) {
  if (all(variables %in% allvariables)) {
  # count surrogates
  s = count.surrogates(trees)
  results.meanAdjAgree = meanAdjAgree(trees, variables, allvariables, candidates, t = t, s$s.a, select.var = select.rel,num.threads = num.threads)
  } else {
    stop("allvariables do not contain the chosen variables")
  }
  } else {
    stop("allvariables do not contain the candidate variables")
  }
  if(select.rel) {
  surr.var = results.meanAdjAgree$surr.var
  varlist = list()
  for (i in 1:nrow(surr.var)) {
    surr.var.var = surr.var[i,]
    if (anyNA(surr.var.var)) {
    surr.var.var = surr.var.var[-which(is.na(surr.var.var))]
    }
    var = names(surr.var.var[surr.var.var == 1])
    name = variables[i]
    varlist[[name]] = var
  }

  var = names(surr.var[surr.var == 1])
  if (save.ranger) {
 return(list(variables = results.meanAdjAgree$variables, surr.res = results.meanAdjAgree$surr.res,
              threshold = results.meanAdjAgree$threshold, var = varlist, ranger = RF))
  } else {
    return(list(variables = results.meanAdjAgree$variables, surr.res = results.meanAdjAgree$surr.res,
                threshold = results.meanAdjAgree$threshold, var = varlist))
  }

  } else {
    if (save.ranger) {
    return(list(variables = results.meanAdjAgree$variables, surr.res = results.meanAdjAgree$surr.res,ranger = RF))

  } else {
    return(list(variables = results.meanAdjAgree$variables, surr.res = results.meanAdjAgree$surr.res))
  }
  }
}
