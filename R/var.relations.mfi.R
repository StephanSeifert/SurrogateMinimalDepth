#' Investigate variable relations of a specific variable with mutual forest impact (corrected mean adjusted agreement).
#'
#'This function corrects the mean adjusted agreement by a permutation approach and generates the relation parameter mutual forest impact. Subsequently p-values are determined and related variables are selected.
#'
#' @param variables variable names (string) for which related variables should be searched for (has to be contained in allvariables)
#' @param candidates vector of variable names (strings) that are candidates to be related to the variables (has to be contained in allvariables)
#' @param p.t p.value threshold for selection of related variables. Default is 0.01.
#' @param select.rel set False if only relations should be calculated and no related variables should be selected.
#' @param method Method  to  compute  p-values.   Use  "janitza"  for  the  method  by  Janitza  et  al. (2016) or "permutation" to utilize permuted relations.
#' @param num.threads number of threads used for determination of relations. Default is number of CPUs available.
#' @inheritParams var.select.smd
#'
#' @return a list containing:
#' \itemize{
#' \item variables: the variables to which relations are investigated.
#' \item surr.res: a matrix with the mutual forest impact values with variables in rows and candidates in columns.
#' \item surr.perm: a matrix with the mutual forest impact values of the permuted variables with variables in rows and candidates in columns.
#' \item p.rel: a list with the obtained p-values for the relation analysis of each variable.
#' \item var.rel: a list with vectors of related variables for each variable.
#' \item ranger: ranger objects.
#' \item method: Method  to  compute  p-values: "janitza" or "permutation".
#' \item p.t: p.value threshold for selection of related variables
#'
#'
#' }
#' @examples
#' # read data
#' data("SMD_example_data")
#' x = SMD_example_data[,2:ncol(SMD_example_data)]
#' y = SMD_example_data[,1]
#' \donttest{
#' # calculate variable relations
#' set.seed(42)
#' res = var.relations.mfi(x = x, y = y, s = 10, num.trees = 100, variables = c("X1","X7"), candidates = colnames(x)[1:100])
#' res$var.rel[[1]]
#' }
#'
#' @export

var.relations.mfi = function(x = NULL, y = NULL, num.trees = 500, type = "regression", s = NULL, mtry = NULL, min.node.size = 1,
                         num.threads = NULL, status = NULL, save.ranger = FALSE, create.forest = TRUE, forest = NULL,
                         save.memory = FALSE, case.weights = NULL,
                         variables, candidates, p.t = 0.01, select.rel = TRUE, method = "janitza") {
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

    # create shadow variables to correct the relation
    x_perm = data.frame(lapply(1:ncol(x),permute.variable,x=x))
    colnames(x_perm) = paste(allvariables,"_perm", sep = "")

    data = data.frame(y, x)
    data_perm = data.frame(y, x_perm)

    if (type == "survival") {
      if (is.null(status)) {
        stop("a status variable named status has to be given for survival analysis")
      }
      data$status = status
      RF = ranger::ranger(data = data,dependent.variable.name = "y",num.trees = num.trees,mtry = mtry, min.node.size = min.node.size,
                          keep.inbag = TRUE, num.threads = num.threads, status.variable.name = "status", save.memory = save.memory,
                          case.weights = case.weights, respect.unordered.factors = "partition")
      data_perm$status = status
      RF_perm = ranger::ranger(data = data_perm,dependent.variable.name = "y",num.trees = num.trees,mtry = mtry,min.node.size = min.node.size,
                          keep.inbag = TRUE, num.threads = num.threads, status.variable.name = "status", save.memory = save.memory,
                          case.weights = case.weights, respect.unordered.factors = "partition")
    }
    if (type == "classification" | type == "regression") {
      RF = ranger::ranger(data = data,dependent.variable.name = "y",num.trees = num.trees,mtry = mtry,min.node.size = min.node.size,
                          keep.inbag = TRUE, num.threads = num.threads, case.weights = case.weights, respect.unordered.factors = "partition")

      RF_perm = ranger::ranger(data = data_perm,dependent.variable.name = "y",num.trees = num.trees,mtry = mtry,min.node.size = min.node.size,
                          keep.inbag = TRUE, num.threads = num.threads, case.weights = case.weights, respect.unordered.factors = "partition")

    }
    trees = getTreeranger(RF = RF,num.trees = num.trees)
    trees.lay = addLayer(trees)
    rm(trees)
    ###AddSurrogates###
    trees.surr = addSurrogates(RF = RF,trees = trees.lay,s = s,Xdata = data[,-1], num.threads = num.threads)
    rm(trees.lay)
    forest = list(trees = trees.surr, allvariables = colnames(data[,-1]))

    # do the same for the permutation forrest
    trees_perm = getTreeranger(RF = RF_perm,num.trees = num.trees)
    trees.lay_perm = addLayer(trees_perm)
    rm(trees_perm)
    ###AddSurrogates###
    trees.surr_perm = addSurrogates(RF = RF_perm,trees = trees.lay_perm,s = s,Xdata = data_perm[,-1], num.threads = num.threads)
    rm(trees.lay_perm)
    forest_perm = list(trees = trees.surr_perm, allvariables = colnames(data_perm[,-1]))
  }

  if (!create.forest) {
    if (is.null(forest)) {
      stop("set create.forest to TRUE or analyze an existing random forest specified by parameter forest")
    }
  }

  if (all(candidates %in% allvariables)) {
  if (all(variables %in% allvariables)) {
  # count surrogates
  s = count.surrogates(forest$trees)
  rel = meanAdjAgree(forest$trees, variables = allvariables, allvariables = allvariables, candidates = allvariables,
                     t = t, s$s.a, select.var = FALSE, num.threads = num.threads)

  allvariables_perm = colnames(x_perm)

  rel_perm = meanAdjAgree(forest_perm$trees, variables = allvariables_perm , allvariables = allvariables_perm , candidates = allvariables_perm,
                          t = t, s$s.a, select.var = FALSE, num.threads = num.threads)

  } else {
    stop("allvariables do not contain the chosen variables")
  }
  } else {
    stop("allvariables do not contain the candidate variables")
  }

  adj.agree = rel$surr.res
  adj.agree.perm = rel_perm$surr.res
  diag(adj.agree) = diag(adj.agree.perm) = 1

  if(anyNA(adj.agree)) {
    no.na = length(which(rowSums(is.na(adj.agree)) != 0 ))
    warning(paste0("Relations for ", no.na, " original variables were not calculated because they were never used as a primary split.
            Affected relations are set to 0. "))
    adj.agree[which(is.na(adj.agree))] = 0
  }

  if(anyNA(adj.agree.perm)) {
    no.na = length(which(rowSums(is.na(adj.agree.perm)) != 0 ))
    warning(paste0("Relations for ", no.na, " permuted variables were not calculated because they were not used as a primary split.
            Affected relations are set to 0. "))
    adj.agree.perm[which(is.na(adj.agree.perm))] = 0
  }
  adj.agree.corr = adj.agree - adj.agree.perm[1:nvar,1:nvar]
  diag(adj.agree.corr) =  diag(adj.agree) = diag(adj.agree.perm) =  NA

  adj.agree.corr.var = adj.agree.corr[variables,candidates]

  if (select.rel) {

    if (method == "janitza") {

      adj.agree.1 = adj.agree.corr
      diag(adj.agree.1) = 1
      ## Mirrored VIMP (# This part is taken from ranger function)
      m1 = adj.agree.1[adj.agree.1 < 0]
      m2 = adj.agree.1[adj.agree.1 == 0]
      null.rel = c(m1, -m1, m2)

      if (length(m1) == 0) {
        stop("No negative importance values found for selection of related variables. Consider the 'permutation' approach.")
      }
      if (length(m1) < 100) {
        warning("Only few negative importance values found for selection of related variables, inaccurate p-values. Consider the 'permutation' approach.")
      }

      rel.p = lapply(1:length(variables),p.relation,
                     null.rel = null.rel,
                     adj.agree.corr = adj.agree.corr.var,
                     candidates = candidates,
                     variables = variables)
      sel.rel = lapply(1:length(variables),select.related,
                       rel.p,
                       p.t)

      names(rel.p) = names(sel.rel) = variables

    }

    if (method == "permutation") {

    null.rel.plus = as.vector(adj.agree.perm)
    null.rel.plus = null.rel.plus[!is.na(null.rel.plus)]

    m1 = null.rel.plus[null.rel.plus > 0]
    m2 = null.rel.plus[null.rel.plus == 0]
    null.rel = c(m1, -m1, m2)


    if (length(null.rel) < 100) {
      warning("Only few null relations used. P-values could be inaccurate.")
    }

    rel.p = lapply(1:length(variables),p.relation,
                   null.rel = null.rel,
                   adj.agree.corr = adj.agree.corr.var,
                   candidates = candidates,
                   variables = variables)

    sel.rel = lapply(1:length(variables),select.related,
                     rel.p,
                     p.t)

    names(rel.p) = names(sel.rel) = variables
  }

  if (save.ranger) {
 return(list(variables = variables, surr.res = adj.agree.corr.var, surr.perm = adj.agree.perm,
              p.rel = rel.p, var.rel = sel.rel, ranger = list(RF = RF,RF_perm = RF_perm), method = method, p.t = p.t))
  } else {
    return(list(variables = variables, surr.res = adj.agree.corr.var, surr.perm = adj.agree.perm, p.rel = rel.p, var.rel = sel.rel, method = method, p.t = p.t))
  }
  } else {
    if (save.ranger) {
    return(list(variables = variables, surr.res = adj.agree.corr.var, surr.perm = adj.agree.perm, ranger = list(RF = RF,RF_perm = RF_perm)))
  } else {
    return(list(variables = variables, surr.res = adj.agree.corr.var, surr.perm = adj.agree.perm))
  }
  }
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

#' p.relation
#'
#' This is an internal function
#'
#' @keywords internal
p.relation = function(l = 1,
                      null.rel,
                      adj.agree.corr,
                      candidates,
                      variables) {
  relations = adj.agree.corr[l,]
  pval <- 1 - ranger:::numSmaller(relations, null.rel) / length(null.rel)
  names(pval) = candidates
  pval[variables[l]] = NA
  return(pval)
}

#' select.related
#'
#' This is an internal function
#'
#' @keywords internal
select.related = function(m=1,
                          rel.p,
                          p.t) {
  rel.var = rel.p[[m]]
  names(which(rel.var <= p.t))
}
