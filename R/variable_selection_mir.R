#' Variable selection with mutual impurity reduction (MIR)
#'
#' This function executes MIR applying \link[ranger]{ranger} for random forests generation and actual impurity reduction and a modified version of \link[rpart]{rpart} to find surrogate variables.
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
#' @param p.t.sel p.value threshold for selection of important variables. Default is 0.01.
#' @param p.t.rel p.value threshold for selection of related variables. Default is 0.01.
#' @param num.permutations number of permutations to determine p-values. Default is 100. (the relations are determined once based on the permuted X data and the utilized AIR values are permuted again for each permutation )
#' @param status status variable, only applicable to survival data. Use 1 for event and 0 for censoring.
#' @param save.ranger set TRUE if ranger object should be saved. Default is that ranger object is not saved (FALSE).
#' @param save.memory Use memory saving (but slower) splitting mode. No effect for survival and GWAS data. Warning: This option slows down the tree growing, use only if you encounter memory problems. (This parameter is transfered to ranger)
#' @param select.var set False if only importance should be calculated and no variables should be selected.
#' @param select.rel set False if only relations should be calculated and no variables should be selected.
#' @param case.weights Weights for sampling of training observations. Observations with larger weights will be selected with higher probability in the bootstrap (or subsampled) samples for the trees.
#' @param method.rel Method  to  compute  p-values for selection of related variables with var.relations.corr. Use  "janitza"  for  the  method  by  Janitza  et  al. (2016) or "permutation" to utilize permuted variables.
#' @param method.sel Method  to  compute  p-values for selection of important variables. Use  "janitza"  for  the  method  by  Janitza  et  al. (2016) (can only be used when corrected variable relations are utilized) or "permutation" to utilize permuted variables.
#' @param corr.rel set FALSE if non-corrected variable relations should be used for calculation of MIR. In this case the method "janitza" should not be used for selection of important variables
#' @param t variable to calculate threshold for non-corrected relation analysis. Default is 5.
#' @param save.rel set FALSE if relation information should not bet saved (default is TRUE)
#'
#'
#' @return list with the following components:
#' \itemize{
#' \item info: list with results containing:
#' \itemize{
#' \item MIR: the calculated variable importance for each variable based on mutual impurity reduction.
#' \item pvalue: the obtained p-values for each variable.
#' \item selected: variables has been selected (1) or not (0).
#' \item relations: a list containing the results of variable relation analysis.
#' \item parameters: a list that contains the parameters s, type, mtry, p.t.sel, p.t.rel and method.sel that were used.
#' }
#' \item var: vector of selected variables.
#'
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
#' res = var.select.mir(x = SMD_example_data[,2:ncol(SMD_example_data)], y = SMD_example_data[,1],s = 10, num.trees = 10)
#' res$var
#' }
#'@references
##' \itemize{
##'   \item Nembrini, S. et al. (2018) The revival of the Gini importance? Bioinformatics, 34, 3711–3718. \url{https://academic.oup.com/bioinformatics/article/34/21/3711/4994791}
##'   \item Seifert, S. et al. (2019) Surrogate minimal depth as an importance measure for variables in random forests. Bioinformatics, 35, 3663–3671. \url{https://academic.oup.com/bioinformatics/article/35/19/3663/5368013}
##'   }
#' @export

var.select.mir = function(x = NULL, y = NULL, num.trees = 500, type = "regression", s = NULL, mtry = NULL, min.node.size = 1,
                          num.threads = NULL, status = NULL, save.ranger = FALSE,
                          save.memory = FALSE, num.permutations = 100, p.t.sel = 0.01, p.t.rel = 0.01, select.var = TRUE, select.rel = FALSE,
                          case.weights = NULL, corr.rel = TRUE, t = 5, method.rel = "permutation", method.sel = "janitza", save.rel = TRUE) {
  if(!is.data.frame(x)){
    stop("x has to be a data frame")
  }
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
                          num.threads = num.threads, status.variable.name = "status", save.memory = save.memory,
                          importance ="impurity_corrected", case.weights = case.weights, respect.unordered.factors = "partition")
      if (corr.rel) {
        rel = var.relations.mfi(x = x, y = y, num.trees = num.trees, type = type, s = s, mtry = mtry, min.node.size = min.node.size,
                                 num.threads = num.threads, status = status, case.weights = case.weights, variables = allvariables,
                                 candidates = allvariables, p.t = p.t.rel, method = method.rel,select.rel = select.rel)
      } else {
        rel = var.relations(x = x, y = y, num.trees = num.trees, type = type, s = s, mtry = mtry, min.node.size = min.node.size,
                            num.threads = num.threads, status = status, case.weights = case.weights, variables = allvariables,
                            candidates = allvariables, t = t, select.rel = select.rel)
      }
    }
    if (type == "classification" | type == "regression") {
      RF = ranger::ranger(data = data,dependent.variable.name = "y",num.trees = num.trees,mtry = mtry,min.node.size = min.node.size,
                          num.threads = num.threads, importance ="impurity_corrected", case.weights = case.weights, respect.unordered.factors = "partition")

      if (corr.rel) {
        rel = var.relations.mfi(x = x, y = y, num.trees = num.trees, type = type, s = s, mtry = mtry, min.node.size = min.node.size,
                                 num.threads = num.threads, case.weights = case.weights, variables = allvariables,
                                 candidates = allvariables, p.t = p.t.rel, method = method.rel,select.rel = select.rel)
      } else {
        rel = var.relations(x = x, y = y, num.trees = num.trees, type = type, s = s, mtry = mtry, min.node.size = min.node.size,
                            num.threads = num.threads, case.weights = case.weights, variables = allvariables,
                            candidates = allvariables, t = t,select.rel = select.rel)
      }

    }



adj.agree = rel$surr.res
diag(adj.agree) = 1

  mir = colSums(adj.agree * RF$variable.importance)

  if (select.var) {
    if (method.sel == "janitza") {
      if (corr.rel) {
      ## Mirrored VIMP (# This part is taken from ranger function)
      m1 = mir[mir< 0]
      m2 = mir[mir == 0]
      null.rel = c(m1, -m1, m2)

      pval <- 1 - ranger:::numSmaller(mir, null.rel) / length(null.rel)
      names(pval) = allvariables
      selected = as.numeric(pval <= p.t.sel)
      names(selected) = names(pval)

      if (length(m1) == 0) {
        stop("No negative importance values found for selection of important variables. Consider the 'permutation' approach.")
      }
      if (length(m1) < 100) {
        warning("Only few negative importance values found for selection of important variables, inaccurate p-values. Consider the 'permutation' approach.")
      }
      } else {
        stop("Janitza approach should only be conducted with corrected relations")
}
    }

    if (method.sel == "permutation") {

      if (corr.rel){
        adj.agree_perm = rel$surr.perm
      } else {
      x_perm = sapply(1:ncol(x),permute.variable,x=x)
      colnames(x_perm) = paste(allvariables,"_perm", sep = "")
      data_perm = data.frame(y, x_perm)
      allvariables_perm = colnames(x_perm)

      if (type == "survival") {
        if (is.null(status)) {
          stop("a status variables has to be given for survival analysis")
        }

        if (corr.rel) {
          rel_perm = var.relations.mfi(x = data.frame(x_perm), y = y, num.trees = num.trees, type = type, s = s, mtry = mtry, min.node.size = min.node.size,
                                   num.threads = num.threads, status = status, case.weights = case.weights, variables = allvariables,
                                   candidates = allvariables, p.t = p.t.rel, method = method.rel, select.rel = select.rel)
        } else {
          rel_perm = var.relations(x = data.frame(x_perm), y = y, num.trees = num.trees, type = type, s = s, mtry = mtry, min.node.size = min.node.size,
                              num.threads = num.threads, status = status, case.weights = case.weights, variables = allvariables,
                              candidates = allvariables, t = t, select.rel = select.rel)
        }

      }
      if (type == "classification" | type == "regression") {

        if (corr.rel) {
          rel_perm = var.relations.mfi(x = data.frame(x_perm), y = y, num.trees = num.trees, type = type, s = s, mtry = mtry, min.node.size = min.node.size,
                                   num.threads = num.threads, case.weights = case.weights, variables = allvariables_perm,
                                   candidates = allvariables_perm, p.t = p.t.rel, method = method.rel,select.rel = select.rel)
        } else {
          rel_perm = var.relations(x = data.frame(x_perm), y = y, num.trees = num.trees, type = type, s = s, mtry = mtry, min.node.size = min.node.size,
                              num.threads = num.threads, case.weights = case.weights, variables = allvariables_perm,
                              candidates = allvariables_perm, t = t,select.rel = select.rel)
        }
        }

      adj.agree_perm = rel_perm$surr.res
      }

      diag(adj.agree_perm) = 0

      null.rel = unlist(lapply(1:num.permutations,calculate.mir.perm,
                               adj.agree_perm = adj.agree_perm,
                               air = RF$variable.importance,
                               allvariables = allvariables))



      pval <- 1 - ranger:::numSmaller(mir, null.rel) / length(null.rel)
      names(pval) = allvariables
      selected = as.numeric(pval <= p.t.sel)
      names(selected) = names(pval)

      }
if(save.rel) {
  info = list(MIR = mir,
              pvalue = pval,
              selected = selected,
              relations = rel,
              AIR = RF$variable.importance,
              parameters = list(s = s, type = type, mtry = mtry, p.t.sel = p.t.sel, p.t.rel = p.t.rel, method.sel = method.sel))
} else {
  info = list(MIR = mir,
              pvalue = pval,
              selected = selected,
              AIR = RF$variable.importance,
              parameters = list(s = s, type = type, mtry = mtry, p.t.sel = p.t.sel, p.t.rel = p.t.rel, method.sel = method.sel))
}
  } else {
    if(save.rel) {
    info = list(MIR = mir,
                relations = rel,
                AIR = RF$variable.importance,
                parameters = list(s = s, type = type, mtry = mtry))
    } else {
      info = list(MIR = mir,
                  AIR = RF$variable.importance,
                  parameters = list(s = s, type = type, mtry = mtry))
    }
  }

  if (save.ranger) {
    results = list(info = info,
                   var = names(info$selected[info$selected == 1]),
                   ranger = RF)
  } else {
    if (select.var) {
    results = list(info = info,
                   var = names(info$selected[info$selected == 1]))
  } else {
    results = list(info = info)
  }
  }
  return(results)
}


#'
#' This is an internal function
#'
#' @keywords internal
calculate.mir.perm = function(r=1, adj.agree_perm, air, allvariables) {
mir.perm = colSums(adj.agree_perm * sample(air,length(air)))
return(mir.perm)
}


