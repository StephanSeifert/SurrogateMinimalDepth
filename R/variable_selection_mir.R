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
#' @param p.t.sel p.value threshold for selection of important variables. Default is 0.01
#' @param p.t.rel p.value threshold for selection of related variables. Default is 0.01
#' @param min.var.p minimum number of permuted variables used to determine p-value for variable selection of important and related variables. Default is 200.
#' @param status status variable, only applicable to survival data. Use 1 for event and 0 for censoring.
#' @param save.ranger set TRUE if ranger object should be saved. Default is that ranger object is not saved (FALSE).
#' @param save.memory Use memory saving (but slower) splitting mode. No effect for survival and GWAS data. Warning: This option slows down the tree growing, use only if you encounter memory problems. (This parameter is transfered to ranger)
#' @param select.var set False if only importance and relations should be calculated and no variables should be selected
#' @param case.weights Weights for sampling of training observations. Observations with larger weights will be selected with higher probability in the bootstrap (or subsampled) samples for the trees.

#' @return list with the following components:
#' \itemize{
#' \item info: list with results containing:
#' \itemize{
#' \item MIR: the calculated variable importance for each variable based on mutual impurity reduction.
#' \item pvalue: the obtained p-values for each variable.
#' \item selected: variables has been selected (1) or not (0).
#' \item relations: a list containing the relation analysis results:
#'  \itemize{
#'  \item corr.maa: a matrix with corrected mean adjusted agreement values. For each variable in the rows the relation parameter of each variable can be found in the columns.
#'  \item p.rel: a list with the obtained p-values for the relation analysis of each variable
#'  \item var.rel: a list with vectors of related variables for each variable
#'  }
#' \item parameters: a list that contains the parameters s, type, mtry, p.t.sel, p.t.rel that were used.
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
                          num.threads = NULL, status = NULL, save.ranger = FALSE,
                          save.memory = FALSE, min.var.p = 200, p.t.sel = 0.01, p.t.rel = 0.01, select.var = TRUE,
                          case.weights = NULL) {

    ## check data
    if (length(y) != nrow(x)) {
      stop("length of y and number of rows in x are different")
    }

    if (any(is.na(x))) {
      stop("missing values are not allowed")
    }

    variables = colnames(x)# extract variables names
    f = ceiling(min.var.p / (ncol(x))) # f determines how often the variables are permuted
    nvar = length(variables)  #count variables
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
    # create shadow variables to correct the relation and to determine the p-value for selection
    x_perm = sapply(1:ncol(x),permute.variable,x=x)
    colnames(x_perm) = paste(variables,"_perm", sep = "")
  if (select.var) {
  if ( ncol(x) < min.var.p) {
    message("More permuted variables than original variables are needed so they are permuted multiple times.")
    x_perm2 = matrix(rep(sapply(1:ncol(x),permute.variable,x=x), (f-1)),nrow = nrow(x), ncol= ncol(x) * (f-1))
    colnames(x_perm2) = rep(paste(variables,"_perm", sep = ""),(f-1))
    x_perm = cbind(x_perm,x_perm2)
  }
  }
    data = data.frame(y, x)
    data_perm = data.frame(y, x_perm)
    if (type == "survival") {
      if (is.null(status)) {
        stop("a status variables has to be given for survival analysis")
      }
      data$status = status
      RF = ranger::ranger(data = data,dependent.variable.name = "y",num.trees = ntree,mtry = mtry,min.node.size = min.node.size,
                          keep.inbag = TRUE, num.threads = num.threads, dependent.variable.name = "status", save.memory = save.memory,
                          importance ="impurity_corrected", case.weights = case.weights)
      RF_perm = ranger::ranger(data = data_perm,dependent.variable.name = "y",num.trees = ntree,mtry = mtry,min.node.size = min.node.size,
                          keep.inbag = TRUE, num.threads = num.threads, dependent.variable.name = "status", save.memory = save.memory,
                          importance ="impurity_corrected", case.weights = case.weights)
    }
    if (type == "classification" | type == "regression") {
      RF = ranger::ranger(data = data,dependent.variable.name = "y",num.trees = ntree,mtry = mtry,min.node.size = min.node.size,
                          keep.inbag = TRUE, num.threads = num.threads, importance ="impurity_corrected", case.weights = case.weights)
      RF_perm = ranger::ranger(data = data_perm,dependent.variable.name = "y",num.trees = ntree,mtry = mtry,min.node.size = min.node.size,
                          keep.inbag = TRUE, num.threads = num.threads, importance ="impurity_corrected", case.weights = case.weights)
    }
    trees = getTreeranger(RF = RF,ntree = ntree)
    trees.lay = addLayer(trees)
    rm(trees)
    ###AddSurrogates###
    trees.surr = addSurrogates(RF = RF,trees = trees.lay,s = s,Xdata = data[,-1], num.threads = num.threads)
    rm(trees.lay)
    forest = list(trees = trees.surr, allvariables = colnames(data[,-1]))

    # do the same for the permutation forrest
    trees_perm = getTreeranger(RF = RF_perm,ntree = ntree)
    trees.lay_perm = addLayer(trees_perm)
    rm(trees_perm)
    ###AddSurrogates###
    trees.surr_perm = addSurrogates(RF = RF_perm,trees = trees.lay_perm,s = s,Xdata = data_perm[,-1], num.threads = num.threads)
    rm(trees.lay_perm)
    forest_perm = list(trees = trees.surr_perm, allvariables = colnames(data_perm[,-1]))


  # determine variable relations
  rel = var.relations(forest = forest,variables = colnames(data[,-1]), candidates = colnames(data[,-1]))
  rel_perm = var.relations(forest = forest_perm,variables = colnames(data_perm[,-1]), candidates = colnames(data_perm[,-1]))

  adj.agree = rel$surr.res
  adj.agree.perm = rel_perm$surr.res
  adj.agree.corr = adj.agree - adj.agree.perm


  if (select.var) {
  rel.p = lapply(1:length(variables),p.relation,
                  adj.agree = adj.agree,
                 adj.agree.perm = adj.agree.perm,
                  adj.agree.corr = adj.agree.corr,
                  variables = variables)

  sel.rel = lapply(1:length(variables),select.related,
                   rel.p,
                   p.t.rel)

  names(rel.p) = names(sel.rel) = variables

  relations = list(corr.maa = adj.agree.corr,
                   p.rel = rel.p,
                   var.rel = sel.rel)
  } else {
    relations = list(corr.maa = adj.agree.corr)
  }
  adj.agree.corr[which(is.na(adj.agree.corr))] = 1
  mir = rowSums(t(t(adj.agree.corr) * RF$variable.importance[1:length(variables)]))

  if (select.var) {
  mir.perm = unlist(lapply(1:f,calculate.mir.perm,
                    adj.agree.perm = adj.agree.perm,
                    air = RF_perm$variable.importance,
                    variables = variables))

  # compute p-values using numSmaller function from ranger
  pval <- 1 - ranger:::numSmaller(mir[1:length(variables)], mir.perm) / length(mir.perm)
  names(pval) = variables
  selected = as.numeric(pval < p.t.sel)
  names(selected) = names(pval)

  info = list(MIR = mir,
              pvalue = pval,
              selected = selected,
              relations = relations,
              parameters = list(s = s, type = type, mtry = mtry, p.t.sel = p.t.sel, p.t.rel = p.t.rel))
  } else {
    info = list(MIR = mir,
                relations = relations,
                parameters = list(s = s, type = type, mtry = mtry))
}

  if (save.ranger) {
    results = list(info = info,
                   var = names(info$selected[info$selected == 1]),
                   forest = forest,
                   ranger = RF)
  }
  else {
    results = list(info = info,
                   var = names(info$selected[info$selected == 1]),
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

#' calculate.mir.perm
#'
#' This is an internal function
#'
#' @keywords internal
calculate.mir.perm = function(r=1, adj.agree.perm, air, variables) {
mir.perm = rowSums(t(t(adj.agree.perm) * air[((r-1) * length(variables) + 1):(r * length(variables))]),na.rm = TRUE)
return(mir.perm)
}

#' calculate.mir.perm
#'
#' This is an internal function
#'
#' @keywords internal
p.relation = function(l = 1,
                      adj.agree,
                      adj.agree.perm,
                      adj.agree.corr,
                      variables) {
 relations = adj.agree.corr[l,]
 null.rel = adj.agree.perm[l,]
 pval <- 1 - ranger:::numSmaller(relations, null.rel) / length(null.rel)
 names(pval) = variables
 pval[l] = NA
 return(pval)
}


#' calculate.mir.perm
#'
#' This is an internal function
#'
#' @keywords internal
select.related = function(m=1,
                          rel.p,
                          p.t) {
  rel.var = rel.p[[m]]
  names(which(rel.var < p.t))
}
