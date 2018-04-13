#' Investigate variable relations of a specific variable with mean adjusted agreement
#'
#'This function uses the mean adjusted agreement to select variables that are related to a defined variable using a threshold T.
#'The parameter t is used to calculate T: t=1 means that every variable with higher probability than "by chance" is identified
#'as "important". t=2 means the probability has to be twice, etc.
#'Based on the threshold a vector is created containing the related variables. In order to get all variable relations
#'
#' @param trees list of trees that was generated by getTreeranger function and layers, surrogate variables, and adjusted agreement values were added by addLayer and getSurrogates functions
#' @param variables variable names (string) for which related variables should be searched for (has to be contained in allvariables)
#' @param allvariables vector of all variable names in the original data set (strings).
#' @param candidates vector of variable names (strings) that are candidates to be related to the variables (has to be contained in allvariables)
#' @param t variable to calculate threshold. Default is 5.

#'
#' @return a list containing:
#' \itemize{
#' \item variables: the variables to which relations are investigated
#' \item surr.res: a matrix with mean adjusted agreement values with variables in rows and candidates in columns
#' \item threshold: the threshold used to select related variables
#' \item var: a list with one vector for each variable containing related variables
#' }
#' @examples
#' # read data
#' data("SMD_example_data")
#'
#' \donttest{
#' ###### use result of SMD variable importance
#' # select variables with smd variable importance (usually more trees are needed)
#' set.seed(42)
#' res = var.select.smd(x=SMD_example_data[,2:ncol(SMD_example_data)],y=SMD_example_data[,1],s=10,ntree=10)
#' # investigate variable relations
#' allvariables=names(res$info$depth)
#' rel=var.relations(trees=res$trees,variables=c("X1","X7"),allvariables=allvariables,candidates=allvariables[1:100],t=5)
#' rel$var
#'
#' ###### investigate variable relations without performing variable selection using SMD
#'  # get trees and variable names
#'  x=SMD_example_data[,2:ncol(SMD_example_data)]
#'  y=SMD_example_data[,1]
#'  allvariables=colnames(x)# extract variables names
#'  nvar=length(allvariables)   # count number of variables
#'  set.seed(42)
#'  RF=ranger::ranger(data=SMD_example_data,dependent.variable.name="y",num.trees=100,keep.inbag = TRUE,mtry=floor(nvar^(3/4)),min.node.size=1)
#'  trees=getTreeranger(RF=RF,ntree=100)
#'  trees.lay=addLayer(trees)
#'  trees.surr=addSurrogates(RF=RF,trees=trees.lay,s=10,Xdata=x)
#'
#'  # investigate variable relations
#'  rel=var.relations(trees=trees.surr,variable=c("X1","X7"),allvariables=allvariables,candidates=allvariables[1:100],t=5)
#'  rel$var
#' }
#'
#' @export

var.relations = function(trees,variables,allvariables,candidates,t=5) {

  if (all(candidates %in% allvariables)){
  if (all(variables %in% allvariables)){
  # count surrogates
  s=count.surrogates(trees)
  results.meanAdjAgree=meanAdjAgree(trees,variables,allvariables,candidates,t=t,s$s.a)
  } else {
    stop("allvariables do not contain the chosen variables")
  }
  } else {
    stop("allvariables do not contain the candidate variables")
  }
  surr.var=results.meanAdjAgree$surr.var
  varlist=list()
  for (i in 1:nrow(surr.var)){
    surr.var.var=surr.var[i,]
    if (anyNA(surr.var.var)){
    surr.var.var=surr.var.var[-which(is.na(surr.var.var))]
    }
    var=names(surr.var.var[surr.var.var==1])
    name=variables[i]
    varlist[[name]]=var
  }

  var=names(surr.var[surr.var==1])
  return(list(variables=results.meanAdjAgree$variable, surr.res=results.meanAdjAgree$surr.res,threshold=results.meanAdjAgree$threshold, var=varlist))
}