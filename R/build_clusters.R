#' Apply cluster analysis to build variable groups
#'
#'This function generates variables groups of relation information that was obtained by \link[SurrogateMinimalDepth]{var.relations} function applying
#'\link[linkcomm]{getLinkCommunities}.
#'
#' @param rel a list containing variables, surr.res, threshold, and var. This is the output of \link[SurrogateMinimalDepth]{var.relations} function.
#' @param hcmethod the hierarchical clustering method that is used. (see \link[linkcomm]{getLinkCommunities})
#'
#' @return a data frame containing the variable names and their associated clusters.
#'
#' @examples
#' # read data
#' data("SMD_example_data")
#'
#' \donttest{
#'  # get trees and variable names
#'  x = SMD_example_data[,2:ncol(SMD_example_data)]
#'  y = SMD_example_data[,1]
#'  allvariables = colnames(x)# extract variables names
#'  nvar = length(allvariables)   # count number of variables
#'  set.seed(42)
#'  RF = ranger::ranger(data = SMD_example_data, dependent.variable.name = "y", num.trees = 10, keep.inbag = TRUE,mtry = floor(nvar^(3/4)), min.node.size = 1)
#'  trees = getTreeranger(RF = RF, num.trees = 10)
#'  trees.lay = addLayer(trees)
#'  trees.surr = addSurrogates(RF = RF, trees = trees.lay, s = 10, Xdata = x, num.threads = NULL)
#'
#'  # investigate variable relations
#'  rel=var.relations(forest = list(trees = trees.surr, allvariables = allvariables), variables = allvariables, candidates = allvariables, t = 10)
#'  groups = build.clusters(rel)
#' }
#'
#' @export


build.clusters = function(rel,hcmethod = "ward.D") {

  # create links matrix
  links = matrix(nrow = length(unlist(rel$var)) ,ncol = 2)
  links[,1] = unlist(sapply(1:length(rel$variables), function(x) {rep(rel$variables[x],length(rel$var[[x]]))}))
  links[,2] = unlist(rel$var)

  # cluster analysis to identify clusters
  link.comm = linkcomm::getLinkCommunities(links,
                                 hcmethod = hcmethod,
                                 plot = FALSE)
  cluster.info = link.comm$nodeclusters
  return(cluster.info)
}
