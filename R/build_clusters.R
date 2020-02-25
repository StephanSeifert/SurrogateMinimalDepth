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
#' x = SMD_example_data[,2:ncol(SMD_example_data)]
#' y = SMD_example_data[,1]
#'
#' \donttest{
#'  # apply SMD to obtain random forest with surrogate variables (usually more trees are needed)
#'  set.seed(42)
#'  res.smd = var.select.smd(x = x, y = y, ntree = 10)
#'
#'  # investigate variable relations
#'  variables = res.smd$forest$allvariables
#'  rel=var.relations(forest = res.smd$forest, variables = variables, candidates = variables, t = 10)
#'  # built groups based on mean adjusted agreement
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
