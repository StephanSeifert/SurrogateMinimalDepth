#'Execute minimal depth variable importance
#'
#'This function determines the minimal depth of variables from a forest that is created by getTreeranger, and addLayer functions.
#'
#' @param variables vector of variable names
#' @param trees list of trees that was generated by getTreeranger and layers functions
#' @return List with the following components:
#' \itemize{
#' \item depth: mean minimal depth for each variable
#' \item selected: variables has been selected (1) or not (0),
#' \item threshold: the threshold that is used for the selection
#' }
#' @export

mindep=function(variables,trees){

  ntree= length(trees)
# This function determines the minimal depth of variables from a tree with layers that is created by getTree and addLayer functions.
# It is based on this paper: Ishwaran et. al. Journal of the American Statistical Accociation 2010 and used the conservative marginal
# approximation to the minimal depth distribution

var.num=length(variables)

#prepare matrix for mindepth
mindepth=matrix(NA,nrow=ntree,ncol=var.num)
colnames(mindepth)=variables

MAX.DEPTH=10000
#maximal Depth of trees and the number of nodes in every layer is saved to calculate treshold in a following step
maxdepth=rep(NA,ntree)
nodesAtDepthMatrix <- matrix(NA, nrow = MAX.DEPTH, ncol = ntree)
# get mindepth for every variable in every tree
for (i in 1:ntree){
nodesAtDepth <- rep(NA, MAX.DEPTH)
tree=trees[[i]]
tree=tree[order(as.numeric(tree[,"layer"])),]
# get layer information of the variables and save it in minimal depth file
depth.tree=rep(NA,length(variables))
  o=1
  while (anyNA(depth.tree) && o<=nrow(tree)){
    if (tree[o,"status"]==1){
    if (is.na(depth.tree[tree[o,"splitvariable"]])) {
      depth.tree[tree[o,"splitvariable"]]=tree[o,"layer"]
      }
  }
  o=o+1
  }
  #variables with no split in the tree get maxdepth+1 as minimal depth
  depth.tree[which(is.na(depth.tree))]=tree[nrow(tree),"layer"]
  #save min and max depth information for every tree
  mindepth[i,]=depth.tree
  maxdepth[i]=tree[nrow(tree),"layer"]
  #find the number of nodes in every layer

  for (u in 1:(maxdepth[i])){

    nodesAtDepth[u]=nrow(subset(tree, tree[,"layer"] == u & tree[,"status"] == 1 ))
  }
  nodesAtDepthMatrix[,i]=nodesAtDepth
}

# create mean values for the minimal depth of different variables
mean.depth=colMeans(mindepth)

# determine the mean depth of an uninformative variable called D star in Ishwaran et. al. Journal of the American Statistical
# Accociation 2010 (treshhold for important variables)

#the following determination of the treshold is taken from max.subtree of the randomForestSRC package
treeHeight =maxdepth
avgTreeHeight = mean(treeHeight, na.rm=TRUE)
maxTreeHeight = max(treeHeight, na.rm=TRUE)
nodes.at.depth.avg = apply(nodesAtDepthMatrix, 1, mean, na.rm = TRUE)
l = nodes.at.depth.avg[1:avgTreeHeight]
minDepthStatObj = randomForestSRC:::minDepthStat(var.num, l = l)
threshold <- minDepthStatObj$first.moment


# Decide if variables are important or unimportant
Importances=rep(NA,var.num)
for (p in 1:var.num){
if (mean.depth[p]<threshold) {
 Importances[p]=1}
  else {
    Importances[p]=0}
}
names(Importances)=variables
results=list(depth=mean.depth,selected=Importances,threshold=threshold)
return(results)
}

