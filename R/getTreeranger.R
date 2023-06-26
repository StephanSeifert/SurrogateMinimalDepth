#'Get a list of structured trees for ranger
#'
#'This functions creates a list of trees for ranger objects similar as getTree function does for random Forest objects.
#'
#' @param RF random forest object created by ranger (with keep.inbag=TRUE)
#' @param num.trees number of trees
#' @return a list with trees. Each row of the list elements corresponds to a node of the respective tree and the columns correspond to:
#' \itemize{
#' \item nodeID: ID of the respective node (important for left and right daughters in the next columns)
#' \item leftdaughter: ID of the left daughter of this node
#' \item rightdaughter: ID of the right daughter of this node
#' \item splitvariable: ID of the split variable
#' \item splitpoint: splitpoint of the split variable (for categorical variables this is a comma separated lists of values, representing the factor levels (in the original order) going to the right)
#' \item status: "0" for terminal and "1" for non-terminal
#' }
#' @export



getTreeranger=function(RF,num.trees) {
  trees=lapply(1:num.trees,getsingletree,RF=RF)

  return(trees)
}


#' getsingletree
#'
#' This is an internal function
#'
#' @keywords internal
getsingletree=function(RF,k=1){
  # here we use the treeInfo function of the ranger package to create extract the trees, in an earlier version this was done with a self implemented function
  tree.ranger = ranger::treeInfo(RF,tree = k)
  ktree=data.frame(as.numeric(tree.ranger$nodeID+1),
              as.numeric(tree.ranger$leftChild+1),
              as.numeric(tree.ranger$rightChild+1),
              as.numeric(tree.ranger$splitvarID+1),
              tree.ranger$splitval,
              tree.ranger$terminal)
  if (is.factor(ktree[,5])) {
  ktree[,5] = as.character(levels(ktree[,5] ))[ktree[,5]]
  }
  ktree[,6] = as.numeric(ktree[,6] == FALSE)

 for (i in 2:4) {
   ktree[,i][is.na(ktree[,i])] = 0
 }
  colnames(ktree)=c("nodeID","leftdaughter","rightdaughter","splitvariable","splitpoint","status")
  return(ktree)
}
