#'Get a list of structured trees for ranger
#'
#'This functions creates a list of trees for ranger objects similar as getTree function does for random Forest objects.
#'
#' @param RF random forest object created by ranger (with keep.inbag=TRUE)
#' @param ntree number of trees
#' @return a list with trees. Each row of the list elements corresponds to a node of the respective tree and the columns correspond to:
#' \itemize{
#' \item nodeID: ID of the respective node (important for left and right daughters in the next columns)
#' \item leftdaughter: ID of the left daughter of this node
#' \item rightdaughter: ID of the right daughter of this node
#' \item splitvariable: ID of the split variable
#' \item splitpoint: splitpoint of the split variable
#' \item status: "0" for terminal and "1" for non-terminal
#' }
#' @export


getTreeranger=function(RF,ntree) {
trees=lapply(1:ntree,getTreeranger_k,RF=RF)
return(trees)
}

#' getTreeranger_k
#'
#' This is an internal function
#'
#' @keywords internal
getTreeranger_k=function(RF,k=1){
  trees=RF$forest
  split.ids=trees$split.varIDs[[k]]+1
  status=rep(1,length(split.ids))
  child.nodes.left=trees$child.nodeIDs[[k]][[1]]+1
  terminal=which(child.nodes.left == 1)
  child.nodes.left[terminal]=0
  child.nodes.right=trees$child.nodeIDs[[k]][[2]]+1
  child.nodes.right[terminal]=0
  status[terminal]=0
  split.values=as.matrix(trees$split.values[[k]])
  id=c(1:length(split.ids))
  ktree=cbind(id,child.nodes.left,child.nodes.right,split.ids,split.values,status)

  colnames(ktree)=c("nodeID","leftdaughter","rightdaughter","splitvariable","splitpoint","status")
  return(ktree)
}
