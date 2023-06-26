#'Add layer information to a forest that was created by getTreeranger
#'
#'This functions adds the layer information to each node in a list with trees that was obtained by getTreeranger.
#'
#' @param trees list of trees created by getTreeranger
#' @return a list with trees. Each row of the list elements corresponds to a node of the respective tree and the columns correspond to:
#' \itemize{
#' \item nodeID: ID of the respective node (important for left and right daughters in the next columns)
#' \item leftdaughter: ID of the left daughter of this node
#' \item rightdaughter: ID of the right daughter of this node
#' \item splitvariable: ID of the split variable
#' \item splitpoint: splitpoint of the split variable
#' \item status: "0" for terminal and "1" for non-terminal
#' \item layer: layer information (0 means root node, 1 means 1 layer below root, etc)
#' }
#' @export

addLayer=function(trees){
  #This function adds the respective layer to the different nodes in a tree. The tree has to be prepared by getTree function
  tree.layer=list()
  num.trees= length(trees)
  for (i in 1:num.trees){
    tree=trees[[i]]
    layer=rep(NA,nrow(tree))
    layer[1]=0
    t=1
    while (anyNA(layer)){
      r=unlist(tree[which(layer==(t-1)),2:3])
      layer[r]=t
      t=t+1
    }
    tree=cbind(tree,layer)
    tree=tree[order(as.numeric(tree[,"layer"])),]
    tree.layer[[i]]=tree
  }
  return(tree.layer)
}

