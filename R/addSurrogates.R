#' Add surrogate information that was created by getTreeranger
#'
#' This function adds surrogate variables and adjusted agreement values to a forest that was created by getTreeranger.
#'
#' @useDynLib SurrogateMinimalDepth, .registration = TRUE, .fixes = "C_"
#'
#' @param RF random forest object created by ranger (with keep.inbag=TRUE).
#' @param trees list of trees created by getTreeranger.
#' @param s Predefined number of surrogate splits (it may happen that the actual number of surrogate splits differes in individual nodes). Default is 1 \% of no. of variables.
#' @param Xdata data without the dependent variable.
#' @return a list with trees containing of lists of nodes with the elements:
#' \itemize{
#' \item nodeID: ID of the respective node (important for left and right daughters in the next columns)
#' \item leftdaughter: ID of the left daughter of this node
#' \item rightdaughter: ID of the right daughter of this node
#' \item splitvariable: ID of the split variable
#' \item splitpoint: splitpoint of the split variable
#' \item status: "0" for terminal and "1" for non-terminal
#' \item layer: layer information (0 means root node, 1 means 1 layer below root, etc)
#' \item surrogate_i: numbered surrogate variables (number depending on s)
#' \item adj_i: adjusted agreement of variable i
#' }
#' @export

addSurrogates=function(RF,trees,s,Xdata){

  if (exists("s")) {
    s=s
  } else {
    s=nvar*0.01
  }
  nvar=ncol(Xdata)            # count number of variables
  ntree= length(trees)
  # categories for variables 0 for continuous (has to be defined for surrogate variables searching)
  ncat = rep(0, nvar)
  #variables to find surrogates (control file similar as in rpart)
  controls = list(minsplit = 0L, minbucket = 0L, cp = 0L,
                  maxcompete = 0L, maxsurrogate = as.integer(s), usesurrogate = 0L,
                  surrogatestyle = 0L, maxdepth = 0L, xval = 0L)

  trees.surr=lapply(1:ntree,getSurrogate,maxsurr=s,surr.par=list(ncat = ncat,
                                                                inbag.counts=RF$inbag.counts,
                                                                Xdata=Xdata,
                                                                controls=controls,
                                                                trees=trees))
  return(trees.surr)
}

#' getSurrogate
#'
#' This is an internal function
#'
#' @keywords internal
getSurrogate=function(surr.par,k=1,maxsurr){
  #weights and trees are extracted
 tree=surr.par$trees[[k]]
 column.names=colnames(tree)
 n.nodes=nrow(tree)
 wt=surr.par$inbag.counts[[k]]
 tree.surr=lapply(1:n.nodes,SurrTree,ncat=surr.par$ncat,wt=wt,Xdata=surr.par$Xdata,controls=surr.par$controls,column.names,tree,maxsurr)
}

#' SurrTree
#'
#' This is an internal function
#'
#' @keywords internal
SurrTree=function(j,ncat,wt,Xdata,controls,column.names,tree,maxsurr){
  node=tree[j,]
  # for non-terminal nodes get surrogates
  if (node["status"]==1){
  pnode=node[4:5]
  #Handover to C
  surrogate.parameters= .Call(C_getSurrogates,
                              ncat = as.integer(ncat),
                              wt=as.numeric(wt),
                              X=as.matrix(Xdata),
                              controls=as.integer(unlist(controls)),
                              as.numeric(pnode))
  if (nrow(surrogate.parameters$isplit)>1){
    surrogates=surrogate.parameters$isplit[2:nrow(surrogate.parameters$isplit),1]
    surr.adj=surrogate.parameters$dsplit[2:nrow(surrogate.parameters$dsplit),3]
    node.new=c(node,surrogates,surr.adj)
    surrogate.names=NULL
    adj.names=NULL
    surrogate.names=sapply(1:length(surrogates),name.surr,surrogate.names)
    adj.names=sapply(1:length(surrogates),name.adj,adj.names)
    names(node.new)=c(column.names,surrogate.names,adj.names)
  }

  if (nrow(surrogate.parameters$isplit)==1){
    node.new=node
  }
  }
  if (node["status"]==0){
    node.new=node
  }
return(node.new)
}

#' name.surr
#'
#' This is an internal function
#'
#' @keywords internal
name.surr=function(i,surrogate.names){
surrogate.names=c(surrogate.names,paste0("surrogate_",i))
return(surrogate.names)
}

#' name.adj
#'
#' This is an internal function
#'
#' @keywords internal
name.adj=function(i,adj.names){
  adj.names=c(adj.names,paste0("adj_",i))
  return(adj.names)
}

