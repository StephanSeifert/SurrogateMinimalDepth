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
#' @param num.threads number of threads used for parallel execution. Default is number of CPUs available.
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

addSurrogates = function(RF,trees,s,Xdata,num.threads) {

  num.trees = length(trees)
  ncat = sapply(sapply(Xdata,levels),length)     # determine number of categories (o for continuous variables)
  names(ncat) = colnames(Xdata)

  if (is.null(num.threads)) {
    num.threads = parallel::detectCores()
  }

  if (any(ncat) > 0) {
  Xdata[,which(ncat > 0)] = sapply(Xdata[,which(ncat > 0)],unclass)
  }

  #variables to find surrogates (control file similar as in rpart)
  controls = list(maxsurrogate = as.integer(s), sur_agree = 0)

  trees.surr = parallel::mclapply(1:num.trees,
                                  getSurrogate,
                                  mc.cores = num.threads,
                                  maxsurr = s,
                                  surr.par = list(inbag.counts = RF$inbag.counts,
                                                         Xdata = Xdata,
                                                      controls = controls,
                                                         trees = trees,
                                                          ncat = ncat))
  return(trees.surr)
}

#' getSurrogate
#'
#' This is an internal function
#'
#' @keywords internal
getSurrogate = function(surr.par, k = 1, maxsurr) {
  #weights and trees are extracted
 tree = surr.par$trees[[k]]
 column.names = colnames(tree)
 n.nodes = nrow(tree)
 wt = surr.par$inbag.counts[[k]]
 tree.surr = lapply(1:n.nodes,
                    SurrTree,
                    wt = wt,
                    Xdata = surr.par$Xdata,
                    controls = surr.par$controls,
                    column.names, tree,maxsurr,
                    ncat = surr.par$ncat)
}
#' SurrTree
#'
#' This is an internal function
#'
#' @keywords internal
SurrTree = function(j,wt,Xdata,controls,column.names,tree,maxsurr,ncat) {
  node = tree[j,]
  # for non-terminal nodes get surrogates
  if (node["status"] == 1) {
  #Handover to C
  var = as.numeric(node[4]) # extract split variable

  if (ncat[var] == 0) { # extract split information: split point for continuous variables and directions for qualitative variables
    split = as.numeric(node[5])
  } else {
    right = as.numeric(strsplit(as.character(node[5]), ",")[[1]])
    directions = rep(-1,ncat[var])
    directions[right] = 1
    split = as.numeric(c(ncat[var],directions))
  }


  surrogate.parameters = .Call(C_getSurrogates,
                               ncat = as.integer(ncat),
                               wt = as.numeric(wt),
                               X = as.matrix(Xdata),
                               controls = as.integer(unlist(controls)),
                               var = as.integer(var),                      # node variables
                               split = as.numeric(split))                    # split info

  if (nrow(surrogate.parameters$isplit) > 1) {
    surrogates = surrogate.parameters$isplit[2:nrow(surrogate.parameters$isplit),1]
    surr.adj = round(surrogate.parameters$dsplit[2:nrow(surrogate.parameters$dsplit),1],2)
    node.new = data.frame(matrix(nrow = 1, ncol = 7 + length(surrogates) + length(surr.adj)))
    node.new[,1:7] = node[1:7]
    node.new[,8:(7 + length(surrogates) + length(surr.adj))] = c(surrogates,surr.adj)
    surrogate.names = NULL
    adj.names = NULL
    surrogate.names = sapply(1:length(surrogates),name.surr,surrogate.names)
    adj.names = sapply(1:length(surrogates),name.adj,adj.names)
    names(node.new) = c(column.names,surrogate.names,adj.names)
  }

  if (nrow(surrogate.parameters$isplit) == 1) {
    node.new = node
  }
  }
  if (node["status"] == 0) {
    node.new = node
  }
return(node.new)
}

#' name.surr
#'
#' This is an internal function
#'
#' @keywords internal
name.surr = function(i,surrogate.names){
surrogate.names = c(surrogate.names,paste0("surrogate_",i))
return(surrogate.names)
}

#' name.adj
#'
#' This is an internal function
#'
#' @keywords internal
name.adj = function(i,adj.names){
  adj.names = c(adj.names,paste0("adj_",i))
  return(adj.names)
}

