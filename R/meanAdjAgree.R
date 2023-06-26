#'Calculate mean adjusted agreement to investigate variables relations
#'
#'This is the main function of var.relations function.
#'
#' @param trees list of trees created by getTreeranger, addLayer and addSurrogate.
#' @param variables vector of variable names.
#' @param allvariables vector of all variable names (strings)
#' @param candidates vector of variable names (strings) that are candidates to be related to the variables (has to be contained in allvariables)
#' @param t variable to calculate threshold. Default is 3.
#' @param s.a average number of surrogate variables (ideally calculated by count.surrogates function).
#' @param select.var set False if only relations should be calculated and no related variables should be selected.
#' @param num.threads number of threads used for parallel execution. Default is number of CPUs available.
#' @return a list containing:
#' \itemize{
#' \item variables: the variables to which relations are investigated
#' \item surr.res: matrix with mean adjusted agreement values and variables investigated in rows and candidate variables in columns
#' \item threshold: the threshold used to create surr.var from surr.res
#' \item surr.var: binary matrix showing if the variables are related (1) or non-related (0) with variables in rows and candidates in columns.
#' }
#' @export


meanAdjAgree=function(trees,variables,allvariables,candidates,t,s.a,select.var,num.threads = NULL){
  num.trees=length(trees)
  index.variables=match(variables,allvariables)
  index.candidates = match(candidates,allvariables)
  if (is.null(num.threads)) {
    num.threads = parallel::detectCores()
  }
  list.res = rlist::list.flatten(parallel::mclapply(trees,
                                                    surr.tree,
                                                    mc.cores = num.threads,
                                                    variables,
                                                    index.variables,
                                                    allvariables,
                                                    index.candidates))

  results.allvar = matrix(unlist(lapply(1:length(index.variables),
                                        mean.index,
                                        list.res,
                                        index.variables)),
                          ncol=length(candidates),nrow=length(variables),byrow = TRUE)
  colnames(results.allvar)=candidates
  rownames(results.allvar)=variables

  if(select.var) {
    # calculate threshold and select variables according to it
    adj.mean=mean(unlist(lapply((1:num.trees),adj.mean.trees,trees)),na.rm = TRUE)
    threshold=((s.a*adj.mean)/(length(allvariables)-1))*t
    SurrVar=ifelse(results.allvar>threshold, 1, 0)
    result=list(surr.res=results.allvar,threshold=threshold,surr.var=SurrVar,variables=variables)
  } else {
    result=list(surr.res=results.allvar,variables=variables)
  }
  return(result)
}

#' mean.index
#'
#' This is an internal function
#'
#' @keywords internal
mean.index=function(i, list.res,index.variables){
  list = list.res[which(names(list.res) == index.variables[i])]
  mean.list = round(Reduce("+",list)/length(list),2)
  if (length(mean.list) > 0) {
  mean.list[index.variables[i]] = NA
  return(mean.list)
  } else {
  return(rep(NA,length(index.variables)))
  }
}

#' surr.tree
#'
#' This is an internal function
#'
#' @keywords internal
surr.tree=function(tree,variables,index.variables,allvariables,index.candidates){
  allvar.num = length(allvariables)
  nonterminal.nodes = tree[which(sapply(tree,"[[","status")==1)]
  relevant.nodes = nonterminal.nodes[sapply(nonterminal.nodes,"[[","splitvariable") %in% index.variables]
  if (length(relevant.nodes) > 0) {
    list.nodes = lapply(1:length(relevant.nodes),adj.node,allvar.num,relevant.nodes,index.candidates)
    splitvar = sapply(relevant.nodes,"[[","splitvariable")
    names(list.nodes) = splitvar
    return(list.nodes)
  }
}


#' adj.node
#'
#' This is an internal function
#'
#' @keywords internal
adj.node = function(i,allvar.num,relevant.nodes,index.candidates) {
  node = relevant.nodes[i]
  adjnode = rep(0,allvar.num)
  surr=unlist(sapply(node,"[",-c(1:7))) # extract surrogates
  if ((length(node[[1]]))>7){
    s=(length(surr))/2
    adjnode[surr[1:s]]=surr[(s+1):(2*s)]
  }
  return(adjnode[index.candidates])
}



#' adj.mean
#'
#' This is an internal function
#'
#' @keywords internal
adj.mean=function(trees){
  adj.trees=sapply(1:length(trees),adj.mean.trees,trees)

}

#' adj.mean.trees
#'
#' This is an internal function
#'
#' @keywords internal
adj.mean.trees=function(t,trees){
  tree=trees[[t]]
  nonterminal.nodes=tree[which(sapply(tree,"[[","status")==1)]
  surr.nonterminal=lapply(nonterminal.nodes,"[",-c(1:7))
  adj.tree=mean(unlist(lapply(1:length(surr.nonterminal),mean.adj.node,surr.nonterminal)),na.rm = TRUE)
  if (adj.tree == "NaN") {
    adj.tree = NA
  }
  return(adj.tree)
}

#' adj.node
#'
#' This is an internal function
#'
#' @keywords internal
mean.adj.node=function(m,surr.nonterminal){
  surr=surr.nonterminal[[m]]
  if (length(surr)!=0){
    num.surr=length(surr)/2
    adj=surr[(num.surr+1):(2*num.surr)]
  }
  if (length(surr)==0){
    adj=NA
  }
  return(adj)
}
