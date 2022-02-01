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
  ntree=length(trees)
  surr.result=rep(NA,length(allvariables))
  index.variables=match(variables,allvariables)
  if (is.null(num.threads)) {
    num.threads = parallel::detectCores()
  }
  results.allvar = t(sapply(1:length(index.variables),maa.p,allvariables,ntree,trees,index.variables,candidates,num.threads))
  colnames(results.allvar)=candidates
  rownames(results.allvar)=variables
  if(select.var) {
  # calculate threshold and select variables according to it
  adj.mean=mean(unlist(lapply((1:ntree),adj.mean.trees,trees)),na.rm = TRUE)
  threshold=((s.a*adj.mean)/(length(allvariables)-1))*t
  SurrVar=ifelse(results.allvar>threshold, 1, 0)
  result=list(surr.res=results.allvar,threshold=threshold,surr.var=SurrVar,variables=variables)
  } else {
    result=list(surr.res=results.allvar,variables=variables)
  }
  return(result)
}


#' maa.p
#'
#' This is an internal function
#'
#' @keywords internal
maa.p=function(p=1,allvariables,ntree,trees,index.variables,candidates,num.threads){
  i=index.variables[p]
  surrMatrix=matrix(unlist(parallel::mclapply(trees[1:ntree],mc.cores = num.threads,surr.tree,allvariables,ntree,i)),ncol=length(allvariables),nrow=ntree,byrow = TRUE)
  colnames(surrMatrix)=allvariables
  means.surr=colMeans(surrMatrix,na.rm=TRUE)
  means.surr[i]=NA
  means.surr.candidate=means.surr[candidates]
  means.surr.candidate[which(means.surr.candidate == "NaN")] = NA
  return(means.surr.candidate)
}


#' surr.var
#'
#' This is an internal function
#'
#' @keywords internal
surr.var=function(i=1,variables,ntree,trees){

  surrMatrix=t(sapply(1:ntree,surr.tree,variables,ntree,trees,i))
  colnames(surrMatrix)=variables
  means.surr=colMeans(surrMatrix,na.rm=TRUE)
  return(means.surr)
}

#' surr.tree
#'
#' This is an internal function
#'
#' @keywords internal
surr.tree=function(tree,variables,ntree,trees,i){
  adjtree=rep(0,length(variables))
  # there are more than one nonterminal nodes with split variable i
  if (length(which(sapply(tree,"[[",4)==i))>1){
    nodes=tree[which(sapply(tree,"[[",4)==i)]
    s=sapply(nodes,length)
    s=(s-7)/2
    surr=lapply(nodes,"[",-c(1:7)) # extract surrogates
    sum=0
    for (o in 1:length(s)) {
      if (s[o]==0) next
      adjtree.k=rep(0,length(variables))
      surr.var=surr[[o]][(1:s[o])]
      surr.adj=surr[[o]][(s[o]+1):(2*s[o])]
        adjtree.k[surr.var]=surr.adj
        adjtree=adjtree+adjtree.k
        sum=sum+1
      }
      adjtree=adjtree/sum
  }
  #there is one nonterminal node with split variable i
  if (length(which(sapply(tree,"[[",4)==i))==1){
    nodes=tree[which(sapply(tree,"[[",4)==i)]
    surr=sapply(nodes,"[",-c(1:7)) # extract surrogates
      if ((length(nodes[[1]]))>7){
        s=(length(surr))/2
        surr.var=surr[1:s]
        surr.adj=surr[(s+1):(2*s)]
        adjtree[surr.var]=surr.adj
      }
    }
  #there is no nonterminal node with split variable i
  if (length(which(sapply(tree,"[[",4)==i))==0){
    adjtree=rep(NA,length(variables))
    surr.mean=NA
  }
  return(adjtree=adjtree)
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
    adj.tree=mean(unlist(lapply(1:length(surr.nonterminal),adj.node,surr.nonterminal)),na.rm = TRUE)
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
  adj.node=function(m,surr.nonterminal){
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

