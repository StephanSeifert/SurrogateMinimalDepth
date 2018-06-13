#' Variable selection with Minimal Depth (MD)
#'
#' This function executes MD applying \link[ranger]{ranger} for random forests generation and is a reimplementation of \link[randomForestSRC]{var.select} from randomForestSRC package.
#'
#' @param x matrix or data.frame of predictor variables with variables in
#'   columns and samples in rows (Note: missing values are not allowed)
#' @param y vector with values of phenotype variable (Note: will be converted to factor if
#'   classification mode is used). For survival forests this is the time variable.
#' @param ntree Number of trees. Default is 500.
#' @param mtry Number of variables to possibly split at in each node. Default is no. of variables^(3/4) as recommended by Ishwaran.
#' @param type Mode of prediction ("regression" or "classification"). Default is regression.
#' @param min.node.size Minimal node size. Default is 1.
#' @param num.threads number of threads used for parallel execution. Default is number of CPUs available.
#' @param status status variable, only applicable to survival data. Use 1 for event and 0 for censoring.
#' @param save.ranger Set TRUE if ranger object should be saved. Default is that ranger object is not saved (FALSE).
#'
#' @return List with the following components:
#' \itemize{
#' \item info: list with results from mindep function:
#' \itemize{
#' \item depth: mean minimal depth for each variable
#' \item selected: variables has been selected (1) or not (0),
#' \item threshold: the threshold that is used for the selection
#' }
#' \item var: vector of selected variables
#'
#' \item trees: list of trees that was created by getTreeranger, and addLayer functions and that was used for minimal depth variable importance
#'
#' \item ranger: ranger object
#'
#'}
#' @examples
#' # read data
#' data("SMD_example_data")
#'
#' \donttest{
#' # select variables (usually more trees are needed)
#' res = var.select.md(x=SMD_example_data[,2:ncol(SMD_example_data)],y=SMD_example_data[,1],ntree=10)
#' res$var
#' }
#'
#'@references
##' \itemize{
##'   \item Ishwaran, H. et al. (2011) Random survival forests for high-dimensional data. Stat Anal Data Min, 4, 115–132. \url{https://onlinelibrary.wiley.com/doi/abs/10.1002/sam.10103}
##'   \item Ishwaran, H. et al. (2010) High-Dimensional Variable Selection for Survival Data. J. Am. Stat. Assoc., 105, 205–217. \url{http://www.ccs.miami.edu/~hishwaran/papers/IKGML.JASA.2010.pdf}
##'   }
#'
#' @export

var.select.md = function(x,y,ntree = 500,type = "regression",mtry=NULL,min.node.size=1,num.threads=NULL,status=NULL,save.ranger=FALSE) {

  if (any(is.na(x))) {
    stop("missing values are not allowed")
  }

  variables=colnames(x)# extract variables names
  nvar=length(variables)   # count number of variables

  if (is.null(mtry)) {
    mtry=floor(nvar^(3/4))
  }

  if (type == "classification") {
    y = as.factor(y)
    if (length(levels(y))>15) {
      stop("Too much classes defined, classification might be the wrong choice")
    }
  }
  if (type == "regression" && class(y) == "factor"){
    stop("use factor variable for y only for classification! ")
  }
  data = data.frame(y, x)
  if (type=="survival"){
    if (is.null(status)){
      stop("a status variables has to be given for survival analysis")
    }
    data$status=status
    RF=ranger::ranger(data=data,dependent.variable.name="y",num.trees=ntree,mtry=mtry,min.node.size=min.node.size,keep.inbag = TRUE,
                      num.threads=num.threads, dependent.variable.name="status")
  }
  if (type=="classification" | type=="regression"){
    RF=ranger::ranger(data=data,dependent.variable.name="y",num.trees=ntree,mtry=mtry,min.node.size=min.node.size,keep.inbag = TRUE,
                      num.threads=num.threads)
  }
  trees=getTreeranger(RF=RF,ntree=ntree)
  trees.lay=addLayer(trees)
  minimaldepth=mindep(variables,trees.lay)
  if(save.ranger){
    results=list(info=minimaldepth,var=names(minimaldepth$selected[minimaldepth$selected == 1]),trees=trees.lay,ranger=RF)
  }
  else {
  results=list(info=minimaldepth,var=names(minimaldepth$selected[minimaldepth$selected == 1]),trees=trees.lay)
  }
  return(results)
}
