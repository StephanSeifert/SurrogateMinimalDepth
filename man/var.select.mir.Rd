% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/variable_selection_mir.R
\name{var.select.mir}
\alias{var.select.mir}
\title{Variable selection with mutual impurity reduction (MIR)}
\usage{
var.select.mir(
  x = NULL,
  y = NULL,
  ntree = 500,
  type = "regression",
  s = NULL,
  mtry = NULL,
  min.node.size = 1,
  num.threads = NULL,
  status = NULL,
  save.ranger = FALSE,
  save.memory = FALSE,
  min.var.p = 200,
  p.t.sel = 0.01,
  p.t.rel = 0.01,
  select.var = TRUE,
  select.rel = FALSE,
  case.weights = NULL,
  corr.rel = TRUE,
  t = 5,
  method.rel = "janitza",
  method.sel = "janitza",
  num.threads.rel = NULL
)
}
\arguments{
\item{x}{matrix or data.frame of predictor variables with variables in
columns and samples in rows (Note: missing values are not allowed)}

\item{y}{vector with values of phenotype variable (Note: will be converted to factor if
classification mode is used). For survival forests this is the time variable.}

\item{ntree}{number of trees. Default is 500.}

\item{type}{mode of prediction ("regression" or "classification"). Default is regression.}

\item{s}{predefined number of surrogate splits (it may happen that the actual number of surrogate splits differs in individual nodes). Default is 1 \% of no. of variables.}

\item{mtry}{number of variables to possibly split at in each node. Default is no. of variables^(3/4) ("^3/4") as recommended by (Ishwaran 2011). Also possible is "sqrt" and "0.5" to use the square root or half of the no. of variables.}

\item{min.node.size}{minimal node size. Default is 1.}

\item{num.threads}{number of threads used for parallel execution. Default is number of CPUs available.}

\item{status}{status variable, only applicable to survival data. Use 1 for event and 0 for censoring.}

\item{save.ranger}{set TRUE if ranger object should be saved. Default is that ranger object is not saved (FALSE).}

\item{save.memory}{Use memory saving (but slower) splitting mode. No effect for survival and GWAS data. Warning: This option slows down the tree growing, use only if you encounter memory problems. (This parameter is transfered to ranger)}

\item{min.var.p}{minimum number of permuted variables used to determine p-value for variable selection of important variables. Default is 200.}

\item{p.t.sel}{p.value threshold for selection of important variables. Default is 0.01.}

\item{p.t.rel}{p.value threshold for selection of related variables. Default is 0.01.}

\item{select.var}{set False if only importance should be calculated and no variables should be selected.}

\item{select.rel}{set False if only relations should be calculated and no variables should be selected.}

\item{case.weights}{Weights for sampling of training observations. Observations with larger weights will be selected with higher probability in the bootstrap (or subsampled) samples for the trees.}

\item{corr.rel}{set FALSE if non-corrected variable relations should be used for calculation of MIR. In this case the method "janitza" should not be used for selection of important variables}

\item{t}{variable to calculate threshold for non-corrected relation analysis. Default is 5.}

\item{method.rel}{Method  to  compute  p-values for selection of related variables with var.relations.corr. Use  "janitza"  for  the  method  by  Janitza  et  al. (2016) or "permutation" to utilize permuted variables.}

\item{method.sel}{Method  to  compute  p-values for selection of important variables. Use  "janitza"  for  the  method  by  Janitza  et  al. (2016) (can only be used when corrected variable relations are utilized) or "permutation" to utilize permuted variables.}

\item{num.threads.rel}{number of threads used for determination of relations. Default is number of CPUs available. (this process can be memory-intensive and it can be preferable to reduce this)}
}
\value{
list with the following components:
\itemize{
\item info: list with results containing:
\itemize{
\item MIR: the calculated variable importance for each variable based on mutual impurity reduction.
\item pvalue: the obtained p-values for each variable.
\item selected: variables has been selected (1) or not (0).
\item relations: a list containing the results of variable relation analysis.
\item parameters: a list that contains the parameters s, type, mtry, p.t.sel, p.t.rel and method.sel that were used.
}
\item var: vector of selected variables.

\item ranger: ranger object.

}
}
\description{
This function executes MIR applying \link[ranger]{ranger} for random forests generation and actual impurity reduction and a modified version of \link[rpart]{rpart} to find surrogate variables.
}
\examples{
# read data
data("SMD_example_data")

\donttest{
# select variables (usually more trees are needed)
set.seed(42)
res = var.select.mir(x = SMD_example_data[,2:ncol(SMD_example_data)], y = SMD_example_data[,1],s = 10, ntree = 10)
res$var
}
}
\references{
\itemize{
  \item Nembrini, S. et al. (2018) The revival of the Gini importance? Bioinformatics, 34, 3711–3718. \url{https://academic.oup.com/bioinformatics/article/34/21/3711/4994791}
  \item Seifert, S. et al. (2019) Surrogate minimal depth as an importance measure for variables in random forests. Bioinformatics, 35, 3663–3671. \url{https://academic.oup.com/bioinformatics/article/35/19/3663/5368013}
  }
}
