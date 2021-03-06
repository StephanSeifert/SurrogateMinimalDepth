% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/var.relations.R
\name{var.relations}
\alias{var.relations}
\title{Investigate variable relations of a specific variable with mean adjusted agreement}
\usage{
var.relations(forest, variables, candidates, t = 5)
}
\arguments{
\item{forest}{a list containing allvariables and trees. Allvariables is a vector of all variable names in the original data set (strings). Trees is a list of trees that was generated by getTreeranger, addLayer, and getSurrogates functions.}

\item{variables}{variable names (string) for which related variables should be searched for (has to be contained in allvariables)}

\item{candidates}{vector of variable names (strings) that are candidates to be related to the variables (has to be contained in allvariables)}

\item{t}{variable to calculate threshold. Default is 5.}
}
\value{
a list containing:
\itemize{
\item variables: the variables to which relations are investigated.
\item surr.res: a matrix with mean adjusted agreement values with variables in rows and candidates in columns.
\item threshold: the threshold used to select related variables.
\item var: a list with one vector for each variable containing related variables.
}
}
\description{
This function uses the mean adjusted agreement to select variables that are related to a defined variable using a threshold T.
The parameter t is used to calculate T: t=1 means that every variable with higher probability than "by chance" is identified
as "important". t=2 means the probability has to be twice, etc.
Based on the threshold a vector is created containing the related variables. In order to get all variable relations
}
\examples{
# read data
data("SMD_example_data")
x = SMD_example_data[,2:ncol(SMD_example_data)]
y = SMD_example_data[,1]

\donttest{
# apply SMD to obtain random forest with surrogate variables (usually more trees are needed)
set.seed(42)
res = var.select.smd(x = x, y = y, s = 10, ntree = 10)

# investigate variable relations
candidates = res$forest[["allvariables"]][1:100]
rel = var.relations(forest = res$forest, variables=c("X1","X7"), candidates = candidates, t = 5)
rel$var

}

}
