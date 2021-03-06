% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/addLayer.R
\name{addLayer}
\alias{addLayer}
\title{Add layer information to a forest that was created by getTreeranger}
\usage{
addLayer(trees)
}
\arguments{
\item{trees}{list of trees created by getTreeranger}
}
\value{
a list with trees. Each row of the list elements corresponds to a node of the respective tree and the columns correspond to:
\itemize{
\item nodeID: ID of the respective node (important for left and right daughters in the next columns)
\item leftdaughter: ID of the left daughter of this node
\item rightdaughter: ID of the right daughter of this node
\item splitvariable: ID of the split variable
\item splitpoint: splitpoint of the split variable
\item status: "0" for terminal and "1" for non-terminal
\item layer: layer information (0 means root node, 1 means 1 layer below root, etc)
}
}
\description{
This functions adds the layer information to each node in a list with trees that was obtained by getTreeranger.
}
