% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allGenerics.R
\name{mapToMatrix}
\alias{mapToMatrix}
\title{mapToMatrix}
\usage{
mapToMatrix(object, results, ID = "unit", sparse = TRUE)
}
\arguments{
\item{object}{a \code{\link{geneSetList}} object}

\item{results}{a \code{\link{rvbResult}} object}

\item{ID}{Name of the column in \code{results} that corresponds to the IDs in the geneSetList. Defaults to 'unit'.}

\item{sparse}{Return a sparse matrix? Defaults to \code{TRUE}.}
}
\description{
Map the units in a \code{\link{rvbResult}} to the geneSets contained in a
\code{\link{geneSetList}} object.
Returns a matrix where the number of rows equals the number of rows in the rvbResults
and the columns equal the number of geneSets in the geneSetList.
}
