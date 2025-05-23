% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allGenerics.R, R/rvatResult.R
\name{densityPlot}
\alias{densityPlot}
\alias{densityPlot,rvatResult-method}
\title{Density plot}
\usage{
densityPlot(
  object,
  geneSet,
  geneSetList,
  showMeans = FALSE,
  INT = FALSE,
  Zcutoffs = NULL,
  title = ""
)
}
\arguments{
\item{object}{a \code{\link{rvbResult}} object}

\item{geneSet}{string with the name of the gene set of interest}

\item{geneSetList}{list with the available gene sets}

\item{showMeans}{Show mean Z-score for genes in and outside the gene set. Defaults to \code{FALSE}.}

\item{INT}{Apply inverse normal transformation to Z-scores? Defaults to \code{FALSE}.}

\item{Zcutoffs}{A vector (length=2, minimum and maximum) of cutoffs to apply to the Z-scores.
Z scores below/above these cutoffs will be set equal to the cutoff.}

\item{title}{the title of the plot}
}
\description{
Generate a density plot for a \code{\link{rvatResult}} highlighting the Z-scores of the selected gene set
}
