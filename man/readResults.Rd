% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rvatResult.R
\name{readResults}
\alias{readResults}
\title{Read association results}
\usage{
readResults(path, header = TRUE, type = NULL, sep = "\\t")
}
\arguments{
\item{path}{File path to data}

\item{header}{A logical value indicating whether the data contains a header. Defaults to \code{TRUE}.}

\item{type}{Result type ('singlevarResult', 'rvbResult', 'gsaResult').
Defaults to \code{NULL} in which case the result type is inferred.}

\item{sep}{The field separator. Defaults to \verb{\\\\t}, which is the default separator using in \code{\link{assocTest}}.}
}
\value{
An object of type \code{\link{rvbResult}} or \code{\link{singlevarResult}}.
}
\description{
Read results generated using the \code{\link{assocTest}} method.
}
