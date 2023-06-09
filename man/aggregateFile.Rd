% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allClasses.R, R/allGenerics.R,
%   R/aggregateFile.R
\docType{class}
\name{aggregateFile}
\alias{aggregateFile}
\alias{listSamples}
\alias{getUnit}
\title{Class to manage interactions with an aggregateFile}
\description{
The \code{\link{aggregate}} method saves genotypes aggregates compressed on disk. The
aggregateFile class manages connecting/interacting with these files and retrieving
units of interest.
}
\section{Build an aggregateFile}{

\itemize{
\item \code{aggregate(x, ...)}: Returns an aggregate of genotypes for each individual.
See \code{\link[=aggregate]{aggregate()}} for details.
}
}

\section{Connect to an aggregateFile}{

\itemize{
\item \code{aggregateFile(path)}: Connect to an aggregateFile object.
}
}

\section{Getters}{

In the following code snippets, x is an aggregateFile object.
\itemize{
\item \code{getUnit(x, unit)}: Retrieve aggregates for specified unit(s).
Use \code{listUnits(x)} to list the units includes in the aggregateFile
Output will be a matrix.
\item \code{listUnits(x)}: Return a vector of all units included in the aggregateFile
\item \code{listSamples(x)}: Return a vector of all sample IDs included in the aggregateFile
}
}

\section{Association testing}{

An aggregateFile can be directly supplied to the \code{\link[=assocTest]{assocTest()}} method.
}

\section{Merging}{

Aggregate files can be merged using the \code{\link{mergeAggregateFiles}} method.
}

\seealso{
\code{\link{mergeAggregateFiles}}

\code{\link{aggregateFileList}}

\code{\link{assocTest-aggregateFile}}

\code{\link{aggregate}}
}
