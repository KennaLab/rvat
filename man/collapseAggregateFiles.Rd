% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allGenerics.R, R/aggregateFile.R
\name{collapseAggregateFiles}
\alias{collapseAggregateFiles}
\alias{collapseAggregateFiles,aggregateFileList-method}
\title{Collapse multiple aggregate files}
\usage{
collapseAggregateFiles(object, output = NULL, verbose = TRUE)
}
\arguments{
\item{object}{an \code{\link{aggregateFileList}} object.}

\item{output}{Output file name (output will be gz compressed text).
Defaults to \code{NULL}, in which case a data.frame will be returned.}

\item{verbose}{Should the function be verbose? Defaults to \code{TRUE}.}
}
\description{
Collapse aggregrateFiles by aggregating values across aggregateFiles. This will result in one aggregate score for
each sample, representing the aggregate value across aggregate files.
The output will be a two-column matrix including sample IDs and aggregate scores respectively.
}
\examples{
library(rvatData)
gdb <- gdb(rvat_example("rvatData.gdb"))

# generate two aggregate files
varsetfile <- varSetFile(rvat_example("rvatData_varsetfile.txt.gz"))
aggregatefile1 <- tempfile()
aggregate(x = gdb,
          varSet = getVarSet(varsetfile, unit = c("SOD1", "FUS"), varSetName = "High"),
          maxMAF = 0.001,
          output = aggregatefile1,
          verbose = FALSE)

aggregatefile2 <- tempfile()
aggregate(x = gdb,
          varSet = getVarSet(varsetfile, unit = c("NEK1"), varSetName = "High"),
          maxMAF = 0.001,
          output = aggregatefile2,
          verbose = FALSE)

# collapse aggregatefiles
aggregatefile <- tempfile()
collapseAggregateFiles(
  aggregateFileList(c(aggregatefile1, aggregatefile2)),
  output = aggregatefile
)
aggregates <- read.table(aggregatefile, header = TRUE)
head(aggregates)
}
