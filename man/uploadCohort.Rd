% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allGenerics.R, R/gdb.R
\name{uploadCohort}
\alias{uploadCohort}
\alias{uploadCohort,gdb-method}
\title{Upload cohort table into gdb.}
\usage{
uploadCohort(object, name, value, sep = "\\t", verbose = TRUE)
}
\arguments{
\item{object}{A \code{\link{gdb}} object.}

\item{name}{Name to assign to cohort.}

\item{value}{Input data frame or a valid file path. Must contain an 'IID' column matching to SM table and a 'sex' column (0=missing,1=male,2=female).}

\item{sep}{Field delimiter (applies only when value is a text file). Defaults to \verb{\\\\t}.}

\item{verbose}{Should the method be verbose? Defaults to \code{TRUE}.}
}
\description{
Function to upload cohort data tables to gdb. These will automatically be reformatted and sorted to match the ordering of samples in the gdb genotype records.
}
\examples{
library(rvatData)
gdb <- create_example_gdb()

# from data.frame
pheno <- read.table(rvat_example("rvatData.pheno"), header = TRUE)
uploadCohort(object = gdb, name = "cohortinfo", value = pheno)

# similarly, a cohort table can be imported directly from file
filepath <- rvat_example("rvatData.pheno")
uploadCohort(object = gdb, name = "cohortinfo2", value = filepath)

}
