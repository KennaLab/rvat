% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allGenerics.R
\name{uploadCohort}
\alias{uploadCohort}
\title{Upload cohort table into gdb.}
\usage{
uploadCohort(object, name, value, sep = "\\t")
}
\arguments{
\item{object}{A \code{\link{gdb}} object.}

\item{name}{Name to assign to cohort.}

\item{value}{Input data frame or a valid file path. Must contain an 'IID' column matching to SM table and a 'sex' column (0=missing,1=male,2=female).}

\item{sep}{Field delimiter (applies only when value is a text file). Defaults to \verb{\\\\t}.}
}
\description{
Function to upload cohort data tables to gdb. These will automatically be reformatted and sorted to match the ordering of samples in the gdb genotype records.
}
