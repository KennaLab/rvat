% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allGenerics.R, R/spatialClust.R
\name{spatialClust}
\alias{spatialClust}
\alias{spatialClust,gdb-method}
\title{Generate variant sets based on spatial clustering}
\usage{
spatialClust(
  object,
  output,
  varSetName,
  unitTable,
  unitName,
  windowSize,
  overlap,
  intersection = NULL,
  where = NULL,
  weightName = "1",
  posField = "POS",
  minTry = 5,
  warning = TRUE
)
}
\arguments{
\item{object}{a \code{\link{gdb}} object}

\item{output}{Output file name (output will be gz compressed text).}

\item{varSetName}{Name to assign varSet grouping. This identifier column is used to allow for subsequent mergeing of multiple varSet files for coordinated analysis of multiple variant filtering/weighting strategies)}

\item{unitTable}{Table containing aggregation unit mappings.}

\item{unitName}{Field to utilize for aggregation unit names.}

\item{windowSize}{Numeric vector to indicate starting fixed window size (number of variants)}

\item{overlap}{Numeric vector to indicate starting fixed window overlap (number of variants, length must match windowSize)}

\item{intersection}{Additional tables to filter through intersection (ie variants absent from intersection tables will not appear in output). Multiple tables should be ',' delimited.}

\item{where}{An SQL compliant where clause to filter output; eg: "CHROM=2 AND POS between 5000 AND 50000 AND AF<0.01 AND (cadd.caddPhred>15 OR snpEff.SIFT='D')".}

\item{weightName}{Field name for desired variant weighting, must be a column within unitTable or other intersection table. Default value of 1 is equivalent to no weighting.}

\item{posField}{Column name to take as variants position. Default is 'POS' which typically corresponds to genomics position. Can be reset to use CDS or other coordinates. "HGVSc" is a recognized identifier and CDS coordinates will be extracted automatically.}

\item{minTry}{Minimum number of variants in varset to perform clustering on. If number of variants < minTry, all variants will be returned as a single cluster.}

\item{warning}{Raise a warning when clusters can't be generated? Defaults to \code{TRUE}.
Defaults to 5.}
}
\description{
Generate weighted variant sets for use in association testing, with partitioning by genomic distances as described (Fier, GenetEpidemiol, 2017).
}
\references{
Loehlein Fier, H. et al. On the association analysis of genome-sequencing data: A spatial clustering approach for partitioning the entire genome into nonoverlapping windows: F ier et al . Genet. Epidemiol. 41, 332–340 (2017).
}
