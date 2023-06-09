% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assocTest.R
\name{assocTest-genoMatrix}
\alias{assocTest-genoMatrix}
\title{assocTest-genoMatrix}
\usage{
\S4method{assocTest}{genoMatrix}(
  object,
  pheno,
  test,
  name = "none",
  continuous = FALSE,
  singlevar = FALSE,
  covar = NULL,
  offset = NULL,
  overwriteAggregate = TRUE,
  geneticModel = c("allelic", "recessive", "dominant"),
  imputeMethod = NULL,
  MAFweights = "none",
  maxitFirth = 1000,
  keep = NULL,
  output = NULL,
  append = FALSE,
  returnDF = FALSE,
  methodResampling = NULL,
  resamplingMatrix = NULL,
  nResampling = 1000,
  outputResampling = FALSE,
  memlimitResampling = NULL,
  minCallrateVar = 0,
  maxCallrateVar = Inf,
  minCallrateSM = 0,
  maxCallrateSM = Inf,
  minMAF = 0,
  maxMAF = 1,
  minMAC = 0,
  maxMAC = Inf,
  minCarriers = 0,
  maxCarriers = Inf,
  minCarrierFreq = 0,
  maxCarrierFreq = Inf,
  verbose = TRUE
)
}
\arguments{
\item{object}{a \code{\link{genoMatrix}} object}

\item{pheno}{colData field to test as response variable, the response variable
can either be binary (0/1) or continuous. If the response variable is continuous set
\code{continuous} to \code{TRUE}.}

\item{test}{Vector of statistical tests to run,
options include firth,glm,lm,scoreSPA,skat,skat_burden,skato,skat_fwe,skat_burden_fwe,
skato_fwe,skat_robust,skato_robust,skat_burden_robust, acatv, acatvSPA. See \code{\link{assocTest}} for details.}

\item{name}{Optional name for the analysis, defaults to "none".}

\item{continuous}{Is the response variable continuous? (TRUE/FALSE). Defaults to \code{FALSE}.}

\item{singlevar}{Run single variant tests? (TRUE/FALSE).
Defaults to \code{FALSE}, in which case collapsing tests are ran.}

\item{covar}{Character vector of covariates. These should be present in the colData slot of the genoMatrix.}

\item{overwriteAggregate}{In case there is already an \code{aggregate} column in the \code{colData} of the genoMatrix
(i.e. \code{aggregate} has been run on the genoMatrix), should it be overwitten? Defaults to \code{TRUE}.}

\item{geneticModel}{Which genetic model to apply? ('allelic', 'recessive' or 'dominant').
Defaults to \code{allelic}.}

\item{imputeMethod}{Which imputation method to apply? ('meanImpute' or 'missingToRef').
Defaults to \code{meanImpute}}

\item{MAFweights}{Apply MAF weighting? Currently Madsen-Browning ('mb') is implemented.
Defaults to 'none'.}

\item{maxitFirth}{Maximum number of iterations to use for estimating firth confidence intervals.}

\item{output}{Output file path for results.
Defaults to \code{NULL}, in which case results are not written to disk, but returned as an \code{\link{rvatResult}} object.}

\item{append}{Relevant if the \code{output} parameter is not \code{NULL}. Should results be appended to \code{output}?
Defaults to \code{FALSE}.}

\item{methodResampling}{Which method to use for resampling? ('permutation' currently implemented)
Defaults to \code{NULL}, in which case no resampling is performed.}

\item{resamplingMatrix}{Pre-calculated resampling matrix (n x p), where n = number of samples, and p number of resamplings.
Can be generated using \code{\link{buildResamplingFile}}.}

\item{nResampling}{Number of resamplings to perform if methodResampling is specified.}

\item{outputResampling}{If \code{TRUE} or a filepath, results for each resampling are returned (or saved to the filepath).
This can be useful if permutations are used to calculated to estimate correlations among genes for example.
Defaults to \code{FALSE} in which case resampling is used to calculate resampled P-values,
results for individual resamplings are not returned.}

\item{memlimitResampling}{Maximum number of resamplings to perform at a time.
Resampling generates a matrix of n x p, where n is the number of samples and p the number of resamplings
thus, for large number of resamplings it can be more efficient to split the permutations in chunks of size \code{memlimitResampling}.
Defaults to \code{NULL} in which case all permutations are performed.}

\item{minCallrateVar}{Minimum genotype rate for variant retention.}

\item{maxCallrateVar}{Maximum genotype rate for variant retention.}

\item{minCallrateSM}{Minimum genotype rate for sample retention.}

\item{maxCallrateSM}{Maximum genotype rate for sample retention.}

\item{minMAF}{Minimum minor allele frequency for variant retention.}

\item{maxMAF}{Maximum minor allele frequency for variant retention.}

\item{minMAC}{Minimum minor allele count for variant retention.}

\item{maxMAC}{Maximum minor allele count for variant retention.}

\item{minCarriers}{Minimum carrier count for variant retention.}

\item{maxCarriers}{Maximum carrier count for variant retention.}

\item{minCarrierFreq}{Minimum carrier frequency for variant retention.}

\item{maxCarrierFreq}{Maximum carrier frequency for variant retention.}

\item{verbose}{Should the function be verbose? (TRUE/FALSE), defaults to \code{TRUE}.}
}
\description{
Run \code{\link{assocTest}} on a \code{\link{genoMatrix}} object. See the main \code{\link{assocTest}} page for details.
}
