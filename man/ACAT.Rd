% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allGenerics.R
\name{ACAT}
\alias{ACAT}
\title{Combine P-values using the ACAT method.}
\usage{
ACAT(
  object,
  aggregate = "test",
  group = c("unit", "cohort", "varSetName", "name", "pheno", "covar", "geneticModel",
    "MAFweight", "test"),
  fixpval = TRUE,
  fixpval_method = c("minmax", "manual", ),
  fixpval_maxP = 0.99,
  fixpval_minP = 1e-32,
  warning = TRUE
)
}
\arguments{
\item{object}{\code{\link{rvatResult}} object}

\item{aggregate}{Variable to ACAT P-values across.
For example, if aggregate = 'test', P-values across statistical tests will be combined using ACAT.
A vector with multiple column names can be specified.
For example, by specifying aggregate = c('test', 'varSetName'), P-values across statistical tests and
varSets will be combined using ACAT.
A list with multiple items can be specified, in which case step-wise ACAT is applied.
For example, by specifying aggregate = list('test', 'varSetName'),
first P-values across statistical tests are combined, and the resulting P-values are
then combined across varSets.}

\item{group}{Variables to group by.
For example, if \code{group = c('unit', 'varSetName')} and \code{aggregate = 'test'},
for each unit-varSetName combination, P-values across statistical tests are combined.
Defaults to c("unit", "cohort", "varSetName","name", "unit", "pheno", "covar", "geneticModel", "MAFweight", "test"),
i.e. all grouping variables in an rvbResult object.
The variable(s) specified in \code{aggregate}, are excluded from the grouping variables.}

\item{fixpval}{Should P-values that are exactly zero or one be fixed? Defaults to \code{TRUE}.
The method used for fixing the P-values can be specified using the \code{fixpval_method} parameter.}

\item{fixpval_method}{Method used to fix p-value (if \code{fixpval = TRUE}). Methods include:
'minmax' = P-values that are exactly 1 are replaced by the maximum value below 1 present in the results;
P-values that are exactly 0 are replaced by the minimum value above 0 in the results.
'manual' = Specify the replacement P-values using \code{fixpval_maxP} and \code{fixpval_minP}.
'Liu' = Method recommended by the author of the ACAT R package
(see: https://github.com/yaowuliu/ACAT),
P-values of 1 are replaced by 1-1/d, where d is the number of p-values combined by ACAT.
Since no recommendation for P-values of 0 is given, these are replaced with the value specified using \code{fixpval_minP}.}

\item{fixpval_maxP}{Replace P-values that are exactly 1 with this P-value if \code{fixpval_method = 'manual'}.}

\item{fixpval_minP}{Replace P-values that are exactly 0 with this P-value if \code{fixpval_method = 'manual'}
or \code{fixpval_method = 'Liu'}.}

\item{warning}{Show warnings? Defaults to \code{TRUE}.}
}
\description{
Combine P-values in an \code{\link{rvbResult}} object using the ACAT method.
For details on the ACAT method see: (Liu \emph{et al.}, 2019).
}
\references{
Liu, Y. et al. ACAT: A Fast and Powerful p Value Combination Method for Rare-Variant Analysis in Sequencing Studies. The American Journal of Human Genetics 104, 410–421 (2019).
}
