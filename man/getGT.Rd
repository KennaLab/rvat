% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allGenerics.R, R/gdb.R
\name{getGT}
\alias{getGT}
\alias{getGT,gdb-method}
\title{Load genotypes from a gdb.}
\usage{
getGT(
  object,
  varSet = NULL,
  VAR_id = NULL,
  ranges = NULL,
  cohort = NULL,
  anno = NULL,
  annoFields = NULL,
  includeVarInfo = FALSE,
  checkPloidy = NULL,
  varSetName = "unnamed",
  unit = "unnamed",
  padding = 250,
  verbose = TRUE,
  strict = TRUE
)
}
\arguments{
\item{object}{an object of class \code{\link[=gdb-class]{gdb}}}

\item{varSet}{\code{\link{varSet}} object. If specified, the \code{VAR_id}, \code{varSetName} and \code{unit} parameters will be ignored.}

\item{VAR_id}{Character vector containing target VAR_id.}

\item{ranges}{Extract variants within specified ranges.
Ranges can be specified as a data.frame, including at least 'CHROM','start', and 'end' columns, or
can be a \code{\link[GenomicRanges:GRanges-class]{GenomicRanges::GRanges}} object.}

\item{cohort}{Optional use of cohort data previously uploaded to the gdb (see \code{uploadCohort}).
If a valid cohort name is provided, then the uploaded data for this cohort is used to filter and annotate the returned genoMatrix object.
If a dataframe is provided, then this is assumed to conform to the SM table constraints required for genoMatrix objects (see genoMatrix).}

\item{anno}{Optional use of variant annotation data previously uplodated to the gdb (see \code{uploadAnno}).
If a valid annotation table name is provided, the variant annotations will be included in the \code{rowData} of the \code{genoMatrix}.
Note that currently only annotation tables that include one row per variant can be included.
The \code{annoFields} parameter can be used to retain only specified fields from the annotation table.}

\item{annoFields}{A vector of field names to retain if the \code{anno} parameter is set.}

\item{includeVarInfo}{Include variant info ('var' table from the gdb) in the \code{genoMatrix}? Defaults to \code{FALSE}.
Note that setting this parameter to \code{TRUE} will override the \code{anno}/\code{annoFields} parameters.}

\item{checkPloidy}{Version of the human genome to use when assigning variant ploidy (diploid, XnonPAR, YnonPAR).
Accepted inputs are GRCh37, hg19, GRCh38, hg38.
If not specified, the genome build in the \code{\link{gdb}} will be used, if available (included in the \code{genomeBuild} parameter was set in \code{\link{buildGdb}}).
Otherwise, if the genome build is not included in the gdb metadata, and no value is provided, then all variants are assigned the default ploidy of "diploid"}

\item{varSetName}{Optional name for the set of variants, for example: 'missense or 'LOF' (ignored if \code{varSet} is specified.)}

\item{unit}{Optional 'unit' name, for example: 'SOD1' or 'ENSG00000142168' (ignored if \code{varSet} is specified.)}

\item{padding}{Number of basepairs to extend the search region beyond the specified genomic ranges to capture variants where the reference allele (REF) overlaps the input ranges,
but the POS of the variant falls outside the ranges. This accounts for variants where the REF allele spans multiple base pairs.}

\item{verbose}{Should the method be verbose? Defaults to \code{TRUE}.}

\item{strict}{Should strict checks be performed? Defaults to \code{TRUE}. Strict tests currently includes
checking whether supplied varSetFile/varSetList/varSet was generated from the same gdb as specified in \code{object}.}
}
\value{
A \code{\link{genoMatrix}} object.
}
\description{
Method to retrieve a \code{\link{genoMatrix}} for variants specified by \code{varSet} or a vector of VAR_ids,
and samples as specified by the \code{cohort} parameter.
The \code{checkPloidy} parameter can be set to \code{GRCh37} (or \code{hg19}) or \code{GRCh38} (or \code{hg38}) to
assign variant ploidy (diploid,XnonPAR,YnonPAR). If not specified, the genome build in the \code{\link{gdb}} will be used, if available (included in the \code{genomeBuild} parameter was set in \code{\link{buildGdb}}).
}
\examples{
library(rvatData)
gdb <- gdb(rvat_example("rvatData.gdb"))

# retrieve genotypes of a set of variants based on their VAR_ids
varinfo <- getAnno(gdb, table = "varinfo", where = "gene_name = 'SOD1' and ModerateImpact = 1")
GT <- getGT(
  gdb,
  VAR_id = varinfo$VAR_id,
  cohort = "pheno")

# retrieve genotypes of a set of variants in a varSet
varsetfile <- varSetFile(rvat_example("rvatData_varsetfile.txt.gz"))
varset <- getVarSet(varsetfile, unit = "NEK1", varSetName = "High")
GT <- getGT(
  gdb,
  varSet = varset,
  cohort = "pheno")
# see ?varSetFile and ?getVarSet for more details

# retrieve genotypes for a genomic interval
GT <- getGT(
  gdb,
  ranges = data.frame(CHROM = "chr21", start = 31659666, end = 31668931),
  cohort = "pheno")

# the `anno` parameter can be specified to include variant annotations in the rowData of the genoMatrix
GT <- getGT(
  gdb,
  ranges = data.frame(CHROM = "chr21", start = 31659666, end = 31668931),
  cohort = "pheno",
  anno = "varInfo",
  annoFields = c("VAR_id", "CHROM", "POS", "REF", "ALT", "HighImpact", "ModerateImpact", "Synonymous")
)
head(rowData(GT))

# the includeVarInfo parameter is a shorthand for include the "var" table
GT <- getGT(
  gdb,
  ranges = data.frame(CHROM = "chr21", start = 31659666, end = 31668931),
  cohort = "pheno",
  includeVarInfo = TRUE
)
head(rowData(GT))

# The `checkPloidy` parameter can be set to the version of the human genome to use
# to assign variant ploidy. (diploid, XnonPAR, YnonPAR). Accepted inputs are GRCh37, hg19, GRCh38, hg38.
# We recommend, however, to set the genome build when building the gdb: the genome build will theb
# be included in the gdb metadata and used automatically. see ?buildGdb for details
varinfo <- getAnno(gdb, table = "varinfo", where = "gene_name = 'UBQLN2' and ModerateImpact = 1")
GT <- getGT(
  gdb,
  VAR_id = varinfo$VAR_id,
  cohort = "pheno",
  includeVarInfo = TRUE,
  checkPloidy = "GRCh38"
)

# see ?genoMatrix for more details on the genoMatric class.

}
