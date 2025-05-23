% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allGenerics.R
\name{buildVarSet}
\alias{buildVarSet}
\title{buildVarSet}
\description{
Generate optionally weighted variant sets using annotation table(s).
See the tutorials for examples.
For building varSets directly from the \code{\link{gdb}}: see \code{\link{buildVarSet-gdb}} for details\cr
For building varSets interactively from a data.frame  seee \code{\link{buildVarSet-data.frame}} for details\cr
}
\examples{

library(rvatData)

# Build a varSetFile including variants with a moderate predicted impact
gdb <- create_example_gdb()
varsetfile_moderate <- tempfile()
buildVarSet(object = gdb, 
            output = varsetfile_moderate,
            varSetName = "Moderate", 
            unitTable = "varInfo", 
            unitName = "gene_name",
            where = "ModerateImpact = 1")

# Build a varSetFile that contains CADD scores
varsetfile_cadd <- tempfile()
buildVarSet(object = gdb, 
            output = varsetfile_cadd,
            varSetName = "CADD", 
            unitTable = "varInfo", 
            unitName = "gene_name",
            weightName = "CADDphred")

# in addition to building a varSetFile from a gdb,
# it can also be build directly from a data.frame
anno <- getAnno(gdb, "varinfo", where = "gene_name in ('SOD1', 'FUS')")
varsetfile_from_df <- tempfile()
buildVarSet(
  anno,
  unitName = "gene_name",
  fields = c("HighImpact"),
  output = varsetfile_from_df
)

# connect to varsetfile and retrieve variant sets
varsetfile <- varSetFile(varsetfile_moderate)
varsets <- getVarSet(varsetfile, unit = c("SOD1", "FUS"))

# see ?getVarSet, ?varSetFile and ?varSetList for more details on connecting and handling varsetfiles.
# see e.g., ?assocTest and ?aggregate for downstream methods that can loop through varsetfiles and varsetlists.


}
