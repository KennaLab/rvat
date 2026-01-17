#' library(rvatData)
#' 
#' # Build a varSetFile including variants with a moderate predicted impact
#' gdb <- create_example_gdb()
#' varsetfile_moderate <- tempfile()
#' buildVarSet(object = gdb, 
#'             output = varsetfile_moderate,
#'             varSetName = "Moderate", 
#'             unitTable = "varInfo", 
#'             unitName = "gene_name",
#'             where = "ModerateImpact = 1")
#' 
#' # connect to the varSetFile
#' varsetfile <- varSetFile(varsetfile_moderate)
#' 
#' # list included units and varSets 
#' units <- listUnits(varsetfile)
#' head(units)
#' varsets <- listVarSets(varsetfile)
#' head(varsets)
#' 
#' # basic operations
#' length(varsetfile)
#' 
#' # metadata
#' metadata(varsetfile)
#' getRvatVersion(varsetfile)
#' getGdbId(varsetfile)
#' 
#' # retrieve varSets 
#' varsets <- getVarSet(varsetfile, unit = c("SOD1", "FUS"))
#' # this returns a varSetList, which is an in-memory representation of varSets
#' # most of the methods that work on a varSetFile also work on a varSetList (see ?varSetLit for details)
#' getVarSet(varsets, unit = "SOD1")
#' 
#' # see e.g., ?assocTest and ?aggregate for downstream methods that can loop through varsetfiles and varsetlists.
