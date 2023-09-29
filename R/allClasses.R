#===============================================================================
# All class definitions for the rvat package
#===============================================================================

#' The genoMatrix class for storing genotype data and associated sample and variant info.
#' 
#' @name genoMatrix
#' @rdname genoMatrix
#' @order 1
#' 
#' @description
#' The genoMatrix class is specifically designed to represent genotype data and associated sample and variant info.
#' It extends the Bioconductor \code{\link[SummarizedExperiment]{SummarizedExperiment}} class to accommodate genotype data.
#' 
#' @import methods
#' @import S4Vectors
#' @import SummarizedExperiment
#'
#' @usage NULL
#' @section Construct a genoMatrix:
#' `getGT()`: Construct a genoMatrix, see [`getGT()`] for details.
#'
#' @section Accessors:
#' The accessors are fully described in [`SummarizedExperiment::SummarizedExperiment-class`], and
#' include `assays(x)`, `rowData(x)`, `colData(x)`, `metadata(x)`, `dim(x)`, `nrow(x)`, `ncol(x)`, `colnames(x)` and `rownames(x)`. 
#' 
#' @section Getters:
#' In the following code snippets, x is a genoMatrix object
#' * `getAF(x)`: Returns alternate allele frequencies. These will be equal to the minor allele frequences after applying `flipToMinor()`.
#' * `getMAF(x)`: Returns variant minor allele frequencies.
#' * `getAC(x)`: Returns alternate allele counts.  These will be equal to the minor allele counts after applying `flipToMinor()`.
#' * `getMAC(x)`: Returns minor allele counts.
#' * `getCR(x, var=TRUE)`: Returns call-rates for variants or samples. 
#' `var=TRUE` by default, set to `FALSE` to return sample call-rates.
#' * `getNCarriers(x)`: Returns number of carriers of the alternate allele for each variant. 
#' Note that when geneticModel = 'recessive' it will return the number of homozygous carriers.
#' * `summariseGeno(x)`: Returns a per variant summary of genotype counts and hwe testing.
#' * `getCarriers(x, VAR_id=NULL,colDataFields=NULL,rowDataFields=NULL,groupBy=NULL)`: 
#' Return sample IDs for carriers of each of the variants in the genoMatrix. 
#' `VAR_id` can be specified to return for the specified subset of variants. 
#' `colDataFields` and `rowDataFields` can be specified include additional fields from `colData(x)` or `rowData(x)` in the output.
#' The `groupBy` parameter can be set to calculate carrier frequency among groups (such as cohort or phenotype), multiple groups can be set.
#' The `aggregate` parameter can be set to TRUE to return mean burden scores among the groupings.
#' 
#' note: Variant ploidy (diploid, XnonPAR, YnonPAR) are handled according to sample sex. 
#' Samples for which sex are not provided are excluded during AF calculation at non-diploid variants.
#' 
#' @section Subsetting:
#' 
#' Subsetting a genoMatrix object is equivalent to subsetting a SummarizedExperiment object,
#' and is fully described in fully described in [`SummarizedExperiment::SummarizedExperiment-class`].
#' Please see the example section for examples.
#' 
#' @section Recode:
#' In the following code snippets, x is a genoMatrix object
#' * `flipToMinor(x)`: Function to flip genotype dosages such that all GT values represent minor allele counts.
#' * `recode(x, geneticModel,imputeMethod,weights,MAFweights)`: Returns a recoded genoMatrix object, genetic model, imputation, and weights can be recoded. See [`recode()`] for details.
#' * `updateGT(x, SM = NULL, anno = NULL)`: Safe replacement of the colData or rowData table within a genoMatrix with a new table.
#' 
#' @section Rare variant testing:
#'  In the following code snippets, x is a genoMatrix object
#' * `aggregate()`: Returns per sample aggregate dosage for a genoMatrix object. 
#' The genoMatrix shouldn't contain missing values, use [`recode()`] to impute missing values.
#' Aggregation is dependent on `geneticModel`, `MAFweights` and `weights`, all can be set using the [`recode()`] method.
#' By default, an updated genoMatrix is returned with an `aggregate` field in `colData`. Set `returnGT` to `FALSE` to return a vector of aggregates.
#' * `assocTest()`: Perform aggregate (burden) and single variant assocation test. See [`assocTest()`] for details.
#' 
#' @seealso \code{\link{assocTest}}
#' @seealso \code{\link{recode}}
NULL


#' @export
setClass("genoMatrix", contains="SummarizedExperiment")


#' Connect and interact with a gdb
#' 
#' @name gdb
#' @rdname gdb
#' @order 1
#' 
#' @description
#' Compressed SQLite (".gdb") representation of genotype data and associated variant annotation and sample info tables
#' These allow for rapid and memory-efficient of loading sample genotype data en associated metadata within R.
#' The slots of the gdb class are inherited entirely from the RSQLite SQLiteConnection class.
#' A host of RVAT methods described here allow for convenient querying and manipulation of a gdb, for complex queries
#' users can also directly perform SQL queries on the gdb as exemplified in the examples and tutorials.
#' 
#' @usage NULL

#' @section Connect:
#' * `gdb(x)`: Connect to a  gdb object  (created using [`buildGdb()`]), where `x` is the filepath.

#' @section Build:
#' In the following code snippets, x is a gdb object.
#' * `buildGdb(x)`: build a gdb directly from a vcf, see [`buildGdb()`] for details.
#'
#' @section Getters:
#' In the following code snippets, x is a gdb object.
#' * `getGT(x, ...)`: Retrieve genotype data from the gdb, returns a [`genoMatrix()`] object. 
#' See [`getGT()`] for details.
#' * `getAnno(x, ...)`: Get an annotation table from gdb. See the [`getAnno()`] documentation for details.
#' * `getCohort(x, ...)`: Get a cohort table from gdb. See the [`getCohort()`] documentation for details.
#' * `listCohort(x)`: List sample cohort tables that have been uploaded to the gdb.
#' * `listAnno(x)`: List variant info tables that have been uploaded to the gdb.
#' 
#' @section Upload and delete sample and variant info tables:
#' In the following code snippets, x is a gdb object.
#' * `uploadCohort(x, ...)`: Upload cohort data tables to gdb. 
#' See the [`uploadCohort()`] documentation for details.
#' * `uploadAnno(x, ...)`: Upload variant annotation data into gdb. 
#' See the [`uploadAnno()`] documentation for details
#' * `mapVariants(x, ...)`: Map variants in the gdb onto features provided in a bed,gff/gtf or ranges file.
#' See [`mapVariants()`] for details.
#' * `dropTable(x, name)`: Drop table with the name specified in `name` from gdb and clear from annotation / cohort metadata tables.
#' 
#' @section Subsetting & Merging:
#' * `subSetGdb`: Create a new gdb that is a subset of the input gdb. See [`subsetGdb()`] for details.
#' * `concatGdb`: Function to concatenate gdb databases. See [`concatGdb()`] for details.
#' 
#' @section Exports:
#' * `writeVcf`: Convert a gdb to a vcf-file. See [`writeVcf()`] for details.
#' 
#' @seealso \code{\link{getGT}}
#' @seealso \code{\link{buildGdb}}
#' @seealso \code{\link{getAnno}}
#' @seealso \code{\link{getCohort}}
#' @seealso \code{\link{uploadCohort}}
#' @seealso \code{\link{uploadAnno}}
#' @seealso \code{\link{mapVariants}}
#' @seealso \code{\link{subsetGdb}}
#' @seealso \code{\link{concatGdb}}
#' @seealso  \code{\link{writeVcf}}
#' @slot path Path to gdb object
#' @keywords gdb
NULL


#' gdb-class
#'
#' @rdname gdb
#' @import DBI
#' @importClassesFrom RSQLite SQLiteConnection
#' @export
setClass("gdb", contains = "SQLiteConnection")


#' Class to manage multiple varSets
#' @name varSetList
#' @rdname varSetList
#' @order 1
#' @usage NULL
#' @description
#' An S4 class to manage varSets. varSets can be generated from annotations using
#' the [`buildVarSet()`] method. Specific units and/or annotations can be retrieved using
#' the `getVarSet` method. A varSetList can be used directly as input in [`assocTest()`] to perform
#' burden/single variant association tests on the varSets included in the varSetList.
#' The on-disk equivalent of a `varSetList` is [`varSetFile()`], from which varSets can be retrieved 
#' in the same way as a varset list using the `getVarSet` method.
#'
#' @section Build a varSetList:
#' * `buildVarSet(x, ...)`: Generate a varSetList or [`varSetFile()`] that stores weighted variant sets
#' for use in association testing. This can be based on 1) annotations uploaded to the gdb or 2) a data.frame including annotations.
#' See [`buildVarSet()`] for details.
#'
#' @section Getters:
#' In the following code snippets, x is a varSetList object.
#' * `getVarSet(x, unit, varSetName)`: Retrieve varSets for specified units and/or varSetNames.
#' * `listUnits(x)`: Return a vector of all units included in the varSetList.
#' * `listVarSets(x)`: Return a vector of all varSetNames included in the varSetList.
#' 
#' @section Subsetting:
#' A varSetList can be subsetted in the same way as a normal R list,
#' however, `getVarSet` is the most convenient way to select varSets (see above).
#' Some examples, where x is a varSetList object:
#' * `x[[i]]`: Return the i'th varSet
#' * `x[i:j]`: Return the i'th till j'th varSets
#' 
#' @section Association testing:
#' A varSetList can be directly supplied to the [`assocTest()`] gdb method, using the 
#' `varSet` parameter. Association tests will then be performed for each varSet included
#' in the varSetList.
#' 
#' @section Miscellaneous:
#' In the following code snippets, x is a varSetList object:
#' * `collapseVarSetList(x, unit = "unnamed", varSetName = "unnamed)`: 
#' Merge all varSets into one varSet. All weights will be set to 1. 
#' Optionally, the unit name and varSetName can be specified.
#' This method is mainly used when you want to load the genotypes for all variants
#' in a varSetList using [`getGT()`]. 
#' * `write(x, file = "data", append = "FALSE")`: Write the varSetList to disk, in 
#' the [`varSetFile`] format.
#' 
#' 
#' @seealso \code{\link{varSetFile}}
#' @seealso \code{\link{buildVarSet}}
#' @seealso \code{\link{getGT}}
#' @seealso \code{\link{assocTest}}
#' @keywords varSetList
NULL

#' Class that represents a set of variants and weights
#'
#' An S4 class to manage an individual varSet record. Usually multiple varSets
#' will be included in a [`varSetList`] or [`varSetFile`] object, which can be generated using the [`buildVarSet`] method.
#' @name varSet
#' @rdname varSet
#' @order 1
#' @slot unit Unit name (such as the gene or transcript ID)
#' @slot varSetName Name of the variant set (such as 'LOF', 'moderate', 'CADD')
#' @slot VAR_id VAR_ids included in the varSet
#' @slot w weights included in the varSet
#' 
#' @section Getters:
#' In the following code snippets, x is a varSet object.
#' * `listVars(x)`: Return a vector of all VAR_ids included in the varSet
#' * `listWeights(x)`: Return a vector of all weights in the varSet
#' 
#' @section Retrieve genotypes:
#' A varSet object can be supplied to the `varSet` parameter in the [`getGT()`] method to load the 
#' variants included in the varSet. 
#' 
#' @seealso \code{\link{varSetFile}}
#' @seealso \code{\link{varSetList}}
#' @seealso \code{\link{getGT}}
#' @seealso \code{\link{buildVarSet}}
#' 
#' @export
setClass("varSet",
         representation(
           unit="character",
           varSetName="character",
           VAR_id="character",
           w="character"
         ))

#' @rdname varSetList
#' @usage NULL
#' @export
setClass("varSetList", 
  representation(
    varSets="list",
    units="character"
  ))


#' Class to manage interactions with a varSetFile
#' @name varSetFile
#' @rdname varSetFile
#' @order 1
#' @usage NULL
#' @description
#' An S4 class to manage interactions with varSetFiles. varSets can be generated from annotations using
#' the [`buildVarSet()`] method. Specific units and/or annotations can be loaded using
#' the `getVarSet` method. A varSetFile can be used as input in [`assocTest()`] to perform
#' burden/single variant association tests on the varSets included in the varSetFile.
#'
#' @section Build a varSetFile:
#' * `buildVarSet(x, ...)`: Generate a varSetList or [`varSetFile()`] that stores weighted variant sets
#' for use in association testing. This can be based on 1) annotations uploaded to the gdb or 2) a data.frame including annotations.
#' See [`buildVarSet()`] for details.
#' 
#' @section Connect to a varSetFile:
#' * `varSetFile(path, memlimit = 5000)`: Connect to a varSetFile object. 
#'
#' @section Getters:
#' In the following code snippets, x is a varSetFile object.
#' * `getVarSet(x, unit, varSetName)`: Retrieve varSets for specified units and/or varSetNames.
#' Output will be of class [`varSetList`].
#' * `listUnits(x)`: Return a vector of all units included in the varSetFile
#' * `listVarSets(x)`: Return a vector of all varSetNames included in the varSetFile
#' 
#' @section Association testing:
#' A varSetFile can be directly supplied to the [`assocTest()`] gdb method, using the 
#' `varSet` parameter. Association tests will then be performed for each varSet included
#' in the varSetFile.
#' 
#' 
#' @seealso \code{\link{varSetList}}
#' @seealso \code{\link{buildVarSet}}
#' @seealso \code{\link{assocTest}}
#' @keywords varSetFile
NULL


#' varSetFile-class
#' @rdname varSetFile
#' @usage NULL
#' @export
setClass("varSetFile",
         representation(
           path="character",
           units="character")
)


# VIRTUAL class rvatResult

#' Class for handling and visualizing results generated using rvat.
#' 
#' @name rvatResult
#' @rdname rvatResult
#' @order 1
#' 
#' @description
#' The rvatResult class is specifically designed to represent association results
#' generated using RVAT (with [`assocTest`] or [`geneSetAssoc`] for example).
#' It extends the BioConductor \code{\link[S4Vectors]{DFrame}} class, and allows for basic
#' operations such as subsetting and merging as well as visualization (manhattan, qqplot, forest plots) and 
#' downstream analyses (e.g. [`ACAT`]). 
#' Different type of results have their own subclasses (rvbResult, singlevarResult, gsaResult) that inherit from
#' rvatResult.
#' 
#' @importClassesFrom S4Vectors DataFrame DFrame
#'
#' @usage NULL
#' @section Constructing:
#' * `rvbResult(object, header = TRUE)`: Intitialize an `rvbResult` object, 	
#' including singlevar results as generated by the [`assocTest`] method. 
#' `object` can be either a data.frame/DataFrame or a filepath pointing to the results.
#' * `singlevarResult(object, header = TRUE)`: Intitialize an `singlevarResult` object, 	
#' including rvb results as generated by the [`assocTest`] method. 
#' `object` can be either a data.frame/DataFrame or a filepath pointing to the results.
#' * `gsaResult(object, header = TRUE)`: Intitialize an `singlevarResult` object, 	
#' including rvb results as generated by the [`geneSetAssoc`] or [`assocTest`] methods. 
#' `object` can be either a data.frame/DataFrame or a filepath pointing to the results
#' * `readResults(path, header = TRUE, type = NULL, sep = "\t")`: Alternatively,
#' results can be read using `readResults` where `type` can be set to `rvbResult`/`singlevarResult`/`gsaResult`.
#' If type is not set, it will be inferred from the input.
#' 
#'
#' @section Subsetting:
#' 
#' Subsetting an rvatResult object is equivalent to subsetting a DataFrame/DFrame object and is fully described
#' in [`S4Vectors::DataFrame-class`]. DataFrame objects behave very 
#' similar to the base `data.frame`, with the most notable exceptions that 
#' row names are optional and it can hold alternative vectors such as run-length encoded  vectors ([`S4Vectors::Rle-class`]).
#' Please see the example section for examples.
#'
#' @section Merging:
#'
#' * `merge(x, y, by)`: Merge an `rvatResult` object with a `data.frame` or `DataFrame` object. 
#' The `by` parameter specifies which variables to join by. 
#' For example, by = c("a" = "b") will match x$a to y$b. 
#' Join by multiple variable by providing a vector of length>1 to `by`.
#' 
#' @section Displaying:
#'In the following code snippets, x is an rvatResult object
#' * `show(x)`: By the default the number of rows and columns will be displayed + the first five rows.
#' * `summary(x, asList = FALSE)`: Shows a summary of the contents of the `rvatResult`.
#'  Optionally the output can be stored as a list by setting `asList = TRUE`.
#' * `topResult(x, n=10)`: Show the top N most significant results (based on the P-value).
#' 
#' @section Visualization:
#' * `qqplot`: Plot a qqplot based on the P-values stored in an rvatResult object. See [`qqplot`] for details
#' * `manhattan`: Generate a manhattan plot for an `rvatResult` object. See [`manhattan`] for details.
#' * `rvatViewer`: Interactive visualization, see [`rvatViewer`] for details.
#'
#' @section Downsteam analyses:
#' * `geneSetAssoc`: An `rvbResult` can be used as input for `geneSetAssoc` to perform 
#' gene set analyses. See [`geneSetAssoc`] for details.
#' * `ACAT`: P-values in an `rvbResult` can be combined using the ACAT method. See [`ACAT`] for details.
#' 
#' @section Writing:
#' * `writeResult(x, ...)`: Write an `rvatResult` to disk, see [`writeResult`] for details.
#'
#' @seealso \code{\link{assocTest}}
#' @seealso \code{\link{geneSetAssoc}}
#' @seealso \code{\link{qqplot}}
#' @seealso \code{\link{manhattan}}
#' 
NULL

#' rvatResult-class
#' @rdname rvatResult
#' @usage NULL
#' @export
setClass("rvatResult",
         representation(
           "VIRTUAL"),
         contains="DFrame")

#' rvbResult-class
#' @rdname rvatResult
#' @usage NULL
#' @export
setClass("rvbResult",
         contains="rvatResult")

#' singlevarResult-class
#' @rdname rvatResult
#' @usage NULL
#' @export
setClass("singlevarResult",
         contains="rvatResult")

#' gsaResult-class
#' @rdname rvatResult
#' @usage NULL
#' @export
setClass("gsaResult",
         contains="rvatResult")


#' Generate and manage permutations
#' 
#' @name resamplingFile
#' @rdname resamplingFile
#' @order 1
#' 
#' @description
#' This class is designed to generate and manage permutations. For some purposes 
#' it is useful to store permutations, for example for empirically establishing the correlation
#' among genes using permutations. In this case identical permutations should be used across genes,
#' which can be achieved by storing the permutations in a resamplingFile.
#' 

#' @usage NULL
#' @section Construct a  resamplingFile:
#' `buildResamplingFile`: Construct a resamplingfile, see [`buildResamplingFile`] for details.
#'
#' @section Connect:
#' * `resamplingFile(x)`: Connect to a resamplingFile object (created using [`buildResamplingFile()`]), where `x` is the filepath.

#' @section Association testing:
#' A resamplingFile can be directly provided to [`assocTest()`], to run resampled association tests.
NULL


#' resamplingFile-class
#' @rdname resamplingFile
#' @usage NULL
#' @export
setClass("resamplingFile",
         representation(
           path="character",
           methodResampling="character",
           nSamples="numeric",
           nResampling="numeric"
         ))

#' Class that represents a geneSet
#'
#' An S4 class to manage an individual geneSet record. Usually multiple geneSets
#' will be included in a [`geneSetList`] or [`geneSetFile`] object, which can be generated using the [`buildGeneSet`] method.
#' @name geneSet
#' @slot geneSetName name of the geneset
#' @slot units units (usually genes) in the geneSet
#' @slot w optional weights (defaults to 1 for each unit)
#' @slot metadata optional metadata for each geneset
#' 
#' @section Getters:
#' In the following code snippets, x is a geneSet object.
#' * `listUnits(x)`: Return a vector of all units included in the geneSet
#' * `listWeights(x)`: Return a vector of all weights in the geneSet
#' 
#' @seealso \code{\link{geneSetFile}}
#' @seealso \code{\link{geneSetList}}
#' @seealso \code{\link{buildGeneSet}}
#' @export
setClass("geneSet", 
         representation(
           geneSetName="character",
           units="character",
           w="character",
           metadata="character"
         ))

#' Class to manage multiple geneSets
#' @name geneSetList
#' @rdname geneSetList
#' @order 1
#' @usage NULL
#' @description
#' An S4 class to manage geneSets. A geneSetList can be generated from gmt-files or a list of gene sets using
#' the [`buildGeneSet()`] method. Specific geneSets can be retrieved using
#' the [`getGeneSet`] method. A geneSetList can be used directly as input in [`geneSetAssoc()`] or [`assocTest`] to perform gene set analyses.
#' The on-disk equivalent of a `geneSetList` is [`geneSetFile()`], from which geneSets can be retrieved 
#' in the same way as a geneSetList using the `getGeneSet` method.
#'
#' @section Build a geneSetList:
#' * `buildGeneSet(x, ...)`: Generate a geneSetList or [`geneSetFile()`] that stores (weighted) gene sets
#' for use geneSetAnalyses. This can be based on 1) GMT-files (https://www.gsea-msigdb.org/gsea/msigdb/) and 2) a list of vectors of units per gene set.
#' See [`buildGeneSet()`] for details.
#'
#' @section Getters:
#' In the following code snippets, x is a geneSetList object.
#' * `geneGeneSet(x, set)`: Retrieve geneSets for specified names.
#' * `listGeneSets(x)`: Return a vector of all gene sets included in the geneSetList
#' * `listUnits(x)`: Return a vector of all units included across the gene sets in the geneSetList
#' 
#' @section Subsetting:
#' A geneSetList can be subsetted in the same way as a normal R list,
#' however, [`getGeneSet`] is the most convenient way to select genesets (see above).
#' Some examples, where x is a geneSetList object:
#' * `x[[i]]`: Return the i'th geneSet
#' * `x[i:j]`: Return the i'th till j'th geneSets.
#' 
#' @section Gene set analyses:
#' A geneSetList can be directly supplied to the [`geneSetAssoc()`] method, using the 
#' `geneSet` parameter, in combination with an [`rvbResult`] object.
#' To perform gene set burden analyses, [`assocTest`] can be used.
#' 
#' @section Miscellaneous:
#' In the following code snippets, x is a geneSetList object:
#' * `write(x, file = "data", append = "FALSE")`: Write the geneSet to disk, in 
#' the [`geneSetFile`] format.
#' * `remapIDs(x,...)`: Remap IDs used in geneSets, see [`remapIDs`] for details.
#' * `dropUnits(x, unit = NULL)`: Remove specified units from geneSets included in the geneSetList.
#' 
#' 
#' @seealso \code{\link{geneSetFile}}
#' @seealso \code{\link{geneSet}}
#' @seealso \code{\link{geneSetAssoc}}
#' @seealso \code{\link{assocTest-aggregateFile}}
#' @keywords geneSetList
NULL

#' @rdname geneSetList
#' @usage NULL
#' @export
setClass("geneSetList", 
         representation(
           geneSets="list",
           geneSetNames="character"
         ))

#' Class to manage interactions with a geneSetFile
#' @name geneSetFile
#' @rdname geneSetFile
#' @order 1
#' @usage NULL
#' @description
#' An S4 class to manage interactions with geneSetFiles. A geneSetFile can be generated from gmt-files or a list of gene
#' sets using the [`buildGeneSet()`] method. Specific geneSets can be retrieved using
#' the [`getGeneSet`] method. A geneSetFile can be used as input in [`geneSetAssoc()`] or [`assocTest()`] 
#' to perform competitive/self-contained gene set analyses.
#'
#' @section Build a geneSetFile:
#' * `buildGeneSet(x, ...)`: Generate a [`geneSetList`] or geneSetFile that stores (weighted) gene sets
#' for use geneSetAnalyses. This can be based on 1) GMT-files (https://www.gsea-msigdb.org/gsea/msigdb/) and 2) a list of #' vectors of units per gene set.
#' See [`buildGeneSet()`] for details.
#' 
#' @section Connect to a geneSetFile:
#' * `geneSetFile(path)`: Connect to a geneSetFile object. 
#'
#' @section Getters:
#' In the following code snippets, x is a geneSetFile object.
#' * `getGeneSet(x, geneSet)`: Retrieve geneSets for specified units and/or geneSetNames.
#' Output will be of class [`geneSetList`].
#' * `listUnits(x)`: Return a vector of all units included in the geneSetFile
#' * `listGeneSets(x)`: Return a vector of all geneSetNames included in the geneSetFile
#' 
#' @section Association testing:
#' A geneSetFile can be directly supplied to the [`geneSetAssoc()`] method, using the 
#' `geneSet` parameter, in combination with an [`rvbResult`] object.
#' To perform gene set burden analyses, [`assocTest`] can be used.
#' 
#' 
#' @seealso \code{\link{geneSetList}}
#' @seealso \code{\link{buildGeneSet}}
#' @seealso \code{\link{geneSetAssoc}}
#' @seealso \code{\link{assocTest-aggregateFile}}
#' @keywords geneSetFile
NULL

#' geneSetFile-class
#' @rdname geneSetFile
#' @usage NULL
#' @export
setClass("geneSetFile",
         representation = representation(
           path="character",
           sets="character"),
         validity = function(object) {
           msg=NULL
           
           # Ensure varSet units have been sorted
           if (mean(object@sets==sort(object@sets))<1)
           {
             msg=c(msg,"Invalid geneSetFile. Records are not correctly sorted by the 'set' field. 
                   Resort geneSet file by set, or else regenerate the GeneSetFile using buildGeneSetFile")
           }
           
           # Return
           if (length(msg)){msg} else{TRUE}
         }
)

#' Class to manage interactions with an aggregateFile
#' @name aggregateFile
#' @rdname aggregateFile
#' @order 1
#' @usage NULL
#' @description
#' The [`aggregate`] method saves genotypes aggregates compressed on disk. The 
#' aggregateFile class manages connecting/interacting with these files and retrieving
#' units of interest.
#'
#' @section Build an aggregateFile:
#' * `aggregate(x, ...)`: Returns an aggregate of genotypes for each individual.
#' See [`aggregate()`] for details.
#' 
#' @section Connect to an aggregateFile:
#' * `aggregateFile(path)`: Connect to an aggregateFile object. 
#'
#' @section Getters:
#' In the following code snippets, x is an aggregateFile object.
#' * `getUnit(x, unit)`: Retrieve aggregates for specified unit(s). 
#' Use `listUnits(x)` to list the units includes in the aggregateFile
#' Output will be a matrix.
#' * `listUnits(x)`: Return a vector of all units included in the aggregateFile
#' * `listSamples(x)`: Return a vector of all sample IDs included in the aggregateFile
#' 
#' @section Association testing:
#' An aggregateFile can be directly supplied to the [`assocTest()`] method.
#' 
#' @section Merging:
#' Aggregate files can be merged using the [`mergeAggregateFiles`] method. 
#' 
#'
#' @seealso \code{\link{mergeAggregateFiles}}
#' @seealso \code{\link{aggregateFileList}}
#' @seealso \code{\link{assocTest-aggregateFile}}
#' @seealso \code{\link{aggregate}}
NULL

#' @name aggregateFile
#' @usage NULL
#' @export
setClass("aggregateFile",
         representation(
           path="character",
           units="character",
           samples="character"))



#' Class to facilitate merging aggregateFiles
#' 
#' @name aggregateFileList
#' @rdname aggregateFileList
#' @order 1
#' @usage NULL
#' @description
#' 
#' Class to facilitate merging aggregateFiles. By providing a vector of aggregateFile filepaths,
#' aggregateFileList will check whether identical samples are included and if duplicated units are included.
#' The aggregateFileList can then be used to merge the aggregateFiles into either a new aggregateFile,
#' or merge all included aggregates into a single aggregate score per sample ([`mergeAggregateFiles`]).
#'
#' @section Initialize an aggregateFileList object:
#' * `aggregateFileList(filelist, checkDups=TRUE)`: Here `filelist` is a vector of [`aggregateFile`]
#' filepaths. `checkDups` is set to `TRUE` by default, in which case an error raised if
#' unit names are duplicated across aggregateFiles.
#'
#' @section Getters:
#' In the following code snippets, x is an aggregateFileList object.
#' * `listUnits(x)`: Return a vector of all units included across aggregateFiles in the aggregateFileList.
#' * `listSamples(x)`: Return a vector of all samples included across aggregateFiles in the aggregateFileList.
#' 
#' @section Merge aggregateFiles:
#' * `mergeAggregateFiles(object,aggregate = TRUE,output = NULL,verbose = TRUE)`: 
#' Merge the aggregatefiles included in the aggregateFileList. 
#' Either aggregate all aggregates into one single aggregate per sample, or merge
#' all aggregatefiles into a new aggregateFile. See [`mergeAggregateFiles`] for details.
#' 
#' @seealso \code{\link{mergeAggregateFiles}}
#' @seealso \code{\link{assocTest-aggregateFile}}
#' @seealso \code{\link{aggregate}}
NULL

#' @name aggregateFileList
#' @usage NULL
#' @export
setClass("aggregateFileList",
         representation(
           paths="character",
           units="character",
           samples="character"))



#' An S4 class to store and handle null models for mixed linear model gene set analysis.
#' 
#' @name nullModelGSA-class
#' @rdname nullModelGSA-class
#' @order 1
#' 
#' @description
#' An S4 class to store and handle null models for mixed linear model gene set analysis.
#' 
#' @usage NULL 
#' @section Generate a GSA null model:
#' `fitNullModelGSA`: Generate a GSA null model, see [`buildCorMatrix()`] for details.
#' 
#' @section Getters:
#' In the following code snippets, x is a nullModelGSA object
#' * `getCovar(x)`: Returns fixed effects 
#' * `getResults(x)`: Returns input results
#' * `getNullModel(x)`: Returns the null model
#' * `listUnits(x)`: Returns the units included
#' 
#' @seealso \code{\link{geneSetAssoc}}
#' @seealso \code{\link{buildCorMatrix}}
#' @seealso \code{\link{assocTest}}
NULL

#' nullModelGSA-class
#' 
#' An S4 class to store and handle null models for gene set analysis
#' @name nullModelGSA-class
#' @slot nullmodel the nullmodel
#' @slot results the rvbResults the nullmodel was based on
#' @slot units units 
#' @slot covar covariates included in the null model
#' @slot method which method was used, currently 'GENESIS' is implemented.
#' @export
setClass("nullModelGSA", 
         representation(
           nullmodel="ANY",
           ID="integer",
           results="rvatResult",
           units="character",
           covar="character",
           method="character"
         ))

#' @rdname nullModelGSA-class
#' @usage NULL
#' @param x \code{\link{nullModelGSA-class}} object
#' @export
setMethod("listUnits", signature = "nullModelGSA",
          definition = function(object) {
            object@units
          })

#' @rdname nullModelGSA-class
#' @usage NULL
#' @export
setMethod("length", signature = "nullModelGSA",
          definition = function(x) {
            length(listUnits(x))
          })


#' @rdname nullModelGSA-class
#' @usage NULL
#' @export
setMethod("show", 
          signature = "nullModelGSA",
          definition = function(object) {
            message(sprintf("nullModelGSA, generated with %s", object@method))
            message(sprintf("Contains a null model including %s units",length(object)))
            message(sprintf("Covar: %s", paste(getCovar(object), collapse=",")))
          })

#' @rdname nullModelGSA-class
#' @usage NULL
#' @export
setMethod("getNullModel", signature = "nullModelGSA",
          definition = function(object) {
            object@nullmodel
          })


#' @rdname nullModelGSA-class
#' @usage NULL
#' @export
setMethod("getResults", signature = "nullModelGSA",
          definition = function(object) {
            object@results
          })

#' @rdname nullModelGSA-class
#' @usage NULL
#' @export
setMethod("getCovar", signature = "nullModelGSA",
          definition = function(object) {
            object@covar
          })