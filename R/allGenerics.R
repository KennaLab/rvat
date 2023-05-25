#===============================================================================
# All generics with multiple methods
#===============================================================================

#' writeVcf
#'
#' Takes input and outputs to file in vcf format.
#'
#' @param object Input [`gdb`].
#' @param output Output file path.
#' @param VAR_id Variants to include in output.
#' @param IID Samples to include in output.
#' @param includeGeno Can be reset to false if individual genotypes are not required in output.
#' Defaults to `TRUE`.
setGeneric("writeVcf", function(object, output, VAR_id, IID, includeGeno=TRUE) standardGeneric("writeVcf"))

#' @rdname gdb
#' @usage NULL
#' @export
setGeneric("listAnno", function(object) standardGeneric("listAnno"))

#' @rdname gdb
#' @usage NULL
#' @export
setGeneric("listCohort", function(object) standardGeneric("listCohort"))

#' @rdname gdb
#' @usage NULL
#' @export
setGeneric("getRvatVersion", function(object) standardGeneric("getRvatVersion"))


#' Extract variants from a gdb based on a set of ranges
#' 
#' Extract variant IDs from a [`gdb`] based on a set of ranges.
#'
#' @param object a [`gdb`] object
#' @param ranges Can be a data.frame, including at least 'CHROM','start', and 'end' columns or
#' a [`GenomicRanges::GRanges`] object. 
#' @param return What to return. Defaults to `VAR_id`, which returns a vector of VAR_ids that map to specified ranges.
#' @export
setGeneric("extractRanges", function(object, ranges, return = "VAR_id") standardGeneric("extractRanges"))


#' Get an annotation table from a gdb
#'
#' Get an annotation table from a \code{\link[=gdb-class]{gdb}} object.
#'
#' @param object an object of class \code{\link[=gdb-class]{gdb}}
#' @param table base table to query
#' @param fields columns to retain
#' @param left left join operations to perform
#' @param inner inner join operations to perform
#' @param VAR_id retain only variants with matching ID
#' @param ranges Extract variants within specified ranges. 
#' Ranges can be specified as a data.frame, including at least 'CHROM','start', and 'end' columns, or
#' can be a [`GenomicRanges::GRanges`] object. 
#' @param where An SQL compliant where clause to filter output; eg: "CHROM=2 AND POS between 5000 AND 50000 AND AF<0.01 AND (cadd.caddPhred>15 OR snpEff.SIFT='D')".

#' @export
setGeneric("getAnno", function(object,table,fields="*",left=c(),inner=c(),VAR_id=c(),ranges=NULL,where=c()) standardGeneric("getAnno"))

#' Get an cohort table from a gdb
#'
#' Get a cohort table from a \code{\link[=gdb-class]{gdb}} object.
#'
#' @param object an object of class \code{\link[=gdb-class]{gdb}}
#' @param cohort name of cohort to get.
#' @param fields columns to retain
#' @export
setGeneric("getCohort", function(object, cohort, fields="*") standardGeneric("getCohort"))

#' Load genotypes from a gdb.
#'
#' Method to retrieve a [`genoMatrix`] for variants specified by `varSet` or a vector of VAR_ids,
#' and samples as specified by the `cohort` parameter.
#' The `checkPloidy` parameter can be set to `GRCh37` (or `hg19`) or `GRCh38` (or `hg38`) to
#' assign variant ploidy (diploid,XnonPAR,YnonPAR).
#' 
#' @param object an object of class \code{\link[=gdb-class]{gdb}}
#' @param varSet [`varSet`] object. If specified, the `VAR_id`, `varSetName` and `unit` parameters will be ignored.
#' @param VAR_id Charaocter vector containing target VAR_id.
#' @param ranges Extract variants within specified ranges. 
#' Ranges can be specified as a data.frame, including at least 'CHROM','start', and 'end' columns, or
#' can be a [`GenomicRanges::GRanges`] object. 
#' @param cohort Optional use of cohort data previously uploaded to the gdb (see uploadCohort). 
#' If a valid cohort name is provided, then the uploaded data for this cohort is used to filter and annotate the returned genoMatrix object. If a dataframe is provided, then this is assumed to conform to the SM table constraints required for genoMatrix objects (see genoMatrix).
#' @param checkPloidy Version of the human genome to use when assigning variant ploidy (diploid,XnonPAR,YnonPAR). 
#' Accepted inputs are GRCh37, hg19, GRCh38, hg38. 
#' If no value is provided then all variants are assigned the default ploidy of "diploid"
#' @param varSetName Optional name for the set of variants, for example: 'missense or 'LOF' (ignored if `varSet` is specified.)
#' @param unit Optional 'unit' name, for example: 'SOD1' or 'ENSG00000142168' (ignored if `varSet` is specified.)
#' @param verbose Should the method be verbose? Defaults to `TRUE`.
#' @return A [`genoMatrix`] object.
#' @export
setGeneric("getGT", function(object, varSet = NULL, VAR_id = NULL, ranges = NULL, cohort = NULL, checkPloidy = NULL, varSetName = "unnamed", unit = "unnamed", verbose = TRUE) standardGeneric("getGT"))


#' Generate subset of gdb, retaining all tables. 
#'
#' Function to allow for generation of a child gdb from a parent gdb, with the option to filter retained variants through table intersection operations and SQL where statements
#'
#' @param object A [`gdb`] object.
#' @param output Output gdb path (output will be a new gdb file).
#' @param intersection Additional tables to filter through intersection (ie variants absent from intersection tables will not appear in output). Multiple tables should be ',' delimited.
#' @param where An SQL compliant where clause to filter output; eg: "CHROM=2 AND POS between 5000 AND 50000 AND AF<0.01 AND (cadd.caddPhred>15 OR snpEff.SIFT='D')".
#' @param skipIndexes Flag to skip generation of indexes for var and dosage table (VAR_id;CHROM, POS,REF,ALT). 
#' Typically only required if you plan to use gdbConcat to concatenate a series of separately generated gdb files before use.
#' Defaults to `FALSE`.
#' @param overWrite Flag indicating whether `output` should be overwritten if it already exists.
#' Defaults to `FALSE`.
#' @export
setGeneric("subsetGdb", function(object, output, intersection=c(),where=c(),skipIndexes=FALSE,overWrite=FALSE) standardGeneric("subsetGdb"))

#' Upload variant annotation into gdb.
#'
#' Function to upload variant annotation data into gdb. Assignment of VAR_id is performed automatically through mapping to the var table, but can be skipped using the skipRemap option. Indexing of the imported table based on VAR_id is also automated but can also be skipped, this could be desirable if it is intended to concatenate many separately generated gdb as otherwise intermediate indexes are generated unecessarily.
#'
#' @param object A [`gdb`] object.
#' @param name Name to assign to annotation table.
#' @param value Input variant annotation table. Can be a data frame object or a path to a text file.
#' @param sep Field delimiter (only used if value is a text file). Defaults to `\\t`.
#' @param skipRemap Flag indicating whether to skip mapping of records to VAR_id using CHROM,POS,REF,ALT. Defaults to `FALSE`.
#' @param skipIndexes Flag indicating whether to skip indexing of imported table. Defaults to `FALSE`.
#' @param ignoreAlleles Flag indicating whether to consider REF and ALT allele during mapping of records to VAR_id or just CHROM,POS. Defaults to `FALSE`.
#' @param mapRef Name of lookup table for VAR_id assignment. Defaults to "var".
#' @export
setGeneric("uploadAnno", function(object,name,value,sep="\t",skipRemap=FALSE,skipIndexes=FALSE,ignoreAlleles=FALSE,mapRef="var") standardGeneric("uploadAnno"))

#' mapVariants
#' 
#' Method to map the variants in a [`gdb`] to a set of ranges or features. 
#' The input can be a set of ranges (CHROM, start, end), a bed-file or a gff/gtf-file. 
#' Variants in the gdb will be mapped onto those ranges and annotated with the features/columns 
#' included in the input file. 
#' For example, variants can be easily mapped upon genomic features downloaded in gff format from ensembl. 
#' The output can be written to disk  (`output` parameter) or directly uploaded to the [`gdb`] (`uploadName` parameter). 
#' 
#' @param object a [`gdb`] object
#' @param ranges Can be 1) a data.frame, including at least 'CHROM','start', and 'end' columns.
#' 2) a [`GenomicRanges::GRanges`] object. 3) a filepath to a ranges file containing at least 'CHROM','start', and 'end' columns.
#' Separator can be specified using the `sep` parameter (defaults to `\\t`).
#' @param gff Path to a gff- or gtf-file. 
#' @param bed Path to a bed-file. Specify extra columns using the `bedCols` parameter.
#' @param bedCols A character vector of names of the extra columns to read from the BED-file. 
#' Optionally the vector can be a named vector to indicate the classes of the columns (i.e. c("gene_id" = "character", "gene_name"="character")).
#' If not named, all extra columns will be read as character columns (see examples).
#' @param fields Feature fields to keep. Defaults to `NULL` in which case all fields are kept.
#' @param uploadName Name of table to upload to the gdb. 
#' If not specified, either specifiy `output` to 
#' write the results to disk, or otherwise the results will be returned in the R session.
#' @param output Optionally, an output file path. Can be used instead of `uploadName` to write the results to disk.
#' @param sep Field separator, relevant if `ranges` is a filepath. Defaults to `\\t`. 
#' @param skipIndexes Flag indicating whether to skip indexing of imported table. 
#' Relevant if `uploadName` is specified, and thus the output table is imported in the gdb.
#' Defaults to `FALSE`.
#' @param overWrite if `uploadName` is specified, should an existing table in the gdb with the same name be overwitten?
#' Defaults to `FALSE`.
#' @param verbose Should the method be verbose? Defaults to `TRUE`.
#' @export
setGeneric("mapVariants", function(object, 
                                   ranges = NULL, 
                                   gff = NULL,
                                   bed = NULL,
                                   bedCols = character(),
                                   fields = NULL, 
                                   uploadName = NULL,
                                   output = NULL, 
                                   sep = "\t",
                                   skipIndexes = FALSE,
                                   overWrite=FALSE,
                                   verbose = TRUE) standardGeneric("mapVariants"))

#' Upload cohort table into gdb.
#'
#' Function to upload cohort data tables to gdb. These will automatically be reformatted and sorted to match the ordering of samples in the gdb genotype records.
#'
#' @param object A [`gdb`] object.
#' @param name Name to assign to cohort.
#' @param value Input data frame or a valid file path. Must contain an 'IID' column matching to SM table and a 'sex' column (0=missing,1=male,2=female).
#' @param sep Field delimiter (applies only when value is a text file). Defaults to `\\t`.
#' @export
setGeneric("uploadCohort", function(object,name,value,sep="\t") standardGeneric("uploadCohort"))

#' Drop table from gdb
#'
#' Drop table from [`gdb`] and clear from annotation / cohort metadata tables.
#' @param object [`gdb`] object.
#' @param name  Name of table to drop.
#' @export
setGeneric("dropTable", function(object, name) standardGeneric("dropTable"))

# gdbUtils --------------------------------------------------------------------------

setGeneric("populateGdbFromVcf", function(object, vcf) standardGeneric("populateGdbFromVcf"))

setGeneric("insertVarRecord", function(object, record) standardGeneric("insertVarRecord"))

setGeneric("insertDosageRecord", function(object,record) standardGeneric("insertDosageRecord"))

setGeneric("addRangedVarinfo", function(object, overwrite = TRUE) standardGeneric("addRangedVarinfo"))

# genoMatrix -------------------------------------------------------------------

setGeneric(".resetSexChromDosage", function(object) standardGeneric(".resetSexChromDosage"))

#' @rdname genoMatrix
#' @usage NULL
#' @export
setGeneric("getAF", function(object) standardGeneric("getAF"))

#' @rdname genoMatrix
#' @usage NULL
#' @export
setGeneric("getMAF", function(object) standardGeneric("getMAF"))

#' @rdname genoMatrix
#' @usage NULL
#' @export
setGeneric("getAC", function(object) standardGeneric("getAC"))

#' @rdname genoMatrix
#' @usage NULL
#' @export
setGeneric("getMAC", function(object) standardGeneric("getMAC"))

#' @rdname genoMatrix
#' @usage NULL
#' @export
setGeneric("getNCarriers", function(object) standardGeneric("getNCarriers"))

#' @rdname genoMatrix
#' @usage NULL
#' @export
setGeneric("getCR", function(object, var = TRUE) standardGeneric("getCR"))

##' Return variant summaries
#' 
#'  Returns a per variant summary of genotype counts, frequencies, call-rates and hwe testing.
#'  Note, the [`gdb`] implementation is described here, `summariseGeno` can also be run directly on a 
#'  [`genoMatrix`] object as described in the  [`genoMatrix`] documentation.
#' 
#' @param object a [`gdb`] object
#' @param cohort If a valid cohort name is provided, then the uploaded data for this cohort is used to filter and annotate the genotypes 
#' If not specified, all samples in the gdb will be loaded.
#' @param varSet a [`varSetList`] or [`varSetFile`] object.
#' @param VAR_id A list of VAR_ids, alternatively the varSet parameter can be specified.
#' If single variant tests are ran, the `memlimit` argument controls how many variants to analyze at a time.
#' @param pheno colData field to test as response variable, although not used within this method,
#' this can be useful to filter samples which have missing data for the response variable.
#' @param memlimit Maximum number of variants to load at once (if `VAR_id` is specified).
#' @param geneticModel Which genetic model to apply? ('allelic', 'recessive' or 'dominant').
#' Defaults to `allelic`.
#' @param checkPloidy Version of the human genome to use when assigning variant ploidy (diploid,XnonPAR,YnonPAR). 
#' Accepted inputs are GRCh37, hg19, GRCh38, hg38. 
#' If no value is provided then all variants are assigned the default ploidy of "diploid".
#' @param keep vector of sample IDs to keep, defaults to `NULL`, in which case all samples are kept.
#' @param output Output file path for results.
#' Defaults to `NULL`, in which case results are not written.
#' @param splitBy Split variant summaries by labels indicated in the field 
#' specified by `splitBy`. 
#' @param minCallrateVar Minimum genotype rate for variant retention.
#' @param maxCallrateVar Maximum genotype rate for variant retention.
#' @param minCallrateSM Minimum genotype rate for sample retention.
#' @param maxCallrateSM Maximum genotype rate for sample retention.
#' @param minMAF Minimum minor allele frequency for variant retention.
#' @param maxMAF Maximum minor allele frequency for variant retention.
#' @param minMAC Minimum minor allele count for variant retention.
#' @param maxMAC Maximum minor allele count for variant retention.
#' @param minCarriers Minimum carrier count for variant retention.
#' @param maxCarriers Maximum carrier count for variant retention.
#' @param minCarrierFreq Minimum carrier frequency for variant retention.
#' @param maxCarrierFreq Maximum carrier frequency for variant retention.
#' 
#' @return 
#' Returns a `data.frame` with the following columns:
#' \itemize{
#'   \item \code{VAR_id}:  VAR_id of the respective variant.
#'   \item \code{AF}:  Allele frequency
#'   \item \code{callRate}:  callRate
#'   \item \code{geno0}:  Number of samples with genotype='0'. 
#'   When `geneticModel`='allelic' or 'dominant' this is the number of individuals that are 
#'   homozygous for the reference allele.
#'   \item \code{geno1}:  Number of samples with genotype='1'. 
#'   When `geneticModel`='allelic' this is the number of individuals that are 
#'   heterozygous for the reference allele.
#'   When `geneticModel` = 'dominant' this represents the number of individuals who carry at least one
#'   alternate allele.
#'   When `geneticModel` = 'recessive' this represents the number of individuals who are homozygous 
#'   for the alternate allele.
#'   \item \code{geno2}: When `geneticModel` = 'allelic', the number of individuals who are homozygous
#'   for the alternate allele. 
#'   }
#' @export
setGeneric("summariseGeno", function(object,...) standardGeneric("summariseGeno"))

#' @rdname genoMatrix
#' @usage NULL
#' @export
setGeneric("updateGT", function(object, SM = NULL, anno = NULL) standardGeneric("updateGT"))

#' @rdname genoMatrix
#' @usage NULL
#' @export
setGeneric("flipToMinor", function(object) standardGeneric("flipToMinor"))

#' Recode genotypes
#'
#' Returns a recoded genoMatrix object, genetic model, imputation, and weights can be recoded.
#' @param object An object of class \code{\link{genoMatrix}}
#' @param geneticModel Accepted values include "allelic", "dominant", "recessive".
#' @param imputeMethod Logic for recoding of missing genotypes. The default, "meanImpute" recodes NA to the mean genotype dosage. "missingToRef" recodes NA to 0.
#' @param weights A vector of weights, with length equal to the number of variants.
#' @param MAFweights Apply MAF-dependent weighting, current options include "none" and "mb" (Madsen-Browning).
#' @export
setGeneric("recode", function(object, geneticModel, imputeMethod, weights, MAFweights) standardGeneric("recode"))

#' @rdname genoMatrix
#' @usage NULL
#' @export
setGeneric("aggregateGT", function(object, returnGT = TRUE) standardGeneric("aggregateGT"))

#' @rdname genoMatrix
#' @usage NULL
#' @export
setGeneric("getCarriers", function(object, VAR_id = NULL, colDataFields = NULL, rowDataFields = NULL, groupBy=NULL, aggregate=FALSE, imputeMethod="meanImpute") standardGeneric("getCarriers"))

# varSet -------------------------------------------------------------------

#' Return the units included in an varSetFile/varSetList/aggregateFile
#' 
#' Returns a vector of all units included in a [`varSetFile`], [`varSetList`] or [`aggregateFile`]
#' @export
setGeneric("listUnits", function(object) standardGeneric("listUnits"))

#' Return names of varsets included in a varSetFile or varSetList
#' 
#' Returns a vector of all varSetNames included in the [`varSetFile`] or [`varSetList`]
#' @export
setGeneric("listVarSets", function(object,...) standardGeneric("listVarSets"))

#' Return names of weights included in a varSet or geneSet
#' 
#' Returns a vector of weights included in the [`varSet`] or [`geneSet`]
#' @export
setGeneric("listWeights", function(object,...) standardGeneric("listWeights"))

#' @rdname varSet
#' @usage NULL
#' @export
setGeneric("listVars", function(object,...) standardGeneric("listVars"))

#' Retrieve varSet records from a varSetFile or varSetList
#'
#' Retrieve varSet record(s) for specified unit(s) and/or varSetName(s).
#' `object` can be either a [`varSetList`] or a [`varSetFile`].
#' @param object A [`varSetList`] or [`varSetFile`].
#' @param unit Vector of unit(s) to retrieve 
#' Use [`listUnits()`] to see which units are included in the varSetList/varSetFile.
#' @param varSetName Vector of unit(s) to retrieve 
#' Use [`listVarSets()`] to see which varSets are included in the varSetList/varSetFile.
#' @export
setGeneric("getVarSet", function(object, unit = NULL, varSetName = NULL) standardGeneric("getVarSet"))

#' Generate weighted variant sets for use in association testing.
#'
#' Generate weighted variant sets using annotation table(s) upload to the gdb.
#' See tutorials for examples.
#'
#' @param gdb a [`gdb`] object.
#' @param varSetName Name to assign varSet grouping. This identifier column is used to allow for subsequent mergeing of multiple varSet files for coordinated analysis of multiple variant filtering/ weighting strategies)
#' @param unitTable Table containing aggregation unit mappings.
#' @param unitName Field to utilize for aggregation unit names.
#' @param intersect Additional tables to filter through intersection (ie variants absent from intersection tables will not appear in output). Multiple tables should be ',' delimited.
#' @param where An SQL compliant where clause to filter output; eg: "CHROM=2 AND POS between 5000 AND 50000 AND AF<0.01 AND (cadd.caddPhred>15 OR snpEff.SIFT='D')".
#' @param weightName Field name for desired variant weighting, must be a column within unitTable or other intersection table. Default value of 1 is equivalent to no weighting.
#' @param output Output file name (output will be gz compressed text).
#' @export
setGeneric("buildVarSet", function(
  object, ...) standardGeneric("buildVarSet"))

#' collapseVarSetList
#' @rdname varSetList
#' @usage NULL
#' @export
setGeneric("collapseVarSetList", function(object,...) standardGeneric("collapseVarSetList"))

#' Generate variant sets based on spatial clustering
#'
#' Generate weighted variant sets for use in association testing, with partitioning by genomic distances as described (Fier, GenetEpidemiol, 2017).
#' 
#' @param object a [`gdb`] object
#' @param varSetName Name to assign varSet grouping. This identifier column is used to allow for subsequent mergeing of multiple varSet files for coordinated analysis of multiple variant filtering/weighting strategies)
#' @param unitTable Table containing aggregation unit mappings.
#' @param unitName Field to utilize for aggregation unit names.
#' @param windowSize Numeric vector to indicate starting fixed window size (number of variants)
#' @param overlap Numeric vector to indicate starting fixed window overlap (number of variants, length must match windowSize)
#' @param intersect Additional tables to filter through intersection (ie variants absent from intersection tables will not appear in output). Multiple tables should be ',' delimited.
#' @param where An SQL compliant where clause to filter output; eg: "CHROM=2 AND POS between 5000 AND 50000 AND AF<0.01 AND (cadd.caddPhred>15 OR snpEff.SIFT='D')".
#' @param weightName Field name for desired variant weighting, must be a column within unitTable or other intersection table. Default value of 1 is equivalent to no weighting.
#' @param posField Column name to take as variants position. Default is 'POS' which typically corresponds to genomics position. Can be reset to use CDS or other coordinates. "HGVSc" is a recognized identifier and CDS coordinates will be extracted automatically.
#' @param minTry Minimum number of variants in varset to perform clustering on. If number of variants <minTry, all variants will be returned as a single cluster.
#' Defaults to 5.
#' 
#' @references
#' Loehlein Fier, H. et al. On the association analysis of genome-sequencing data: A spatial clustering approach for partitioning the entire genome into nonoverlapping windows: F ier et al . Genet. Epidemiol. 41, 332–340 (2017).
#' @export
setGeneric("spatialClust", function(object,output,varSetName,unitTable,unitName,windowSize,overlap,intersection=c(),where=c(),weightName=1, posField="POS",minTry=5) standardGeneric("spatialClust"))

# assocTest -------------------------------------------------------------------

#' Perform association tests on binary or quantitative traits
#'
#' @description
#' Main method for performing association tests on binary or quantitative traits. 
#' Performs either aggregate tests or single variant tests, with a wide range of 
#' available statistical methods (see 'Details' for more information). 
#' Various parameters allow for implementation of for example recessive/dominant analyses,
#' MAF-weighted burden tests and permutation-based association tests.
#' Variant and sample filters (e.g. maximum MAF, minimum number of carriers) can be specified.
#' assocTest can be run either directly on a [`genoMatrix`] object or alternatively 
#' on a [`gdb`] object, in which case assocTest can loop through multiple variant sets by providing
#' a [`varSetFile`] or [`varSetList`] as input. 
#' Moreover, assocTest can be run on sets of genes ('gene set burden'), 
#' for which gene burden scores have been generated using [`aggregate`] and stored in an [`aggregateFile`].
#' For running assocTest on a [`genoMatrix`] object see: [`assocTest-genoMatrix`] for details\cr
#' for running assocTest on a [`gdb`]  object see: [`assocTest-gdb-method`] for details
#' for running assocTest on a [`aggregateFile`]  object see: [`assocTest-aggregateFile`] for details
#' 
#' @details
#' **Aggregate tests**
#' For aggregate tests, the following tests are implemented:
#' 
#' *binary/quantitative traits*:
#' * `skat`: SKAT test as implemented in the `SKAT` R package (Wu \emph{et. al}, 2011). 
#' * `skat_burden`: Burden test as implemented in the `SKAT` R package (Wu \emph{et. al}, 2011). 
#' * `skato`: SKAT-O (Optimal Unified Test) as implemented in the `SKAT` R package (Wu \emph{et. al}, 2011). 
#' * `acatv`: ACAT-V test as implemented in the `ACAT` R package (Liu \emph{et al.}, 2019)
#'
#' *quantitative traits*:
#' * `lm`: linear model
#' 
#' *binary traits*:
#' * `firth`: firth logistic regression (Firth, 1993). The maximum number of iterations
#' can be specified using the `maxitFirth` parameter 
#' * `glm`: logistic regression
#' * `nbinom`: Negative binomial test
#' * `skat_robust`: robust SKAT test (robust in the presence of an unbalanced case/control ratio) as implemented in the `SKAT` R package (Zhao \emph{et al.}, 2020).
#' * `skat_burden_robust`: robust burden test (robust in the presence of an unbalanced case/control ratio) as implemented in the `SKAT` R package (Zhao \emph{et al.}, 2020).
#' * `skato_robust`: robust SKAT-O test (robust in the presence of an unbalanced case/control ratio) (Zhao \emph{et al.}, 2020)
#' * `acatvSPA`: adjusted ACAT-v test, using a score test using saddlepoint approximation.
#' 
#' **Single variant tests**
#' 
#' For single variants, the following tests are implemented:
#' 
#' *quantitative traits*:
#' * `lm`: linear model
#' 
#' *binary traits:*
#' * `scoreSPA`: score test with saddle point approximation.
#' * `firth`: firth logistic regression (Firth, 1993). The maximum number of iterations
#' can be specified using the `maxitFirth` parameter.
#' * `glm`: logistic regression
#' 
#' @return 
#' The object returned is of class [`rvbResult`] or [`singlevarResult`]. 
#' The columns are described below: (note that some columns are specific to either `rvbResult` or `singlevarResult`
#' respectively)
#' \itemize{
#'   \item \code{cohort}:  Name of the cohort
#'   \item \code{varSetName}:  varSetName included in the provided `varSet`, if provided, or from the `genoMatrix` metadata.
#'   \item \code{name}:  if specified, contains the `name` argument specified in assocTest. 
#'   \item \code{unit}:  (rvbResult) name of the unit tested (from `varSet` or from `genoMatrix` metadata)
#'   \item \code{VAR_id}:  (singlevarResult) VAR_id of the respective variant.
#'   \item \code{pheno}:  Phenotype tested, as specified by the `pheno` argument.
#'   \item \code{covar}:  Covariates included, as specified by the `covar` argument.
#'   \item \code{geneticModel}:  genetic model used, as specified by the `geneticModel` argument.
#'   \item \code{MAFweight}: MAF-weighting applied, as specified by the `MAFweights` argument.
#'   \item \code{test}:  statistical test used, as specified by the `test` argument.
#'   \item \code{nvar}:  (rvbResult) number of variants included in the aggregate test.
#'   \item \code{caseMAC}:  (singlevarResult) Minor allele count among cases 
#'   \item \code{ctrlMAC}:  (singlevarResult) Minor allele count among controls.
#'   \item \code{caseCarriers}:  (rvbResult) The number of cases that carry at least one variant among variants in the variant set.
#'   \item \code{ctrlCarriers}:  (rvbResult) The number of controls that carry at least one variant among variants in the variant set.
#'   \item \code{meanCaseScore}:  (rvbResult) Mean burden score among cases.
#'   \item \code{meanCtrlScore}:  (rvbResult) Mean burden score among controls.
#'   \item \code{caseN}:  Number of cases.
#'   \item \code{ctrlN}:  Number of controls.
#'   \item \code{caseCallRate}:  Variant call-rate in cases. For aggregate tests this is the average call-rate across variants.
#'   \item \code{ctrlCallRate}:  Variant call-rate in controls For aggregate tests this is the average call-rate across variants.
#'   \item \code{P}:  P-value for specified statistical test. 
#'   \item \code{effect}:  Effect estimate for specified statistical test. For `glm` and `firth` tests, 
#'   the `effect` represents the odds-ratio.
#'   Note that SKAT tests, ACAT-v tests, SPA tests and negative binomial tests don't yield effect estimates.
#'   \item \code{effectSE}:  Standard error of effect estimate of specified statistical test. 
#'   Note that SKAT tests, ACAT-v tests, SPA tests and negative binomial tests don't yield effect estimates.
#'   \item \code{effectCIlower}:  Lower confidence interval of effect estimate of specified statistical test. 
#'   Note that SKAT tests, ACAT-v tests, SPA tests and negative binomial tests don't yield effect estimates.
#'   \item \code{effectCIupper}:  Upper confidence interval of effect estimate of specified statistical test. 
#'   Note that SKAT tests, ACAT-v tests, SPA tests and negative binomial tests don't yield effect estimates.
#' }
#' @usage NULL
#' @references
#' Wu, M. C. et al. Rare-Variant Association Testing for Sequencing Data with the Sequence Kernel Association Test. The American Journal of Human Genetics 89, 82–93 (2011).
#' 
#' Liu, Y. et al. ACAT: A Fast and Powerful p Value Combination Method for Rare-Variant Analysis in Sequencing Studies. The American Journal of Human Genetics 104, 410–421 (2019).
#' 
#' Firth, D. Bias Reduction of Maximum Likelihood Estimates. Biometrika 80, 27–38 (1993).
#' 
#' Zhao, Z. et al. UK Biobank Whole-Exome Sequence Binary Phenome Analysis with Robust Region-Based Rare-Variant Test. The American Journal of Human Genetics 106, 3–12 (2020).

#' @export
setGeneric("assocTest", function(object,
                                 pheno,
                                 test,
                                 ...
) standardGeneric("assocTest"))


# resampling -------------------------------------------------------------------

# rvatResult -------------------------------------------------------------------

#' Write an rvatResult to disk.
#' 
#' Write an [`rvatResult`] object to disk.
#' 
#' @param object \code{\link{rvatResult}} object
#' @param file either a character string naming a file or a connection open for writing. "" indicates output to the console.
#' @param append logical. Only relevant if file is a character string. 
#' If TRUE, the output is appended to the file. If FALSE, any existing file of the name is destroyed.
#' @param quote a logical value (TRUE or FALSE) or a numeric vector. 
#' If TRUE, any character or factor columns will be surrounded by double quotes. 
#' If a numeric vector, its elements are taken as the indices of columns to quote. 
#' In both cases, row and column names are quoted if they are written. 
#' If FALSE, nothing is quoted.
#' @param sep the field separator string. Values within each row of x are separated by this string.
#' @param eol the character(s) to print at the end of each line (row). 
#' For example, eol = "`\\r\\n`" will produce Windows' line endings on a Unix-alike OS, 
#' and eol = "`\\r`" will produce files as expected by Excel:mac 2004.
#' @param na the string to use for missing values in the data.
#' @param dec the string to use for decimal points in numeric or complex columns: must be a single character.
#' @param col.names either a logical value indicating whether the column names of x are to be written along with x, 
#' or a character vector of column names to be written. 
#' See the section on ‘CSV files’ for the meaning of col.names = NA.
#' @param qmethod a character string specifying how to deal with embedded double quote characters when quoting strings. 
#' Must be one of "escape" (default for write.table), i
#' n which case the quote character is escaped in C style by a backslash, or "double" (default for write.csv and write.csv2), 
#' in which case it is doubled. You can specify just the initial letter.
#' @param fileEncoding character string: if non-empty declares the encoding to be used on a file (not a connection) so the character data can be re-encoded as they are written.
#' @export
setGeneric("writeResult", function(object, 
                                   file = "",
                                   append = FALSE,
                                   quote = FALSE,
                                   sep = "\t",
                                   eol = "\n",
                                   na = "NA",
                                   dec = ".",
                                   col.names = TRUE,
                                   qmethod = c("escape", "double"),
                                   fileEncoding = "") standardGeneric("writeResult"))

#' Generate a qq-plot
#' 
#' Plot a qqplot based on the P-values stored in an \code{\link{rvatResult}} object.
#' @param object \code{\link{rvatResult}} object.
#' @param title Optional title.
#' @param label column that significant results should be labeled with. 
#' @param threshold Provide alternative fwe threshold (bonferroni is applied by default.)
#' @param labelThreshold Provide alternative fwe threshold (bonferroni is applied by default.)
#' @param cex Font size, defaults to 16.
#' @param lambda1000 Show lambda1000? Defaults to `FALSE`. 
#' If `TRUE`, the maximum number of cases and controls as specified in the \code{\link{rvatResult}} object 
#' will be used to calculate lambda1000, unless the number of cases and controls are specified
#' using the `case` and `control` parameters.
#' @param case Number of cases in association tests.
#' @param control Number of controls in association tests.
#' 
#' @export
setGeneric("qqplot", function(object, title = "", label = "label", threshold = NULL, showThreshold = TRUE, labelThreshold = NULL,cex = 16, lambda1000 = FALSE, case = NULL, control = NULL) standardGeneric("qqplot"))

#' Generate manhattan plot
#' 
#' Generate a manhattan plot for an \code{\link{rvatResult}} object.
#' @param object \code{\link{rvatResult}} object.
#' @param highlight vector of units/VAR_ids that should be highlighted in red
#' Either give a vector directly or use the method 'manhattanGeneSet'
#' @param label Column name of the labels that are used for significant results
#' @param threshold Fwe threshold (bonferroni applied by default).
#' @param labelThreshold P-value threshold to display labels.
#' @param contigs Update contig lengths from GRCh37 defaults.
#' @param title Optional title.
#' @export
setGeneric("manhattan", function(object, highlight = NULL, label = "label", threshold = NULL, labelThreshold = NULL, labelRepel=FALSE, contigs=c(), title="") standardGeneric("manhattan"))

#' Generate a manhattan plot with a set of genes highlighted
#' 
#' Generate a manhattan plot where the genes of the specified gene set are highlighted and the P-values of each GSA of that gene set is given
#' 
#' @param object a \code{\link{rvbResult}} object
#' @param geneSet a string with the name of a \code{\link{geneSet}} object
#' @param geneSetList a list of \code{\link{geneSet}} objects
#' @param gsaObject a \code{\link{gsaResult}} object
#' @export
setGeneric("manhattanGeneSet", function(object, geneSet, geneSetList, gsaObject, label = 
                                          "label",threshold = NULL, labelThreshold = NULL, contigs=c(), title="") 
  standardGeneric("manhattanGeneSet")) 

#' Density plot
#' 
#' Generate a density plot for a \code{\link{rvatResult}} highlighting the Z-scores of the selected gene set
#' 
#' @param object a \code{\link{rvbResult}} object
#' @param geneSet string with the name of the gene set of interest
#' @param geneSetList list with the available gene sets
#' @param showMeans Show mean Z-score for genes in and outside the gene set. Defaults to `FALSE`.
#' @param INT Apply inverse normal transformation to Z-scores? Defaults to `FALSE`. 
#' @param Zcutoffs A vector (length=2, minimum and maximum) of cutoffs to apply to the Z-scores. 
#' Z scores below/above these cutoffs will be set equal to the cutoff.
#' @param title the title of the plot
#' @export
setGeneric("densityPlot", function(object, geneSet, geneSetList, showMeans = FALSE, INT = FALSE, Zcutoffs = NULL, title = "") standardGeneric("densityPlot"))

#' forestplot
#' Makes a forest plot of the effect size of all occurrences of the chosen gene set
#' 
#' @param object \code{\link{rvatResult}} object
#' @param unit string or vector with names of the units (rvbResult or singlevarResult) or gene sets (gsaResult) of interest
#' @param class the type of rvatResult object: rvbResult, singlevarResult, or gsaResult.
#' @export
setGeneric("forestplot",function(object, unit, class = NULL) standardGeneric("forestplot"))

#' Combine P-values using the ACAT method.
#'
#' Combine P-values in an [`rvbResult`] object using the ACAT method.
#' For details on the ACAT method see: (Liu \emph{et al.}, 2019).
#' @param object \code{\link{rvatResult}} object
#' @param aggregate Variable to ACAT P-values across. 
#' For example, if aggregate = 'test', P-values across statistical tests will be combined using ACAT.
#' A vector with multiple column names can be specified. 
#' For example, by specifying aggregate = c('test', 'varSetName'), P-values across statistical tests and 
#' varSets will be combined using ACAT. 
#' A list with multiple items can be specified, in which case step-wise ACAT is applied.
#' For example, by specifying aggregate = list('test', 'varSetName'), 
#' first P-values across statistical tests are combined, and the resulting P-values are 
#' then combined across varSets. 
#' @param group Variables to group by. 
#' For example, if `group = c('unit', 'varSetName')` and `aggregate = 'test'`, 
#' for each unit-varSetName combination, P-values across statistical tests are combined.
#' Defaults to c("unit", "cohort", "varSetName","name", "unit", "pheno", "covar", "geneticModel", "MAFweight", "test"),
#' i.e. all grouping variables in an rvbResult object. 
#' The variable(s) specified in `aggregate`, are excluded from the grouping variables.
#' @param fixpval Should P-values that are exactly zero or one be fixed? Defaults to `TRUE`. 
#' The method used for fixing the P-values can be specified using the `fixpval_method` parameter.
#' @param fixpval_method Method used to fix p-value (if `fixpval = TRUE`). Methods include:
#' 'minmax' = P-values that are exactly 1 are replaced by the maximum value below 1 present in the results; 
#' P-values that are exactly 0 are replaced by the minimum value above 0 in the results.
#' 'manual' = Specify the replacement P-values using `fixpval_maxP` and `fixpval_minP`.
#' 'Liu' = Method recommended by the author of the ACAT R package 
#' (see: https://github.com/yaowuliu/ACAT), 
#' P-values of 1 are replaced by 1-1/d, where d is the number of p-values combined by ACAT.
#' Since no recommendation for P-values of 0 is given, these are replaced with the value specified using `fixpval_minP`.
#' @param fixpval_maxP Replace P-values that are exactly 1 with this P-value if `fixpval_method = 'manual'`.
#' @param fixpval_minP Replace P-values that are exactly 0 with this P-value if `fixpval_method = 'manual'` 
#' or `fixpval_method = 'Liu'`.
#' @param warning Show warnings? Defaults to `TRUE`. 
#' @references
#' Liu, Y. et al. ACAT: A Fast and Powerful p Value Combination Method for Rare-Variant Analysis in Sequencing Studies. The American Journal of Human Genetics 104, 410–421 (2019).

#' @export
setGeneric("ACAT", function(object,
                            aggregate = "test",
                            group = c("unit", "cohort", "varSetName","name", "pheno", "covar", "geneticModel", "MAFweight", "test"),
                            fixpval = TRUE,
                            fixpval_method = c("minmax", "manual", ),
                            fixpval_maxP = 0.99,
                            fixpval_minP = 1e-32,
                            warning = TRUE) standardGeneric("ACAT"))

#' @rdname rvatResult
#' @usage NULL
#' @export
setGeneric("topResult", function(object, n = 10) standardGeneric("topResult"))

## dump method
setGeneric(".dump", function(object,
                            cohort = "SM",
                            what = c("aggregate", "varSummary"),
                            varSet = NULL,
                            VAR_id = NULL,
                            pheno = NULL,
                            memlimit = 1000,
                            geneticModel = "allelic",
                            imputeMethod = "meanImpute",
                            MAFweights = "none",
                            checkPloidy = NULL,
                            keep = NULL,
                            output = NULL,
                            signif = 3,
                            splitBy = NULL,
                            minCallrateVar = 0,
                            maxCallrateVar = Inf,
                            minCallrateSM = 0,
                            maxCallrateSM =Inf,
                            minMAF = 0,
                            maxMAF = 1,
                            minMAC = 0,
                            maxMAC = Inf,
                            minCarriers = 0,
                            maxCarriers = Inf,
                            minCarrierFreq  = 0,
                            maxCarrierFreq = Inf) standardGeneric(".dump"))


# aggregateFile ----------------------------------------------------------------

#' @rdname aggregateFile
#' @usage NULL
#' @export
setGeneric("listSamples", function(object) standardGeneric("listSamples"))

#' @rdname aggregateFile
#' @usage NULL
#' @export
setGeneric("getUnit", function(object, unit) standardGeneric("getUnit"))

#' Merge multiple aggregate files
#' 
#' Merge aggregrateFiles. Two types of merging can be performed
#' 1) `aggregate=TRUE`: Aggregate values across aggregateFiles. This will result in one aggregate score for
#' each sample, representing the aggregate value across aggregate files. 
#' The output will be a two-column matrix including sample IDs and aggregate scores respectively.
#' 2) `aggregate=FALSE`: Generate a new `aggregateFile` including all aggregates across aggregateFiles.
#'
#' @param object an [`aggregateFileList`] object.
#' @param aggregate Aggregate values? Defaults to `TRUE`.
#' @param output Output file name (output will be gz compressed text). 
#' Defaults to `NULL`, in which case a data.frame will be returned.
#' @param verbose Should the function be verbose ? Defaults to `TRUE`.
#' @export
setGeneric("mergeAggregateFiles", function(
  object,
  collapse = TRUE,
  output = NULL,
  verbose = TRUE
) standardGeneric("mergeAggregateFiles"))


# geneSetAssoc --------------------------------------------------------

#' Return the geneSets included in an geneSetFile/geneSetList
#' 
#' Returns a vector of all units included in a [`geneSetFile`] or [`geneSetList`]
#' @export
setGeneric("listGeneSets", function(object) standardGeneric("listGeneSets"))

#' mapToMatrix
#' 
#' Map the units in a [`rvbResult`] to the geneSets contained in a
#' [`geneSetList`] object. 
#' Returns a matrix where the number of rows equals the number of rows in the rvbResults
#' and the columns equal the number of geneSets in the geneSetList.
#' 
#' @param object  a [`geneSetList`] object
#' @param results a [`rvbResult`] object
#' @param ID Name of the column in `results` that corresponds to the IDs in the geneSetList. Defaults to 'unit'.
#' @param sparse Return a sparse matrix? Defaults to `TRUE`.
setGeneric("mapToMatrix", function(object, results, ID = "unit", sparse = TRUE) standardGeneric("mapToMatrix"))

#' Get geneset(s) from a geneSetList/geneSetFile
#' 
#' Subset a [`geneSetList`] or [`geneSetFile`] object based on one or more geneSet names and/or
#' based on specific units (i.e. subset the geneSets which contain the specified units.)
#' 
#' @param object a [`geneSetList`] or [`geneSetFile`] object
#' @param geneSet a vector of geneSets to subset.
#' @param unit a vector of units to subset
#' @export
setGeneric("getGeneSet", function(object, geneSet = NULL, unit = NULL) standardGeneric("getGeneSet"))

#' @rdname dropUnits
#' @usage NULL
#' @export
setGeneric("dropUnits", function(object, unit = NULL) standardGeneric("dropUnits"))

#' @rdname geneSetList
#' @usage NULL
#' @export
setGeneric("remapIDs", function(object, 
                                dict, 
                                targets = NULL, 
                                duplicate_ids = c("keep_all", "random")) standardGeneric("remapIDs"))

#' @describeIn geneSetFile-class as.geneSetList
#' 
#' Convert a geneSetFile to a geneSetList object
#' 
#' @param object a \code{\link[geneSetFile]{geneSetFile-class}} object
#' 
setGeneric("as.geneSetList", function(object) standardGeneric("as.geneSetList"))

setGeneric("checkDuplicates", function(object, stop = TRUE) standardGeneric("checkDuplicates"))

#' geneSetAssoc
#' 
#' Perform self-contained or competitive gene set analyses. 
#' 
#' @param object an [`rvbResult`] object.
#' @param geneSet a [`geneSetList`] or [`geneSetFile`] object.
#' @param scoreMatrix a matrix (rows = genes, columns = features)
#' @param cormatrix a correlation matrix with row and column names corresponding to the units in the rvbResult. 
#' @param condition Perform conditional analyses. Input can be a 1) a [`geneSetList`] or [`geneSetFile`], in which case
#' the genesets specified in the `geneSet` parameter will all be conditioned on the gene sets provided here.
#' 2) a vector of gene set names present in `geneSet`, all genesets specified in the geneSetList/geneSetFile
#' will be conditioned on the genesets specified here.
#' @param covar a vector of covariates
#' @param test Vector of tests to perform. 
#' Currently implemented tests are the competitive tests lm,mlm,fisher and self-contained tests ttest,ztest and ACAT.
#' @param threshold A vector of thresholds for cutoff-based tests (fisher's exact test / glm).
#' @param Zcutoffs A vector (length=2, minimum and maximum) of cutoffs to apply to the Z-scores. 
#' Z scores below/above these cutoffs will be set equal to the cutoff.
#' @param INT Apply inverse normal transformation to Z-scores? Defaults to `FALSE`. 
#' @param minSetSize Exclude genesets with size < minSetSize
#' @param maxSetSize Exclude genesets with size > maxSetSize
#' @param oneSided Calculate a one-sided P-value? Defaults to `TRUE`.
#' @param memlimit Maximum number of genesets to process in one go.
#' @param REML Use REML for mixed linear models? Defaults to `TRUE`.
#' @param ID ID column in the rvbResult that corresponds with the IDs used in the geneSetList.
#' Defaults to 'unit'.
#' @param output Optional: save results to specified path
#' @export
setGeneric("geneSetAssoc", function(object,
                                    geneSet = NULL,
                                    scoreMatrix = NULL,
                                    cormatrix = NULL,
                                    condition = NULL,
                                    covar = NULL,
                                    test = c("lm", "glm", "mlm", "fisher", "ttest", "ztest", "ACAT"),
                                    threshold = NULL,
                                    Zcutoffs = NULL,
                                    INT = FALSE,
                                    scoreCutoffs = NULL,
                                    minSetSize = 1,
                                    maxSetSize = Inf,
                                    oneSided = TRUE,
                                    memlimit = 1000, 
                                    REML = TRUE,
                                    ID = "unit",
                                    output = NULL
) standardGeneric("geneSetAssoc"))