# init ----------------------------------------------------------------------
check_expected_args <- function(expected, args, func_name) {
  if(!all(names(args) %in% expected)) {
    stop(sprintf("Unexpected argument(s) for %s: %s", 
                 func_name, paste(paste0("--", names(args)[!names(args) %in% expected]), collapse=",")),
         call. = FALSE
    )
  }
}

check_required_args <- function(required, required_one_of=NULL, args, func_name) {
  if(!all(required %in% names(args))) {
    stop(sprintf("The following required arguments for %s are missing: %s", 
                 func_name,
                 paste(paste0("--", required[!required %in% names(args)]), collapse=",")),
         call. = FALSE
    )
  }
  if(!is.null(required_one_of)) {
    if(sum(required_one_of %in% names(args)) == 0) {
      stop(sprintf("One of the following arguments should be specified: %s", 
                   paste(paste0("--", required_one_of), collapse=",")),
           call. = FALSE
      )
    }
  }
}

stopQuietly <- function(...) {
  blankMsg <- sprintf("\r%s\r", paste(rep(" ", getOption("width")-1L), collapse=" "));
  stop(simpleError(blankMsg));
} 

check_help <- function(args,  help, func_name) {
  if(length(args) == 1) {
    message(help[[func_name]])
    stopQuietly()
  }
}

help=list(
  buildGdb = "
  buildGdb
  
  Creates a new gdb file and returns a connection object of type gdb-class. The gdb can be structured and populated using a provided vcf, bgen or plink file. 
  If no input variant file is provided then only an empty gdb is created.
  
  Usage:
    Rscript rvat.R --buildGdb --vcf={vcf} --output={output} [options]
 
  Arguments:
    output           Path for output gdb file
    vcf              Input vcf file used to structure and populate gdb. Warning this function makes the following of assumptions: 1) strict adherence to vcf format (GT subfield first element in genotype firelds), 2) multiallelic records have been split, 3) desired genotype QC has already been applied (DP,GQ filters), 4) GT values conform to the set {0/0,0/1,1/0,1/1,./.,0|0,0|1,1|0,1|1,.|.}. Multiallelic parsing and genotype QC can be performed using vcftools and/or accompanying parser scripts included on the rvat github.
    skipIndexes      Flag to skip generation of indexes for var and dosage table (VAR_id;CHROM, POS,REF,ALT). Typically only required if you plan to use gdbConcat to concatenate a series of separately generated gdb files before use
    overWrite        overwrite if `output` already exists? Defaults to `FALSE`, in which case an error is raised.
  ",
  concatGdb="
concatGdb

Function to concatenate gdb databases. Only retains content of base tables (SM, var, dosage).

Usage:
  Rscript rvat.R --concatGdb --targets={targets} --output={output} [options]

Arguments:
  targets        File listing full paths of gdb to concatenate
  output         Output gdb file path.
  skipRemap      Flag to skip reseting of VAR_id to row id after concatenation.
  skipIndexes    Flag to skip generation of standard var and dosage table indexes (VAR_id;CHROM, POS,REF,ALT)

",
  uploadAnno="
uploadAnno

Function to upload variant annotation data into gdb.
Assignment of VAR_id is performed automatically through mapping to the var table, but can be skipped using the skipRemap option. 
Indexing of the imported table based on VAR_id is also automated but can also be skipped, 
this could be desirable if it is intended to concatenate many separately generated gdb as otherwise intermediate indexes are generated unecessarily.

Usage:
  Rscript rvat.R --uploadAnno --gdb={gdb} --name={name} --value={value} [options]

Arguments:
  gdb            gdb file path
  name           Name to assign to annotation table
  value          Full path to annotation table text file
  sep            Field delimiter
  skipRemap      Flag indicating whether to skip mapping of records to VAR_id using CHROM,POS,REF,ALT
  skipIndexes    Flag indicating whether to skip indexing of imported table.
  ignoreAllele   Flag indicating whether to consider REF and ALT allele during mapping of records to VAR_id or just CHROM,POS.
  mapRef         Look up table for VAR_id assignment
",
mapVariants="
mapVariants

Method to map the variants in a gdb to a set of ranges or features. 
The input can be a set of ranges (CHROM, start, end), a bed-file or a gff/gtf-file. 
Variants in the gdb will be mapped onto those ranges and annotated with the features/columns 
included in the input file. 
For example, variants can be easily mapped upon genomic features downloaded in gff format from ensembl. 
The output can be written to disk  (`output` parameter) or directly uploaded to the gdb (`uploadName` parameter). 

Usage:
  Rscript rvat.R --mapVariants --gdb={gdb} --gff={gff} --uploadName=ensembl

Arguments:
  gdb            gdb file path
  ranges         a filepath to a ranges file containing at least 'CHROM','start', and 'end' columns.
                 Separator can be specified using the `sep` parameter (defaults to `\\t`).
  gff            Path to a gff- or gtf-file. 
  bed            Path to a bed-file. Specify extra columns using the `bedCols` parameter.
  bedCols        A comma-delimited list of names of the extra columns to read from the BED-file. 
  fields         A comma-delimited list of feature fields to keep. By default all fields are kept.
  uploadName     Name of table to upload to the gdb. If not specified, specifiy `output` to write results to disk.
  output         Optionally, an output file path. Can be used instead of `uploadName` to write the results to disk.
  sep            Field separator, relevant if `ranges` is specified. Defaults to '\t'. 
  skipIndexes    skipIndexes Flag indicating whether to skip indexing of imported table. 
                 Relevant if `uploadName` is specified, and thus the output table is imported in the gdb.
                 By default the table is indexed on VAR_id.
  overWrite      If `uploadName` is specified, should an existing table in the gdb with the same name be overwitten?
                 By default the method is aborted when the table already exists in the gdb.
  verbose        verbose Should the method be verbose? Defaults to `TRUE`.

",
uploadCohort="
uploadCohort

Function to upload cohort data tables to gdb

Usage:
  Rscript rvat.R --uploadCohort --gdb={gdb} --name={name} --value={value} [options]

Arguments:
  gdb            gdb file path
  name           Name to assign to cohort
  value          Full path to cohort annotation file. Must contain an 'IID' column matching to SM table and a 'sex' column (0=missing,1=male,2=female)
  sep            Field delimiter

",
listAnno="
listAnno

Dump variant annotation metadata

Usage:
  Rscript rvat.R --listAnno --gdb={gdb}

Arguments:
  gdb            gdb file path

",
listCohort="
listCohort

Dump cohort annotation metadata

Usage:
  Rscript rvat.R --listCohort --gdb={gdb}

Arguments:
  gdb           gdb file path

",
dropTable="
dropTable

drop table from gdb file and update metadata.

Usage:
Rscript rvat.R --dropTable --gdb={gdb} --name={name}

Arguments:
  gdb           gdb file path
  name          Name of table to drop

",
subsetGdb="
subsetGdb

Function to allow for generation of a child gdb from a parent gdb with the option to filter retained variants through table intersections and SQL where statements.

Usage:
  Rscript rvat.R --subsetGdb --gdb={gdb} --output={output} [options]

Arguments:
  gdb            gdb file path
  output         Output file name (output will be a new gdb file).
  intersect      Additional tables to filter through intersection (ie variants absent from intersection tables will not appear in output). Multiple tables should be ',' delimited.
  where          An SQL compliant where clause to filter output; eg: \"CHROM=2 AND POS between 5000 AND 50000 AND AF<0.01 AND (cadd.caddPhred>15 OR snpEff.SIFT='D')\"
  skipIndexes    Flag to skip generation of indexes for var and dosage table (VAR_id;CHROM, POS,REF,ALT). Typically only required if you plan to use gdbConcat to concatenate a series of separately generated gdb files before use.
  overWrite      Flag indicating whether `output` should be overwritten if it already exists.
",
buildVarSet="
buildVarSet

Generate weighted variant sets for use in association testing.

Usage:
  Rscript rvat.R --buildVarSet --gdb={gdb} --unitTable={unitTable} --unitName={unitName} --output={output} [options]

Arguments:
  gdb           gdb file path
  varSetName    Name to assign varSet grouping. This identifier column is used to allow for subsequent mergeing of multiple varSet files for coordinated analysis of multiple variant filtering/ weighting strategies)
  unitTable     Table containing aggregation unit mappings
  unitName      Field to utilize for aggregation unit names
  intersect     Additional tables to filter through intersection (ie variants absent from intersection tables will not appear in output). Multiple tables should be ',' delimited
  where         An SQL compliant where clause to filter output; eg: \"CHROM=2 AND POS between 5000 AND 50000 AND AF<0.01 AND (cadd.caddPhred>15 OR snpEff.SIFT='D')\".
  weightName    Field name for desired variant weighting, must be a column within unitTable or other intersection table. Default value of 1 is equivalent to no weighting.
  output        Output file name (output will be gz compressed text)
",

spatialClust="
spatialClust
Generate weighted variant sets for use in association testing, with partitioning by genomic distances as described (Fier, GenetEpidemiol, 2017).
Usage:
  Rscript rvat.R --spatialClust --gdb={gdb} --unitTable={unitTable} --unitName={unitName} --windowSize=25,50 --overlap=10,20 --output={output} [options]
Arguments:
  gdb           gdb file path
  varSetName    Name to assign varSet grouping. This identifier column is used to allow for subsequent mergeing of multiple varSet files for coordinated analysis of multiple variant filtering/ weighting strategies)
  unitTable     Table containing aggregation unit mappings
  unitName      Field to utilize for aggregation unit names
  windowSize    Starting fixed window sizes (number of variants), comma-delimited.
  overlap       Starting fixed window overlap (number of variants, length must match windowSize), comma-delimited.
  intersect     Additional tables to filter through intersection (ie variants absent from intersection tables will not appear in output). Multiple tables should be ',' delimited
  where         An SQL compliant where clause to filter output; eg: \"CHROM=2 AND POS between 5000 AND 50000 AND AF<0.01 AND (cadd.caddPhred>15 OR snpEff.SIFT='D')\".
  weightName    Field name for desired variant weighting, must be a column within unitTable or other intersection table. Default value of 1 is equivalent to no weighting.
  posField      Column name to take as variants position. Default is 'POS' which typically corresponds to genomics position. Can be reset to use CDS or other coordinates. 'HGVSc' is a recognized identifier and CDS coordinates will be extracted automatically.
  minTry        Minimum number of variants in varset to perform clustering on. If number of variants <minTry, all variants will be returned as a single cluster. Default = 5.
  output        Output file name (output will be gz compressed text)
",

assocTest="
  assocTest
  
  Test for assocation between phenotype of interest and (aggregated) variants
  
  Usage:
    Rscript rvat.R --assocTest --gdb={gdb} --varSet={varSet} --cohort={cohort} --pheno={pheno} --test={test} --output={output} [options]
  
  Arguments:
    gdb               gdb file path
    pheno             Cohort field to test as response variable, the response variable can either be binary (0/1) or continuous. If the response variable is continuous, include the `--continuous` argument. Multiple phenotypes can be specified using ',' to delimit them.
    test              Statistical tests to run (',' delimited), options include firth,glm,lm,nbinom,skat,skat_burden,skato,skat_robust,skato_robust,skat_burden_robust, acatv, acatvSPA
    cohort            If a valid cohort name is provided, then the uploaded data for this cohort is used to filter and annotate the genoMatrix object. 
                      If not specified, all samples in the gdb will be loaded.
    varSet            varSetFile path. Format as per varSet: unit|VAR_id|weight|varSetName where VAR_id and weight are in turn ',' delimited lists.
    VAR_id            A list of VAR_ids, alternatively the varSet parameter can be specified.
                      If single variant tests are ran, the `memlimit` argument controls how many variants to analyze at a time.
    name              Optional name for the analysis, defaults to 'none'
    continuous        Flag that indicates that the phenotype(s) is/are continuous.
    singlevar         Flag that indicates that single variant tests should be run.
    covar             Regression covariates, ',' delimited. Multiple sets of covariates can be specified, delimited by '/'.
    geneticModel      Genetic models to test ('allelic', 'recessive', 'dominant'), multiple models can be specified (',' delimited).
    imputeMethod      Imputation method, either 'meanImpute' or 'missingToRef'. If not specified, single variant tests are not imputed, whereas burden tests are mean imputed.
    MAFweights        MAF weighting method. Currently Madsen-Browning ('mb') is implemented, by default no MAF weighting is applied.
    maxitFirth        Maximum number of iterations to use for estimating firth confidence intervals. Defaults to 1000.
    checkPloidy       Version of the human genome to use when assigning variant ploidy (diploid,XnonPAR,YnonPAR). 
                      Accepted inputs are 'GRCh37', 'hg19', 'GRCh38', 'hg38'. 
                      If no value is provided then all variants are assigned the default ploidy of 'diploid'.
    keep              Filepath to list of samples to retain in analysis (by default all samples are kept).
    output            Output file path for results.
    methodResampling  Which method to use for resampling? ('permutation' currently implemented). 
                      Defaults to `NULL`, in which case no resampling is performed. 
    resamplingFile    file path to a resamplingFile (see `buildResamplingFile` function)
    nResampling       Number of resamplings (note this is ignored if a resamplingFile is specified)
    outputResampling  Output resamplings? If not specified, it's assumed that resampled P-values should be calculated.
    memlimitResampling  Number of resamplings to do at a time.
    minCallrateVar    Minimum genotype rate for variant retention.
    maxCallrateVar    Maximum genotype rate for variant retention.
    minCallrateSM     Minimum genotype rate for sample retention.
    maxCallrateSM     Maximum genotype rate for sample retention.
    minMAF            Minimum minor allele frequency for variant retention.
    maxMAF            Maximum minor allele frequency for variant retention.
    minMAC            Minimum minor allele count for variant retention.
    maxMAC            Maximum minor allele count for variant retention.
    minCarriers       Minimum carrier count for variant retention.
    maxCarriers       Maximum carrier count for variant retention.
    minCarrierFreq    Minimum carrier frequency for variant retention.
    maxCarrierFreq    Maximum carrier frequency for variant retention.
    memlimit          Maximum number of variants to load at once (if --VAR_id is specified).
    seed              Seed to test if running permutation analyses.
  ",
summariseGeno="
  summariseGeno
  
  Usage:
    Rscript rvat.R --summariseGeno --gdb={gdb} --VAR_id={VAR_id} --cohort={cohort} --output={output} [options]
  
  Arguments:
    gdb           gdb file path
    cohort        cohort data previously uploaded to the gdb (see uploadCohort).
    varSet        varSetFile path. Format as per varSet: unit|VAR_id|weight|varSetName where VAR_id and weight are in turn ',' delimited lists
    VAR_id        A list of VAR_ids, alternatively the varSet parameter can be specified.
                  If single variant tests are ran, the `memlimit` argument controls how many variants to analyze at a time.
    pheno         Cohort field to test as response variable, although not used within this method, this can be useful to filter samples which have missing data for the response variable.
    memlimit      Maximum number of variants to load at once (if `VAR_id` is specified).
    geneticModel  Genetic model to apply ('allelic', 'recessive', 'dominant'). Defaults to 'allelic'.
    checkPloidy   Version of the human genome to use when assigning variant ploidy (diploid,XnonPAR,YnonPAR). Accepted inputs are GRCh37, hg19, GRCh38, hg38. 
                  If no value is provided then all variants are assigned the default ploidy of 'diploid'.
    keep          Filepath to list of samples to retain in analysis (by default all samples are kept).
    output        Output file path for results.
    splitBy       Split variant summaries by labels indicated in the field specified by `splitBy`. 
    minCallrateVar Minimum genotype rate for variant retention.
    maxCallrateVar Maximum genotype rate for variant retention.
    minCallrateSM Minimum genotype rate for sample retention.
    maxCallrateSM Maximum genotype rate for sample retention.
    minMAF Minimum minor allele frequency for variant retention.
    maxMAF Maximum minor allele frequency for variant retention.
    minMAC Minimum minor allele count for variant retention.
    maxMAC Maximum minor allele count for variant retention.
    minCarriers Minimum carrier count for variant retention.
    maxCarriers Maximum carrier count for variant retention.
    minCarrierFreq Minimum carrier frequency for variant retention.
    maxCarrierFreq Maximum carrier frequency for variant retention.
  ",
aggregate="
  aggregate
  
  Usage:
    Rscript rvat.R --aggregate --gdb={gdb} --varSet={varSet} --cohort={cohort} --output={output} [options]
  
  Arguments:
    gdb           gdb file path
    cohort        cohort data previously uploaded to the gdb (see uploadCohort).
    varSet        varSetFile path. Format as per varSet: unit|VAR_id|weight|varSetName where VAR_id and weight are in turn ',' delimited lists
    VAR_id        A list of VAR_ids, alternatively the varSet parameter can be specified.
                  If single variant tests are ran, the `memlimit` argument controls how many variants to analyze at a time.
    pheno         Cohort field to test as response variable, although not used within this method, this can be useful to filter samples which have missing data for the response variable.
    memlimit      Maximum number of variants to load at once (if `VAR_id` is specified).
    geneticModel  Genetic model to apply ('allelic', 'recessive', 'dominant'). Defaults to 'allelic'.
    imputeMethod  Which imputation method to apply? ('meanImpute' or 'missingToRef'). Defaults to 'meanImpute'.
    MAFweights    Apply MAF weighting? Currently Madsen-Browning ('mb') is implemented. Defaults to 'none'.
    checkPloidy   Version of the human genome to use when assigning variant ploidy (diploid,XnonPAR,YnonPAR). Accepted inputs are GRCh37, hg19, GRCh38, hg38. 
                  If no value is provided then all variants are assigned the default ploidy of 'diploid'.
    keep          Filepath to list of samples to retain in analysis (by default all samples are kept).
    output        Output file path for results.
    minCallrateVar Minimum genotype rate for variant retention.
    maxCallrateVar Maximum genotype rate for variant retention.
    minCallrateSM Minimum genotype rate for sample retention.
    maxCallrateSM Maximum genotype rate for sample retention.
    minMAF Minimum minor allele frequency for variant retention.
    maxMAF Maximum minor allele frequency for variant retention.
    minMAC Minimum minor allele count for variant retention.
    maxMAC Maximum minor allele count for variant retention.
    minCarriers Minimum carrier count for variant retention.
    maxCarriers Maximum carrier count for variant retention.
    minCarrierFreq Minimum carrier frequency for variant retention.
    maxCarrierFreq Maximum carrier frequency for variant retention.
  ",
mergeAggregateFiles = "
mergeAggregateFiles
  
  merge aggregate files 
  
  Usage:
    Rscript rvat.R --mergeAggregateFiles --filelist={filelist} --output={output} [options]
  
  Arguments:
    filelist   Filepath to list of aggregateFiles.
    collapse   Aggregate values? Defaults to `TRUE`.
    output     Output file name (output will be gz compressed text). 
    verbose    Should the function be verbose ? Defaults to `TRUE`.
    checkDups  check if dupicate units are present in the aggregateFiles? Defaults to `TRUE`.
",
geneSetAssoc="

  geneSetAssoc
  
  Perform gene set analyses
  
  Usage:
    Rscript rvat.R --geneSetAssoc --geneSetFile={geneSetFile} --rvbResult={result} --output={output} [options]
  
  Arguments:
    geneSetfile   file path to a geneSetFile object (generated using the buildGeneSet method)
    rvbResult     file path to an rvbResult object (generated using the assocTest method)
    nullmodel     file path (rds-file) to a nullModelGSA object (generated using the fitNullModel method)
    cormatrix     file path (rds-file) to a correlation matrix
    covar         covariates (comma-delimited)
    test          which tests to perform, competitive analyses include 'lm', 'fisher', and 'mlm'.
                  self-contained analyses include 'ttest', 'ztest' and 'ACAT'.
    threshold     threshold to use for count-based tests (fisher and glm). Defaults to bonferroni threshold.
    Zcutoffs      Cutoffs to apply to the Z-scores (minimum,maximum). Z scores below/above these cutoffs will be set equal to the cutoff.
    minSetSize    Exclude genesets with size < minSetSize
    maxSetSize    Exclude genesets with size > maxSetSize
    oneSided      Calculate a one-sided P-value? Defaults to `TRUE`.
    memlimit      Maximum number of genesets to process in one go, defaults to 1000.
    ID            ID column in the rvbResult that corresponds with the IDs used in the geneSetList. Defaults to 'unit'.
    output        output file path
  ",
buildResamplingFile="
  buildResamplingFile

  Generate a file containing resamplings.

  Usage:
    Rscript rvat.R --buildResamplingFile --nSamples={nSamples} --nResamplings={nResamplings} --methodRsampling={methodResampling} --output={output} [options]
  
  Arguments:
    nSamples  Number of samples
    nResampling  Number of resamplings
    memlimit Chunk sizes
    methodResampling  Resampling method, currently 'permutation' is implemented.
    output  File path (.gz extension) to write output to.
    seed Random seed. Defaults to 10.
",
vcfInfo2Table="
  vcfInfo2Table

  Convert vcf info field to table format. Requires valid vcf where INFO fields are specified in header.

  Usage:
    Rscript rvat.R --vcfInfo2Table --vcf={vcf} --output={output} --splitMultiallelic={splitMultiallelic}
  
  Arguments:
    vcf     Number of samples
    output   output path
    splitMultiallelic Returns one row per alternative allele instead of one row per variant. Default=TRUE.
"
)

message(sprintf("Start time: %s", Sys.time()))
message(sprintf("RVAT version: %s", packageVersion("rvat")))

# Construct general help message
help$general="Rare variant analysis toolkit (rvat)

Provide a valid function call for detailed help."
for (i in names(help))
{
  if (i!="general"){help$general=sprintf("%s
  --%s",help$general,i)}
}


# Parse arguments
args=rvat:::getArgs(help=help$general)

# Validate program selection
if (sum(names(help) %in% names(args))!=1){stop(help$general)}


# importVcf --------------------------------------------------------------------

if (!is.null(args[["buildGdb"]]))
{
  check_help(args = args, help = help, func_name = "buildGdb")
  required <- c("output")
  expected <- c("buildGdb","vcf", "output","skipIndexes", "skipVarRanges", "overWrite","memlimit")
  check_required_args(required = required, args = args, func_name = "buildGdb")
  check_expected_args(expected = expected, args = args, func_name = "buildGdb")
  
  if (is.null(args[["skipIndexes"]]))  {skipIndexes = FALSE } else{ skipIndexes = TRUE }
  if (is.null(args[["skipVarRanges"]]))  {skipVarRanges = FALSE } else{ skipVarRanges = TRUE }
  if (is.null(args[["vcf"]])) {vcf = c()} else {vcf = args[["vcf"]]}
  if (is.null(args[["overWrite"]])) {overWrite = FALSE} else {overWrite = as.logical(args[["overWrite"]])}
  if (is.null(args[["memlimit"]])) {memlimit = 1000} else {memlimit = as.numeric(args[["memlimit"]])}
  
  rvat::buildGdb(vcf = vcf,
                 output = args[["output"]],
                 skipIndexes = skipIndexes,
                 skipVarRanges = skipVarRanges,
                 overWrite = overWrite,
                 memlimit = memlimit
                 )
}

# concatGdb --------------------------------------------------------------------
if (!is.null(args[["concatGdb"]]))
{
  check_help(args = args, help = help, func_name = "concatGdb")
  required=c("targets","output")
  expected=c("concatGdb","targets","output","skipRemap","skipIndexes")
  check_required_args(required = required, args = args, func_name = "concatGdb")
  check_expected_args(expected = expected, args = args, func_name = "concatGdb")
  
  if (is.null(args[["skipRemap"]])) {skipRemap <- FALSE} else{skipRemap <- TRUE}
  if (is.null(args[["skipIndexes"]])) {skipIndexes <- FALSE} else{skipIndexes <- TRUE}
  
  rvat::concatGdb(args[["targets"]],
                  args[["output"]],
                  skipRemap = skipRemap,
                  skipIndexes  = skipIndexes)
}

# uploadAnno -------------------------------------------------------------------

if (!is.null(args[["uploadAnno"]]))
{
  check_help(args = args, help = help, func_name = "uploadAnno")
  required <- c("gdb","name","value")
  expected <- c("uploadAnno","gdb","name","value","sep","skipRemap","skipIndexes","ignoreAlleles","mapRef")
  check_required_args(required = required, args = args, func_name = "uploadAnno")
  check_expected_args(expected = expected, args = args, func_name = "uploadAnno")
  
  gdb <- rvat::gdb(args[["gdb"]])
  if (is.null(args[["sep"]])) {sep <- "\t"} else{sep <- args[["sep"]]}
  if (is.null(args[["skipRemap"]])){skipRemap <- FALSE} else{skipRemap <- TRUE}
  if (is.null(args[["skipIndexes"]])){skipIndexes <- FALSE} else{skipIndexes <- TRUE}
  if (is.null(args[["ignoreAlleles"]])){ignoreAlleles <- FALSE} else{ignoreAlleles <- TRUE}
  if (is.null(args[["mapRef"]])){mapRef <- "var"} else{mapRef <- args[["mapRef"]]}
  
  rvat::uploadAnno(gdb,
                   args[["name"]],
                   args[["value"]],
                   sep = sep,
                   skipIndexes = skipIndexes,
                   skipRemap = skipRemap,
                   ignoreAlleles = ignoreAlleles,
                   mapRef = mapRef)
}

# mapVariants -------------------------------------------------------------------
if (!is.null(args[["mapVariants"]]))
{
  check_help(args = args, help = help, func_name = "mapVariants")
  required <- c("gdb")
  expected <- c("mapVariants","gdb","ranges","gff","bed","fields","uploadName","output","bedCols", "sep", "skipIndexes", "overWrite", "verbose")
  check_required_args(required = required, args = args, func_name = "mapVariants")
  check_expected_args(expected = expected, args = args, func_name = "mapVariants")
  
  gdb <- rvat::gdb(args[["gdb"]])
  if (is.null(args[["sep"]])) {sep <- "\t"} else{sep <- args[["sep"]]}
  if(is.null(args[["fields"]])) fields <- NULL else fields <- unlist(strsplit(args[["fields"]],split=","))
  if (is.null(args[["skipIndexes"]])){skipIndexes <- FALSE} else{skipIndexes <- TRUE}
  if (is.null(args[["overWrite"]])) {overWrite = FALSE} else {overWrite = as.logical(args[["overWrite"]])}
  if (is.null(args[["bedCols"]])) {bedCols = character()} else {bedCols <- unlist(strsplit(args[["bedCols"]],split=","))}
  
  
  rvat::mapVariants(gdb,
                    ranges = args[["ranges"]],
                    gff = args[["gff"]],
                    bed = args[["bed"]],
                    bedCols = bedCols,
                    fields = fields,
                    uploadName = args[["uploadName"]],
                    output = args[["output"]],
                    sep = sep,
                    overWrite = overWrite,
                    skipIndexes = skipIndexes,
                    verbose = if(is.null(args[["verbose"]])) TRUE else as.logical(args[["verbose"]])
                   )
}


# mapToCDS -------------------------------------------------------------------
if (!is.null(args[["mapToCDS"]]))
{
  check_help(args = args, help = help, func_name = "mapToCDS")
  required <- c("gdb", "gff", "output")
  expected <- c("mapToCDS","gdb","gff","exonPadding","output","gene_id","transcript_id","verbose")
  check_required_args(required = required, args = args, func_name = "mapToCDS")
  check_expected_args(expected = expected, args = args, func_name = "mapToCDS")
  
  gdb <- rvat::gdb(args[["gdb"]])
  if(is.null(args[["gene_id"]])) gene_id <- NULL else gene_id <- readLines(args[["gene_id"]])
  if(is.null(args[["transcript_id"]])) transcript_id <- NULL else transcript_id <- readLines(args[["transcript_id"]])
  if(is.null(args[["biotype"]])) biotype <- NULL else biotype <- readLines(args[["biotype"]])
  
  rvat::mapToCDS(gdb,
                 gff = args[["gff"]],
                 exonPadding = if(!is.null(args[["exonPadding"]])) as.numeric(args[["exonPadding"]]) else 12,
                 output = args[["output"]],
                 gene_id = gene_id,
                 transcript_id = transcript_id,
                 biotype = biotype,
                 verbose = if(is.null(args[["verbose"]])) TRUE else as.logical(args[["verbose"]])
  )
}


# uploadCohort -------------------------------------------------------------------
if (!is.null(args[["uploadCohort"]]))
{
  check_help(args = args, help = help, func_name = "uploadCohort")
  required=c("gdb","name","value")
  expected=c("uploadCohort","gdb","name","value","sep")
  check_required_args(required = required, args = args, func_name = "uploadCohort")
  check_expected_args(expected = expected, args = args, func_name = "uploadCohort")
  
  gdb <- rvat::gdb(args[["gdb"]])
  if (is.null(args[["sep"]])) {sep <- "\t"} else{sep <- args[["sep"]]}
  rvat::uploadCohort(gdb,
                     args$name,
                     args$value,
                     sep=sep)
}

# listAnno ---------------------------------------------------------------------
if (!is.null(args[["listAnno"]]))
{
  check_help(args = args, help = help, func_name = "listAnno")
  required=c("gdb")
  expected=c("listAnno","gdb")
  check_required_args(required = required, args = args, func_name = "listAnno")
  check_expected_args(expected = expected, args = args, func_name = "listAnno")
  
  gdb <- rvat::gdb(args[["gdb"]])
  rvat::listAnno(gdb)
}

# listCohort -------------------------------------------------------------------
if (!is.null(args[["listCohort"]]))
{
  check_help(args = args, help = help, func_name = "listCohort")
  required=c("gdb")
  expected=c("listCohort","gdb")
  check_required_args(required = required, args = args, func_name = "listCohort")
  check_expected_args(expected = expected, args = args, func_name = "listCohort")
  
  gdb <- rvat::gdb(args[["gdb"]])
  rvat::listCohort(gdb)
}

# dropTable --------------------------------------------------------------------
if (!is.null(args[["dropTable"]]))
{
  check_help(args = args, help = help, func_name = "dropTable")
  required=c("gdb","name")
  expected=c("dropTable","gdb","name")
  check_required_args(required = required, args = args, func_name = "dropTable")
  check_expected_args(expected = expected, args = args, func_name = "dropTable")
  
  gdb <- rvat::gdb(args[["gdb"]])
  rvat::dropTable(gdb,
                  args[["name"]])
}

# subsetGdb --------------------------------------------------------------------
if (!is.null(args[["subsetGdb"]]))
{
  check_help(args = args, help = help, func_name = "subsetGdb")
  required=c("gdb","output")
  expected=c("subsetGdb","gdb","output","intersection","where","skipIndexes", "overWrite")
  
  check_required_args(required = required, args = args, func_name = "subsetGdb")
  check_expected_args(expected = expected, args = args, func_name = "subsetGdb")
  
  gdb <- rvat::gdb(args[["gdb"]])
  if (is.null(args[["intersection"]])) {intersection=c()} else{intersection = args[["intersection"]]}
  if (is.null(args[["where"]])){where = c()} else{where = args[["where"]]}
  if (is.null(args[["skipIndexes"]])) {skipIndexes = FALSE} else{skipIndexes = TRUE}
  if (is.null(args[["overWrite"]])) {overWrite = FALSE} else{overWrite = TRUE}
  
  rvat::subsetGdb(gdb,
                  output = args[["output"]],
                  intersection = intersection,
                  where = where,
                  skipIndexes = skipIndexes,
                  overWrite = overWrite 
                  )
}

# buildVarSet ------------------------------------------------------------------
if (!is.null(args[["buildVarSet"]]))
{
  check_help(args = args, help = help, func_name = "buildVarSet")
  required=c("gdb","varSetName","unitTable","unitName","output")
  expected=c("buildVarSet","gdb","output","varSetName","unitTable","unitName","intersection","weightName","where")
  
  check_required_args(required = required, args = args, func_name = "buildVarSet")
  check_expected_args(expected = expected, args = args, func_name = "buildVarSet")
  
  gdb <- rvat::gdb(args[["gdb"]])
  if (is.null(args[["intersection"]])){intersection <- c()} else{intersection <- args[["intersection"]]}
  if (is.null(args[["where"]])) {where <- c()} else{where <- args[["where"]]}
  if (is.null(args[["weightName"]])) {weightName <- 1} else{weightName <- args[["weightName"]]}
  
  rvat::buildVarSet(
                  object = gdb,
                  output = args[["output"]],
                  varSetName = args[["varSetName"]],
                  unitTable = args[["unitTable"]],
                  unitName = args[["unitName"]],
                  intersection = intersection,
                  where = where,
                  weightName = weightName)
}

# spatialClust ------------------------------------------------------------------
if (!is.null(args[["spatialClust"]]))
{
  check_help(args = args, help = help, func_name = "spatialClust")
  required=c("gdb","varSetName","unitTable","unitName","output", "windowSize", "overlap")
  expected=c("spatialClust","gdb","output","varSetName","unitTable","unitName","intersection","weightName","where","windowSize","overlap", "posField", "minTry")
  
  check_required_args(required = required, args = args, func_name = "spatialClust")
  check_expected_args(expected = expected, args = args, func_name = "spatialClust")
  
  gdb <- rvat::gdb(args[["gdb"]])
  if (is.null(args[["intersection"]])){intersection <- c()} else{intersection <- args[["intersection"]]}
  if (is.null(args[["where"]])) {where <- c()} else{where <- args[["where"]]}
  if (is.null(args[["weightName"]])) {weightName <- 1} else{weightName <- args[["weightName"]]}
  windowSize <- as.numeric(unlist(strsplit(args[["windowSize"]], split=",")))
  overlap <- as.numeric(unlist(strsplit(args[["overlap"]], split=",")))
  
  rvat::spatialClust(
    object = gdb,
    output = args[["output"]],
    varSetName = args[["varSetName"]],
    unitTable = args[["unitTable"]],
    unitName = args[["unitName"]],
    windowSize = windowSize,
    overlap = overlap,
    intersection = intersection,
    where = args[["where"]],
    weightName = weightName,
    posField = if(!is.null(args[["posField"]])) args[["posField"]] else "POS",
    minTry = if(!is.null(args[["minTry"]])) as.numeric(args[["minTry"]]) else 5
  )
}

# assocTest --------------------------------------------------------------------
if (!is.null(args[["assocTest"]]) && is.null(args[["aggregateFile"]]))
{
  check_help(args = args, help = help, func_name = "assocTest")
  
  required=c("gdb","pheno", "test")
  required_one_of=c("output", "outputResampling")
  expected=c("assocTest", "gdb","varSet", "VAR_id", "cohort","pheno", "keep", "test", "output",
             "name", "continuous", "singlevar", "covar", "offset", "geneticModel",
             "imputeMethod", "MAFweights", "maxitFirth","checkPloidy", "append", "skatResampling",
             "resamplingFile", "nResampling", "outputResampling", "methodResampling",
             "memlimitResampling",
             "minCallrateVar", "maxCallrateVar", "minCallrateSM", "maxCallrateSM",
             "minMAF", "maxMAF", "minMAC", "maxMAC", "minCarriers", "maxCarriers", "minCarrierFreq", "maxCarrierFreq", "memlimit", "seed")
  check_required_args(required = required, required_one_of = required_one_of, args = args, func_name = "assocTest")
  check_expected_args(expected = expected, args = args, func_name = "assocTest")
  # 
  gdb <- rvat::gdb(args[["gdb"]])
  if(!is.null(args[["varSet"]])) varSet <- rvat::varSetFile(args[["varSet"]]) else varSet <- NULL
  if(!is.null(args[["VAR_id"]])) VAR_id <- readLines(args[["VAR_id"]]) else VAR_id <- NULL
  if(!is.null(args[["keep"]])) keep <- readLines(args[["keep"]]) else keep <- NULL
  if(!is.null(args[["resamplingFile"]])) resamplingFile <- rvat::resamplingFile(args[["resamplingFile"]]) else resamplingFile <- NULL
  if(!is.null(args[["seed"]])) set.seed(as.numeric(args[["seed"]]))
  
  ## Run 
  rvat::assocTest(
    object = gdb, 
    varSet = varSet, 
    VAR_id = VAR_id,
    cohort = if(is.null(args[["cohort"]])) "SM" else args[["cohort"]], 
    output = args[["output"]],
    pheno = unlist(strsplit(args[["pheno"]], split=",")), 
    keep = keep,
    test = unlist(strsplit(args[["test"]],split=",")), 
    name = if (is.null(args[["name"]])) "none" else args[["name"]],
    continuous = if(!is.null(args[["continuous"]])) TRUE else FALSE,
    singlevar = if(!is.null(args[["singlevar"]])) TRUE else FALSE,
    covar = if(length(args[["covar"]]) > 0) lapply(unlist(strsplit(args[["covar"]],split="/")), 
                                                   FUN = function(x) unlist(strsplit(x, split=","))) else NULL,
    offset = args[["offset"]],
    geneticModel = if(is.null(args[["geneticModel"]])) "allelic" else unlist(strsplit(args[["geneticModel"]],split=",")),
    imputeMethod = if (is.null(args[["imputeMethod"]])) NULL else args[["imputeMethod"]],
    MAFweights = if (is.null(args[["MAFweights"]])) "none" else unlist(strsplit(args[["MAFweights"]],split=",")),
    checkPloidy = if(!is.null(args[["checkPloidy"]])) args[["checkPloidy"]] else NULL,
    maxitFirth = if (is.null(args[["maxitFirth"]])) 1000 else as.numeric(args[["maxitFirth"]]),
    resamplingFile = resamplingFile,
    nResampling = if (is.null(args[["nResampling"]])) 1000 else args[["nResampling"]],
    outputResampling = if (is.null(args[["outputResampling"]])) NULL else args[["outputResampling"]],
    methodResampling = if (is.null(args[["methodResampling"]])) NULL else args[["methodResampling"]],
    minCallrateVar = if (is.null(args[["minCallrateVar"]])) 0 else as.numeric(args[["minCallrateVar"]]),
    maxCallrateVar = if (is.null(args[["maxCallrateVar"]])) Inf else as.numeric(args[["maxCallrateVar"]]),
    minCallrateSM = if (is.null(args[["minCallrateSM"]])) 0 else as.numeric(args[["minCallrateSM"]]),
    maxCallrateSM = if (is.null(args[["maxCallrateSM"]])) Inf else as.numeric(args[["maxCallrateSM"]]),
    minMAF = if (is.null(args[["minMAF"]])) 0 else as.numeric(args[["minMAF"]]),
    maxMAF = if (is.null(args[["maxMAF"]])) 1 else as.numeric(args[["maxMAF"]]),
    minMAC = if (is.null(args[["minMAC"]])) 0 else as.numeric(args[["minMAC"]]),
    maxMAC = if (is.null(args[["maxMAC"]])) Inf else as.numeric(args[["maxMAC"]]),
    minCarriers = if (is.null(args[["minCarriers"]])) 0 else as.numeric(args[["minCarriers"]]),
    maxCarriers = if (is.null(args[["maxCarriers"]])) Inf else as.numeric(args[["maxCarriers"]]),
    minCarrierFreq = if (is.null(args[["minCarrierFreq"]])) 0 else as.numeric(args[["minCarrierFreq"]]),
    maxCarrierFreq = if (is.null(args[["maxCarrierFreq"]])) Inf else as.numeric(args[["maxCarrierFreq"]]),
    memlimit = if (is.null(args[["memlimit"]])) Inf else as.numeric(args[["memlimit"]])
  )
}

# summariseGeno -----------------------------------------------------------
if (!is.null(args[["summariseGeno"]]))
{
  check_help(args = args, help = help, func_name = "summariseGeno")
  
  required=c("gdb", "output")
  expected=c("summariseGeno", "gdb", "varSet", "VAR_id","cohort","pheno", "memlimit", 
             "keep", "output","geneticModel", "checkPloidy",
             "splitBy",
             "minCallrateVar", "maxCallrateVar", "minCallrateSM", "maxCallrateSM",
             "minMAF", "maxMAF", "minMAC", "maxMAC", "minCarriers", "maxCarriers", "minCarrierFreq", "maxCarrierFreq")
  check_required_args(required = required, args = args, func_name = "summariseGeno")
  check_expected_args(expected = expected, args = args, func_name = "summariseGeno")
  
  gdb <- rvat::gdb(args[["gdb"]])
  varSet <- if(!is.null(args[["varSet"]])) rvat::varSetFile(args[["varSet"]]) else NULL
  VAR_id <- if(!is.null(args[["VAR_id"]])) readLines(args[["VAR_id"]]) else NULL
  if(!is.null(args[["keep"]])) keep <- readLines(args[["keep"]]) else keep <- NULL
  
  rvat::summariseGeno(
    object = gdb, 
    cohort = args[["cohort"]], 
    varSet = varSet,
    VAR_id = VAR_id,
    pheno = if(!is.null(args[["pheno"]])) args[["pheno"]] else NULL, 
    memlimit = if(!is.null(args[["memlimit"]])) as.numeric(args[["memlimit"]]) else 1000, 
    output = args[["output"]],
    keep = keep,
    splitBy = if(!is.null(args[["splitBy"]])) args[["splitBy"]] else NULL, 
    geneticModel = if(is.null(args[["geneticModel"]])) "allelic" else unlist(strsplit(args[["geneticModel"]],split=",")),
    checkPloidy = if(!is.null(args[["checkPloidy"]])) args[["checkPloidy"]] else NULL,
    minCallrateVar = if (is.null(args[["minCallrateVar"]])) 0 else as.numeric(args[["minCallrateVar"]]),
    maxCallrateVar = if (is.null(args[["maxCallrateVar"]])) Inf else as.numeric(args[["maxCallrateVar"]]),
    minCallrateSM = if (is.null(args[["minCallrateSM"]])) 0 else as.numeric(args[["minCallrateSM"]]),
    maxCallrateSM = if (is.null(args[["maxCallrateSM"]])) Inf else as.numeric(args[["maxCallrateSM"]]),
    minMAF = if (is.null(args[["minMAF"]])) 0 else as.numeric(args[["minMAF"]]),
    maxMAF = if (is.null(args[["maxMAF"]])) 1 else as.numeric(args[["maxMAF"]]),
    minMAC = if (is.null(args[["minMAC"]])) 0 else as.numeric(args[["minMAC"]]),
    maxMAC = if (is.null(args[["maxMAC"]])) Inf else as.numeric(args[["maxMAC"]]),
    minCarriers = if (is.null(args[["minCarriers"]])) 0 else as.numeric(args[["minCarriers"]]),
    maxCarriers = if (is.null(args[["maxCarriers"]])) Inf else as.numeric(args[["maxCarriers"]]),
    minCarrierFreq = if (is.null(args[["minCarrierFreq"]])) 0 else as.numeric(args[["minCarrierFreq"]]),
    maxCarrierFreq = if (is.null(args[["maxCarrierFreq"]])) Inf else as.numeric(args[["maxCarrierFreq"]])
  )
}

# aggregate --------------------------------------------------------------------

if (!is.null(args[["aggregate"]]))
{
  check_help(args = args, help = help, func_name = "aggregate")
  
  required=c("gdb", "output")
  expected=c("aggregate", "gdb", "varSet", "VAR_id","cohort","pheno", "memlimit", 
             "keep", "output","geneticModel","imputeMethod", "MAFweights", "checkPloidy", "signif",
             "minCallrateVar", "maxCallrateVar", "minCallrateSM", "maxCallrateSM",
             "minMAF", "maxMAF", "minMAC", "maxMAC", "minCarriers", "maxCarriers", "minCarrierFreq", "maxCarrierFreq")
  check_required_args(required = required, args = args, func_name = "aggregate")
  check_expected_args(expected = expected, args = args, func_name = "aggregate")
  
  gdb <- rvat::gdb(args[["gdb"]])
  varSet <- if(!is.null(args[["varSet"]])) rvat::varSetFile(args[["varSet"]]) else NULL
  VAR_id <- if(!is.null(args[["VAR_id"]])) readLines(args[["VAR_id"]]) else NULL
  if(!is.null(args[["keep"]])) keep <- readLines(args[["keep"]]) else keep <- NULL
  
  rvat::aggregate(
    x = gdb, 
    varSet = varSet,
    VAR_id = VAR_id,
    cohort = if(is.null(args[["cohort"]])) "SM" else args[["cohort"]], 
    pheno = if(!is.null(args[["pheno"]])) args[["pheno"]] else NULL, 
    memlimit = if(!is.null(args[["memlimit"]])) as.numeric(args[["memlimit"]]) else 1000, 
    output = args[["output"]],
    keep = keep,
    signif = if(!is.null(args[["signif"]])) as.numeric(args[["signif"]]) else 3, 
    geneticModel = if(is.null(args[["geneticModel"]])) "allelic" else unlist(strsplit(args[["geneticModel"]],split=",")),
    imputeMethod = if (is.null(args[["imputeMethod"]])) "meanImpute" else args[["imputeMethod"]],
    MAFweights = if (is.null(args[["MAFweights"]])) "none" else unlist(strsplit(args[["MAFweights"]],split=",")),
    checkPloidy = if(!is.null(args[["checkPloidy"]])) args[["checkPloidy"]] else NULL,
    minCallrateVar = if (is.null(args[["minCallrateVar"]])) 0 else as.numeric(args[["minCallrateVar"]]),
    maxCallrateVar = if (is.null(args[["maxCallrateVar"]])) Inf else as.numeric(args[["maxCallrateVar"]]),
    minCallrateSM = if (is.null(args[["minCallrateSM"]])) 0 else as.numeric(args[["minCallrateSM"]]),
    maxCallrateSM = if (is.null(args[["maxCallrateSM"]])) Inf else as.numeric(args[["maxCallrateSM"]]),
    minMAF = if (is.null(args[["minMAF"]])) 0 else as.numeric(args[["minMAF"]]),
    maxMAF = if (is.null(args[["maxMAF"]])) 1 else as.numeric(args[["maxMAF"]]),
    minMAC = if (is.null(args[["minMAC"]])) 0 else as.numeric(args[["minMAC"]]),
    maxMAC = if (is.null(args[["maxMAC"]])) Inf else as.numeric(args[["maxMAC"]]),
    minCarriers = if (is.null(args[["minCarriers"]])) 0 else as.numeric(args[["minCarriers"]]),
    maxCarriers = if (is.null(args[["maxCarriers"]])) Inf else as.numeric(args[["maxCarriers"]]),
    minCarrierFreq = if (is.null(args[["minCarrierFreq"]])) 0 else as.numeric(args[["minCarrierFreq"]]),
    maxCarrierFreq = if (is.null(args[["maxCarrierFreq"]])) Inf else as.numeric(args[["maxCarrierFreq"]])
  )
}

# mergeAggregateFiles ----------------------------------------------------------

if (!is.null(args[["mergeAggregateFiles"]]))
{
  check_help(args = args, help = help, func_name = "mergeAggregateFiles")
  
  required=c("filelist", "output")
  expected=c("mergeAggregateFiles", "filelist", "collapse", "output", "verbose", "checkDups")
  check_required_args(required = required, args = args, func_name = "mergeAggregateFiles")
  check_expected_args(expected = expected, args = args, func_name = "mergeAggregateFiles")
  
  files <- readLines(args[["filelist"]])
  checkDups <- if(is.null(args[["checkDups"]])) TRUE else as.logical(args[["checkDups"]])
  aggregateFileList <- rvat::aggregateFileList(files, checkDups = checkDups)
  
  rvat::mergeAggregateFiles(
    object = aggregateFileList,
    collapse = if(!is.null(args[["collapse"]])) as.logical(args[["collapse"]]) else TRUE,
    output = args[["output"]],
    verbose = if(!is.null(args[["verbose"]])) TRUE else FALSE
  )
}

# assocTest-aggregateFile ------------------------------------------------------

if (!is.null(args[["assocTest"]]) && !is.null(args[["aggregateFile"]]))
{
  check_help(args = args, help = help, func_name = "aggregateGeneSetAssoc")
  
  required=c("geneSetFile", "aggregateFile", "gdb","cohort","pheno", "test", "output")
  expected=c("assocTest", "geneSetFile", "gdb" , "cohort", "pheno", "test",
             "aggregateFile", "varSetFile", "name", "continuous", "covar", "substractCovar", "keep",
             "output")
  check_required_args(required = required, args = args, func_name = "aggregateGeneSetAssoc")
  check_expected_args(expected = expected, args = args, func_name = "aggregateGeneSetAssoc")
  # 
  gdb <- rvat::gdb(args[["gdb"]])
  genesetfile <- rvat::geneSetFile(args[["geneSetFile"]])
  aggregatefile <- rvat::aggregateFile(args[["aggregateFile"]])
  
  if(!is.null(args[["keep"]])) keep <- readLines(args[["keep"]]) else keep <- NULL
  
  rvat::assocTest(
    object = aggregatefile,
    gdb = gdb, 
    cohort = args[["cohort"]], 
    pheno = unlist(strsplit(args[["pheno"]], split=",")), 
    test = unlist(strsplit(args[["test"]],split=",")), 
    geneSet = genesetfile,
    name = if (is.null(args[["name"]])) "none" else args[["name"]], 
    continuous =  if(!is.null(args[["continuous"]])) TRUE else FALSE,
    covar = if(length(args[["covar"]]) > 0) lapply(unlist(strsplit(args[["covar"]],split="/")), 
                                                   FUN = function(x) unlist(strsplit(x, split=","))) else NULL, 
    substractCovar = if(!is.null(args[["substractCovar"]])) args[["substractCovar"]] else NULL,
    maxitFirth = if (is.null(args[["maxitFirth"]])) 1000 else as.numeric(args[["maxitFirth"]]),
    keep = keep,
    output = args[["output"]]
  )
}


# geneSetAssoc --------------------------------------------------------
if (!is.null(args[["geneSetAssoc"]]))
{
  check_help(args = args, help = help, func_name = "geneSetAssoc")
  
  required=c("geneSetFile", "output")
  required_one_of=c("rvbResult", "nullmodel")
  expected=c("geneSetAssoc", "nullmodel", "rvbResult", "geneSetFile", "cormatrix", "covar", "test",
             "threshold", "Zcutoffs", "minSetSize", "maxSetSize", "oneSided", "memlimit",
             "REML", "ID", "output")
  check_required_args(required = required, required_one_of = required_one_of, args = args, func_name = "geneSetAssoc")
  check_expected_args(expected = expected, args = args, func_name = "geneSetAssoc")
   
  genesetfile <- rvat::geneSetFile(args[["geneSetFile"]])
  if(!is.null(args[["cormatrix"]])) cormatrix <- readRDS(args[["cormatrix"]]) else cormatrix <- NULL
  
  if(is.null(args[["nullmodel"]])) {
    result <- rvat::rvbResult(args[["rvbResult"]])
    rvat::geneSetAssoc(
      object = result,
      geneSetList = NULL,
      geneSetFile = genesetfile,
      cormatrix = cormatrix,
      covar = if(!is.null(args[["covar"]])) unlist(strsplit(args[["covar"]], split=",")) else NULL,
      test = unlist(strsplit(args[["test"]],split=",")), 
      threshold = if(!is.null(args[["threshold"]])) as.numeric(args[["threshold"]]) else NULL,
      Zcutoffs = if(!is.null(args[["Zcutoffs"]])) as.numeric(unlist(strsplit(args[["Zcutoffs"]], split=","))) else NULL,
      minSetSize = if(!is.null(args[["minSetSize"]])) as.numeric(args[["minSetSize"]]) else 1,
      maxSetSize = if(!is.null(args[["maxSetSize"]])) as.numeric(args[["maxSetSize"]]) else Inf,
      oneSided =  if(!is.null(args[["oneSided"]])) as.logical(args[["oneSided"]]) else TRUE,
      memlimit = if(!is.null(args[["memlimit"]])) as.numeric(args[["memlimit"]]) else 1000, 
      REML = if(!is.null(args[["REML"]])) as.logical(args[["REML"]]) else TRUE,
      ID = if(!is.null(args[["ID"]])) args[["ID"]] else "unit",
      output = args[["output"]]
    )
  } else {
    nullmodel <- readRDS(args[["nullmodel"]])
    
    rvat::geneSetAssoc(
      object = nullmodel,
      geneSetList = NULL,
      geneSetFile = genesetfile,
      test = unlist(strsplit(args[["test"]],split=",")), 
      threshold = if(!is.null(args[["threshold"]])) as.numeric(args[["threshold"]]) else NULL,
      minSetSize = if(!is.null(args[["minSetSize"]])) as.numeric(args[["minSetSize"]]) else 1,
      maxSetSize = if(!is.null(args[["maxSetSize"]])) as.numeric(args[["maxSetSize"]]) else Inf,
      oneSided =  if(!is.null(args[["oneSided"]])) as.logical(args[["oneSided"]]) else TRUE,
      memlimit = if(!is.null(args[["memlimit"]])) as.numeric(args[["memlimit"]]) else 1000, 
      ID = if(!is.null(args[["ID"]])) args[["ID"]] else "unit",
      output = args[["output"]]
    )
  }
}


# buildResamplingFile ------------------------------------------------------------------

if (!is.null(args[["buildResamplingFile"]]))
{
  check_help(args = args, help = help, func_name = "buildResamplingFile")
  required=c("nSamples", "nResampling", "output")
  expected=c("buildResamplingFile", "nSamples", "nResampling", "output", "methodResampling", "memlimit", "seed")
  
  check_required_args(required = required, args = args, func_name = "buildResamplingFile")
  check_expected_args(expected = expected, args = args, func_name = "buildResamplingFile")
  
  if(is.null(args[["seed"]])) set.seed(10) else set.seed(as.numeric(args[["seed"]]))
  
  rvat::buildResamplingFile(
    nSamples = as.numeric(args[["nSamples"]]),
    nResampling = as.numeric(args[["nResampling"]]),
    memlimit = if(is.null(args[["memlimit"]])) 1000 else as.numeric(args[["memlimit"]]),
    methodResampling = if(is.null(args[["methodResampling"]])) "permutation" else args[["methodResampling"]],
    output = args[["output"]]
  )
}

if (!is.null(args[["vcfInfo2Table"]]))
{
  check_help(args = args, help = help, func_name = "vcfInfo2Table")
  required=c("vcf", "output")
  expected=c("vcfInfo2Table", "vcf", "output", "splitMultiallelic")
  
  check_required_args(required = required, args = args, func_name = "vcfInfo2Table")
  check_expected_args(expected = expected, args = args, func_name = "vcfInfo2Table")
  
  rvat::vcfInfo2Table(
    vcf = args[["vcf"]],
    output = args[["output"]],
    splitMultiallelic = if(is.null(args[["splitMultiallelic"]])) TRUE else as.logical(args[["splitMultiallelic"]])
  )
} 

message("\nFinished!")

if(length(warnings()) > 0) {
  message("Overview of warnings:")
  summary(warnings())
}

message(sprintf("End time: %s", Sys.time()))

