rvat_cli_help  <- list(
  buildGdb = "
buildGdb

  Creates a new gdb file. The gdb can be structured and populated using a provided vcf file. 

Usage:
  Rscript rvat.R --buildGdb --vcf={vcf} --output={output} [options]

Arguments:
  output           Path for output gdb file
  vcf              Input vcf file used to structure and populate gdb. Warning this function makes the following of assumptions: 1) strict adherence to vcf format (GT subfield first element in genotype firelds), 2) multiallelic records have been split, 3) desired genotype QC has already been applied (DP,GQ filters), 4) GT values conform to the set {0/0,0/1,1/0,1/1,./.,0|0,0|1,1|0,1|1,.|.}. Multiallelic parsing and genotype QC can be performed using vcftools and/or accompanying parser scripts included on the rvat github.
  skipIndexes      Flag to skip generation of indexes for var and dosage table (VAR_id;CHROM, POS,REF,ALT). Typically only required if you plan to use gdbConcat to concatenate a series of separately generated gdb files before use
  skipVarRanges    Flag to skip generation of ranged var table. Typically only required if you plan to use gdbConcat to concatenate a series of separately generated gdb files before use
  overWrite        Overwrite if output already exists? Defaults to FALSE, in which case an error is raised.
  genomeBuild      Optional genome build to include in the gdb metadata. If specified, it will be used to set ploidies (diploid, XnonPAR, YnonPAR) if the genome build is implemented in RVAT (currently: GRCh37, hg19, GRCh38, hg38).
  memlimit         Maximum number of vcf records to parse at a time, defaults to 1000.
  ",
  
  concatGdb = "
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
  
  subsetGdb="
subsetGdb

  Function to allow for generation of a child gdb from a parent gdb, with the option to filter retained variants through table intersection operations and SQL where statements.

Usage:
  Rscript rvat.R --subsetGdb --gdb={gdb} --output={output} [options]

Arguments:
  gdb            gdb file path
  output         Output file name (output will be a new gdb file).
  intersection   Additional tables to filter through intersection (ie variants absent from intersection tables will not appear in output). Multiple tables should be ',' delimited.
  where          An SQL compliant where clause to filter output; eg: \"CHROM=2 AND POS between 5000 AND 50000 AND AF<0.01 AND (cadd.caddPhred>15 OR snpEff.SIFT='D')\"
  VAR_id         Retain only variants with matching VAR_id.
  tables         Tables to retain from the gdb, multiple tables, Multiple tables should be ',' delimited. By default all tables will be included in the output gdb.
  skipIndexes    Flag to skip generation of indexes for var and dosage table (VAR_id;CHROM, POS,REF,ALT). Typically only required if you plan to use gdbConcat to concatenate a series of separately generated gdb files before use.
  overWrite      Flag indicating whether `output` should be overwritten if it already exists.
",
  
  
  uploadAnno = "
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
  sep            Field delimiter. Defaults to \t.
  skipRemap      Flag indicating whether to skip mapping of records to VAR_id using CHROM,POS,REF,ALT
  skipIndexes    Flag indicating whether to skip indexing of imported table.
  ignoreAlleles  Flag indicating whether to consider REF and ALT allele during mapping of records to VAR_id or just CHROM,POS.
  keepUnmapped   Flag indicating whether to keep records which cannot be mapped to the gdb.
  mapRef         Name of lookup table for VAR_id assignment. Defaults to 'var'.
",
  
  uploadCohort="
uploadCohort

  Function to upload cohort data tables to gdb. These will automatically be reformatted and sorted to match the ordering of samples in the gdb genotype records.

Usage:
  Rscript rvat.R --uploadCohort --gdb={gdb} --name={name} --value={value} [options]

Arguments:
  gdb            gdb file path
  name           Name to assign to cohort
  value          Full path to cohort annotation file. Must contain an 'IID' column matching to SM table and a 'sex' column (0=missing,1=male,2=female)
  sep            Field delimiter, defaults to '\t'.

",
  
  listAnno="
listAnno

  List variant info tables that have been uploaded to the gdb.

Usage:
  Rscript rvat.R --listAnno --gdb={gdb}

Arguments:
  gdb            gdb file path

",
  listCohort="
listCohort
  
  List sample cohort tables that have been uploaded to the gdb.

Usage:
  Rscript rvat.R --listCohort --gdb={gdb}

Arguments:
  gdb           gdb file path

",
  dropTable="
dropTable

  Drop table from gdb and clear from annotation / cohort metadata tables.

Usage:
  Rscript rvat.R --dropTable --gdb={gdb} --name={name}

Arguments:
  gdb           gdb file path
  name          Name of table to drop

",
  
  mapVariants="
mapVariants

  Method to map the variants in a gdb to a set of ranges or features. 
  The input can be a set of ranges (CHROM, start, end), a bed-file or a gff/gtf-file. 
  Variants in the gdb will be mapped onto those ranges and annotated with the features/columns 
  included in the input file. 
  For example, variants can be easily mapped upon genomic features downloaded in gff format from ensembl. 
  The output can be written to disk  (--output flag) or directly uploaded to the gdb (--uploadName flag). 

Usage:
  Rscript rvat.R --mapVariants --gdb={gdb} --gff={gff} --uploadName={uploadName}

Arguments:
  gdb            gdb file path
  ranges         A filepath to a ranges file containing at least 'CHROM','start', and 'end' columns.
                 Separator can be specified using the `sep` parameter (defaults to '\t').
  gff            Path to a gff- or gtf-file. 
  bed            Path to a bed-file. Specify extra columns using the --bedCols flag
  bedCols        A comma-delimited list of names of the extra columns to read from the BED-file. 
  fields         A comma-delimited list of feature fields to keep. By default all fields are kept.
  uploadName     Name of table to upload to the gdb. If not specified, specifiy --output to write results to disk.
  output         Optionally, an output file path. Can be used instead of --uploadName to write the results to disk.
  sep            Field separator, relevant if --ranges is specified. Defaults to '\t'. 
  skipIndexes    skipIndexes Flag indicating whether to skip indexing of imported table. 
                 Relevant if --uploadName is specified, and thus the output table is imported in the gdb.
                 By default the table is indexed on VAR_id.
  overWrite      If --uploadName is specified, should an existing table in the gdb with the same name be overwitten?
                 By default the method is aborted when the table already exists in the gdb.
",
  buildVarSet="
buildVarSet

  Generate optionally weighted variant sets using annotation table(s) uploaded to the gdb.

Usage:
  Rscript rvat.R --buildVarSet --gdb={gdb} --unitTable={unitTable} --unitName={unitName} --output={output} [options]

Arguments:
  gdb           gdb file path
  varSetName    Name to assign varSet grouping. This identifier column is used to allow for subsequent mergeing of multiple varSet files for coordinated analysis of multiple variant filtering/ weighting strategies)
  unitTable     Table containing aggregation unit mappings
  unitName      Field to utilize for aggregation unit names
  intersection  Additional tables to filter through intersection (i.e. variants absent from intersection tables will not appear in output). Multiple tables should be ',' delimited
  where         An SQL compliant where clause to filter output; eg: \"CHROM=2 AND POS between 5000 AND 50000 AND AF<0.01 AND (cadd.caddPhred>15 OR snpEff.SIFT='D')\".
  weightName    Field name for desired variant weighting, must be a column within unitTable or other intersection table. Default value of 1 is equivalent to no weighting.
  output        Output file name (output will be gz compressed text)
",
  spatialClust="
spatialClust

  Generate weighted variant sets for use in association testing, with partitioning by genomic distances as described (Fier, GenetEpidemiol, 2017).

Usage:
  Rscript rvat.R --spatialClust --gdb={gdb} --output={output} --unitTable={unitTable} --unitName={unitName} --windowSize={windowSize} --overlap={overlap} [options]
  
Arguments:
  gdb           gdb file path
  output        Output file name (output will be gz compressed text)
  varSetName    Name to assign varSet grouping. This identifier column is used to allow for subsequent mergeing of multiple varSet files for coordinated analysis of multiple variant filtering/ weighting strategies)
  unitTable     Table containing aggregation unit mappings
  unitName      Field to utilize for aggregation unit names
  windowSize    Starting fixed window sizes (number of variants), comma-delimited.
  overlap       Starting fixed window overlap (number of variants, length must match windowSize), comma-delimited.
  intersection  Additional tables to filter through intersection (i.e. variants absent from intersection tables will not appear in output). Multiple tables should be ',' delimited
  where         An SQL compliant where clause to filter output; eg: \"CHROM=2 AND POS between 5000 AND 50000 AND AF<0.01 AND (cadd.caddPhred>15 OR snpEff.SIFT='D')\".
  weightName    Field name for desired variant weighting, must be a column within unitTable or other intersection table. Default value of 1 is equivalent to no weighting.
  posField      Column name to take as variants position. Default is 'POS' which typically corresponds to genomics position. Can be reset to use CDS or other coordinates. 'HGVSc' is a recognized identifier and CDS coordinates will be extracted automatically.
  minTry        Minimum number of variants in varset to perform clustering on. If number of variants < minTry, all variants will be returned as a single cluster. Default = 5.
 
",
  summariseGeno="
summariseGeno

  Returns a per variant summary of genotype counts, frequencies, call-rates and hwe testing.

Usage:
  Rscript rvat.R --summariseGeno --gdb={gdb} --VAR_id={VAR_id} --cohort={cohort} --output={output} [options]

Arguments:
  gdb             gdb file path
  cohort          If a valid cohort name is provided, then the uploaded data for this cohort is used to filter and annotate the genotypes If not specified, all samples in the gdb will be loaded.
  varSet          varSetFile path. Alternatively the --VAR_id flag can be specified.
  VAR_id          A list of VAR_ids, alternatively the varSet parameter can be specified.
                  The `memlimit` argument controls how many variants to analyze at a time.
  pheno           Cohort field with response variable, although not used within this method, this can be useful to filter samples which have missing data for the response variable.
  memlimit        Maximum number of variants to load at once (if --VAR_id is specified).
  geneticModel    Genetic model to apply ('allelic', 'recessive', 'dominant'). Defaults to 'allelic'.
  checkPloidy     Version of the human genome to use when assigning variant ploidy (diploid, XnonPAR, YnonPAR). Accepted inputs are GRCh37, hg19, GRCh38, hg38.
                  If not specified, the genome build in the gdb will be used, if available (included if the `genomeBuild` parameter was set in --buildGdb).
                  Otherwise, if the genome build is not included in the gdb metadata, and no value is provided, then all variants are assigned the default ploidy of 'diploid'.
  keep            Filepath to list of samples to retain in analysis (by default all samples are kept).
  output          Output file path for results.
  splitBy         Split variant summaries by labels indicated in the specified field.
  minCallrateVar  Minimum genotype rate for variant retention.
  maxCallrateVar  Maximum genotype rate for variant retention.
  minCallrateSM   Minimum genotype rate for sample retention.
  maxCallrateSM   Maximum genotype rate for sample retention.
  minMAF          Minimum minor allele frequency for variant retention.
  maxMAF          Maximum minor allele frequency for variant retention.
  minMAC          Minimum minor allele count for variant retention.
  maxMAC          Maximum minor allele count for variant retention.
  minCarriers     Minimum carrier count for variant retention.
  maxCarriers     Maximum carrier count for variant retention.
  minCarrierFreq  Minimum carrier frequency for variant retention.
  maxCarrierFreq  Maximum carrier frequency for variant retention.
  not-strict      Flag to turn off strict checks. Strict checks currently includes checking whether supplied varSetFile/varSetList was generated from the same gdb as the gdb provided in --gdb.
  ",
  
  aggregate="
aggregate

  Returns an aggregate of genotypes for each individual. 
  Specified genetic model, weights, MAF-weighting are taken into account when aggregating.
  Aggregates are written to disk in the aggregateFile format, which can be used as input
  for assocTest to perform gene set burden analyses.

Usage:
  Rscript rvat.R --aggregate --gdb={gdb} --varSet={varSet} --cohort={cohort} --output={output} [options]

Arguments:
  gdb             gdb file path
  cohort          cohort data previously uploaded to the gdb (see uploadCohort).
  varSet          varSetFile path. Alternatively the --VAR_id flag can be specified.
  VAR_id          A list of VAR_ids, alternatively the varSet parameter can be specified.
                  If single variant tests are ran, the --memlimit parameter controls how many variants to analyze at a time.
  pheno           Cohort field with response variable, although not used within this method, this can be useful to filter samples which have missing data for the response variable.
  memlimit        Maximum number of variants to load at once (if --VAR_id is specified).
  geneticModel    Genetic model to apply ('allelic', 'recessive', 'dominant'). Defaults to 'allelic'.
  imputeMethod    Which imputation method to apply? ('meanImpute' or 'missingToRef'). Defaults to 'meanImpute'.
  MAFweights      MAF weighting method. Currently Madsen-Browning ('mb') is implemented, by default no MAF weighting is applied.
  checkPloidy     Version of the human genome to use when assigning variant ploidy (diploid, XnonPAR, YnonPAR). Accepted inputs are GRCh37, hg19, GRCh38, hg38.
                  If not specified, the genome build in the gdb will be used, if available (included if the `genomeBuild` parameter was set in --buildGdb).
                  Otherwise, if the genome build is not included in the gdb metadata, and no value is provided, then all variants are assigned the default ploidy of 'diploid'.
  keep            Filepath to list of samples to retain in analysis (by default all samples are kept).
  output          Output file path for results.
  signif          Number of significant digits to store. Defaults to 6.
  minCallrateVar  Minimum genotype rate for variant retention.
  maxCallrateVar  Maximum genotype rate for variant retention.
  minCallrateSM   Minimum genotype rate for sample retention.
  maxCallrateSM   Maximum genotype rate for sample retention.
  minMAF          Minimum minor allele frequency for variant retention.
  maxMAF          Maximum minor allele frequency for variant retention.
  minMAC          Minimum minor allele count for variant retention.
  maxMAC          Maximum minor allele count for variant retention.
  minCarriers     Minimum carrier count for variant retention.
  maxCarriers     Maximum carrier count for variant retention.
  minCarrierFreq  Minimum carrier frequency for variant retention.
  maxCarrierFreq  Maximum carrier frequency for variant retention.
  not-strict      Flag to turn off strict checks. Strict checks currently includes checking whether supplied varSetFile/varSetList was generated from the same gdb as the gdb provided in --gdb.
  ",
  mergeAggregateFiles = "
mergeAggregateFiles

  Merge aggregrateFiles, this will generate a new aggregateFile including all aggregates across provided aggregateFiles.

Usage:
  Rscript rvat.R --mergeAggregateFiles --filelist={filelist} --output={output} [options]

Arguments:
  filelist       Filepath to a list of aggregateFiles.
  output         Output file name (output will be an aggregateFile). 
  not-checkDups  Flag to indicate that no error should be raised if unit names are duplicated across aggregateFiles.
",
  collapseAggregateFiles = "
collapseAggregateFiles

  Collapse aggregrateFiles by aggregating values across aggregateFiles. 
  This will result in one aggregate score for each sample, representing the aggregate value across aggregate files. 
  The output will be a two-column matrix including sample IDs and aggregate scores respectively.

Usage:
  Rscript rvat.R --collapseAggregateFiles --filelist={filelist} --output={output} [options]

Arguments:
  filelist       Filepath to a list of aggregateFiles.
  output         Output file name.
  not-checkDups  Flag to indicate that no error should be raised if unit names are duplicated across aggregateFiles.
",
  assocTest="
assocTest

  Main method for performing association tests on binary or quantitative traits. 
  Performs either aggregate tests or single variant tests, with a wide range of available statistical methods.
  Various parameters allow for implementation of for example recessive/dominant analyses, 
  MAF-weighted burden tests and permutation-based association tests. 
  Variant and sample filters (e.g. maximum MAF, minimum number of carriers) can be specified.
  
  assocTest can be run using a gdb or an aggregateFile object, both usages are outlined below.
  
Usage:
  Rscript rvat.R --assocTest --gdb={gdb} --varSet={varSet} --cohort={cohort} --pheno={pheno} --test={test} --output={output} [options]

Arguments:
  gdb                 gdb file path
  pheno               Cohort field to test as response variable, the response variable can either be binary (0/1) or continuous. 
                      If the response variable is continuous, include the --continuous flag. Multiple phenotypes can be specified using ',' to delimit them.
  test                Statistical tests to run (',' delimited), options include firth,glm,lm,nbinom,skat,skat_burden,skato,skat_robust,skato_robust,skat_burden_robust, acatv, acatvSPA, and acatvfirth.
  cohort              If a valid cohort name is provided, then the uploaded data for this cohort is used to filter and annotate the genoMatrix object. 
                      If not specified, all samples in the gdb will be loaded.
  varSet              varSetFile path.
  VAR_id              A list of VAR_ids, alternatively the --varSet parameter can be specified.
                      If single variant tests are ran, the --memlimit argument controls how many variants to analyze at a time.
  name                Optional name for the analysis, defaults to 'none'
  continuous          Flag that indicates that the phenotype(s) is/are continuous.
  singlevar           Flag that indicates that single variant tests should be run.
  covar               Covariates, ',' delimited. Multiple sets of covariates can be specified, delimited by '/'.
  geneticModel        Genetic models to test ('allelic', 'recessive', 'dominant'), multiple models can be specified (',' delimited).
  imputeMethod        Imputation method, either 'meanImpute' or 'missingToRef'. If not specified, single variant tests are not imputed, whereas burden tests are mean imputed.
  MAFweights          MAF weighting method. Currently Madsen-Browning ('mb') is implemented, by default no MAF weighting is applied. Multiple MAFweights can be specified (comma-delimited), in which case each will be analyzed separately.
  maxitFirth          Maximum number of iterations to use for estimating firth confidence intervals. Defaults to 1000.
  checkPloidy         Version of the human genome to use when assigning variant ploidy (diploid, XnonPAR, YnonPAR). Accepted inputs are GRCh37, hg19, GRCh38, hg38.
                      If not specified, the genome build in the gdb will be used, if available (included if the `genomeBuild` parameter was set in --buildGdb).
                      Otherwise, if the genome build is not included in the gdb metadata, and no value is provided, then all variants are assigned the default ploidy of 'diploid'.
  keep                Filepath to list of samples to retain in analysis (by default all samples are kept).
  output              Output file path for results.
  methodResampling    Which method to use for resampling? ('permutation' currently implemented). 
                      Defaults to `NULL`, in which case no resampling is performed. 
  resamplingFile      file path to a resamplingFile (see `buildResamplingFile` function)
  nResampling         Number of resamplings (note this is ignored if a resamplingFile is specified)
  outputResampling    If specified (a filepath), resampling results will be stored in the respective filepath.
  memlimitResampling  Maximum number of resamplings to perform at a time. 
                      Resampling generates a matrix of n x p, where n is the number of samples and p the number of resamplings thus, 
                      for large number of resamplings it can be more efficient to split the permutations in chunks of size memlimitResampling. 
  minCallrateVar      Minimum genotype rate for variant retention.
  maxCallrateVar      Maximum genotype rate for variant retention.
  minCallrateSM       Minimum genotype rate for sample retention.
  maxCallrateSM       Maximum genotype rate for sample retention.
  minMAF              Minimum minor allele frequency for variant retention.
  maxMAF              Maximum minor allele frequency for variant retention.
  minMAC              Minimum minor allele count for variant retention.
  maxMAC              Maximum minor allele count for variant retention.
  minCarriers         Minimum carrier count for variant retention.
  maxCarriers         Maximum carrier count for variant retention.
  minCarrierFreq      Minimum carrier frequency for variant retention.
  maxCarrierFreq      Maximum carrier frequency for variant retention.
  memlimit            Maximum number of variants to load at once (if --VAR_id is specified).
  seed                Seed to set when running permutation analyses.
  not-strict          Flag to turn off strict checks. Strict checks currently includes checking whether supplied varSetFile/varSetList was generated from the same gdb as the gdb provided in --gdb.
  

Usage:
  Rscript rvat.R --assocTest --aggregateFile={aggregateFile} --geneSet={geneSet} --cohort={cohort} --pheno={pheno} --test={test} --output={output} [options]

Arguments:
  aggregateFile       Path to an aggregateFile.
  pheno               Cohort field to test as response variable, the response variable can either be binary (0/1) or continuous. 
                      If the response variable is continuous, include the --continuous flag. Multiple phenotypes can be specified using ',' to delimit them.
  test                Statistical tests to run (',' delimited), options include firth,glm,lm,nbinom.
  geneSet             Path to a geneSetFile.
  gdb                 gdb file path.
  cohort              If a valid cohort name is provided, then the uploaded data for this cohort is used to filter and annotate the genoMatrix object. 
                      If not specified, all samples in the gdb will be loaded.
  name                Optional name for the analysis, defaults to 'none'
  continuous          Flag that indicates that the phenotype(s) is/are continuous.
  covar               Covariates, ',' delimited. Multiple sets of covariates can be specified, delimited by '/'.
  substractCovar      Covariate from which aggregate should be substracted. Useful when adjusting for total variant counts, by specifying the total variant count variable here, 
                      the aggregate score of the gene set tested will be substracted from the total count variable.
  dropUnits           Optional list of units to exclude.
  maxitFirth          Maximum number of iterations to use for estimating firth confidence intervals. Defaults to 1000.
  keep                Filepath to list of samples to retain in analysis (by default all samples are kept).
  output              Output file path for results.
  not-strict          Flag to turn off strict checks. Strict checks currently includes checking whether supplied varSetFile/varSetList was generated from the same gdb as the gdb provided in --gdb.
  ",
  
  buildResamplingFile="
buildResamplingFile

  Currently the 'permutation' approach is implemented. 
  The resampling matrix will be written to disk, and can be connected to with resamplingFile. 
  Can be used in combination with assocTest to perform resampling tests.

Usage:
  Rscript rvat.R --buildResamplingFile --nSamples={nSamples} --nResamplings={nResamplings} --output={output} [options]

Arguments:
  nSamples            Number of samples
  nResampling         Number of resamplings. Defaults to 1000.
  memlimit            Chunk sizes
  methodResampling    Resampling method, currently 'permutation' is implemented.
  output              File path (.gz extension) to write output to.
  seed                Set a seed for reproducibility.
",
  
  buildGeneSet = "
buildGeneSet

  Build a geneSetFile for use in gene set analyses (geneSetAssoc or assocTest-aggregateFile).
  Currently these can be build directly from GMT-files. Can also be build directly from a list interactively in R.

Usage:
  Rscript rvat.R --buildGeneSet --gmtpath={gmtpath} --output={output} [options]

Arguments:
  gmtpath         Path to a gmt-file
  output          Output file path (geneSetFile format).
  sep             Separator used in input file. Defaults to '\t'.
  ",
  
  buildCorMatrix = "
buildCorMatrix

  Build a block-wise burden correlation matrix, in order to correct for gene-gene correlations in geneSetAssoc. 
  Burden scores should be stored in an aggregateFile object (see aggregate). 
  The size of the blocks are controlled using the maxDist parameter, 
  all gene-gene correlations beyond the block are set to zero. 
  This function is based on previous work (https://github.com/opain/TWAS-GSEA).

Usage:
  Rscript rvat.R --buildCorMatrix --rvbResult={rvbResult} --aggregateFile={aggregateFile} --output={output} [options]

Arguments:
  rvbResult         File path to an rvbResult object (generated using the assocTest method).
  aggregateFile     Path to an aggregateFile.
  memlimit          Maximum number of units to load from the aggregateFile at a time. Defaults to 1000.
  minR2             R2 values < minR2 will be set to zero (leading to increased sparsity).
  not-makePD        Flag to skip forcing correlation matrix to be positive definitive.
  not-absolute      Flag to not make the matrix absolute.
  maxDist           A distance larger than maxDist defines a new block. Defaults to 2.5Mb (5Mb window)
  output            Output file path (R's RDS format).
  ",
  
  geneSetAssoc = "

geneSetAssoc

  Perform self-contained or competitive gene set analyses.

Usage:
  Rscript rvat.R --geneSetAssoc --geneSet={geneSet} --rvbResult={result} --output={output} [options]

Arguments:
  rvbResult       file path to an rvbResult object (generated using the assocTest method)
  geneSet.        file path to a geneSetFile object (generated using the buildGeneSet method)
  scoreMatrix     A matrix (rows = genes, columns = features) can be provided to perform enrichment analyses on continuous values. 
                  These can be used to perform e.g. cell-type enrichment analyses. 
                  Should be a filepath to a matrix stored in R's RDS format.
  cormatrix       File path to a correlation matrix (stored in R's RDS format) with row and column names corresponding to the units in the rvbResult. 
                  Needs to be specified in order to run the 'mlm' (mixed linear model) test. 
                  A burden score correlation matrix can be generated using the --buildCorMatrix method. 
                  The mixed mixed linear model is run using the GENESIS R package.
  condition       Perform conditional analyses. 
                  Input can be a 1) a geneSetFile (set --condition_type=geneSet), in which case the genesets specified in the geneSet parameter will all be conditioned on the gene sets provided here. 
                  2) a comma-delimited list of genesets (set --condition_type=vector) present in specified geneSetFile, 
                  all genesets specified in the geneSetList/geneSetFile will be conditioned on the genesets specified here.
  condition_type  If the --condition argument is specified, use this argument to specifiy the format (geneSet,vector or matrix).
  covar           covariates (comma-delimited)
  test            Tests to perform (comma-delimited). Currently implemented tests are the competitive tests lm,mlm,fisher and self-contained tests ttest,ztest and ACAT.
  threshold       Thresholds for cutoff-based tests (fisher's exact test / glm). Multiple thresholds can be specified comma-delimited.
  Zcutoffs        Cutoffs to apply to the Z-scores (minimum,maximum). Z scores below/above these cutoffs will be set equal to the cutoff.
  INT             Flag to perform inverse normal transformation to Z-scores.
  scoreCutoffs.   If --scoreMatrix is specified, this parameter can be set to cap scores in the scorematrix. It should be a comma-delimited list of 2 values (minimum and maximum sd). 
  minSetSize      Exclude genesets with size < minSetSize.
  maxSetSize      Exclude genesets with size > maxSetSize.
  twoSided        Flag to perform two-sided tests rather than the default one-sided tests.
  memlimit        Maximum number of genesets to process in one go, defaults to 1000.
  ID              ID column in the rvbResult that corresponds with the IDs used in the geneSetList. Defaults to 'unit'.
  output          Output file path
  ",
  vcfInfo2Table="
vcfInfo2Table

  Convert vcf info field to table format. Requires valid vcf where INFO fields are specified in header.

Usage:
  Rscript rvat.R --vcfInfo2Table --vcf={vcf} --output={output}

Arguments:
  vcf                     Vcf file path
  output                  Output path
  not-splitMultiallelic   Flag to not return one row per alternative allele instead of one row per variant.
"
)