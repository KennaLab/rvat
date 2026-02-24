library(rvatData)
gdb <- create_example_gdb()

# retrieve genotypes of a set of variants based on their VAR_ids
varinfo <- getAnno(
  gdb,
  table = "varinfo",
  where = "gene_name = 'SOD1' and ModerateImpact = 1"
)
GT <- getGT(
  gdb,
  VAR_id = varinfo$VAR_id,
  cohort = "pheno"
)

# retrieve genotypes of a set of variants in a varSet
varsetfile <- varSetFile(rvat_example("rvatData_varsetfile.txt.gz"))
varset <- getVarSet(varsetfile, unit = "NEK1", varSetName = "High")
GT <- getGT(
  gdb,
  varSet = varset,
  cohort = "pheno"
)
# see ?varSetFile and ?getVarSet for more details

# retrieve genotypes for a genomic interval
GT <- getGT(
  gdb,
  ranges = data.frame(CHROM = "chr21", start = 31659666, end = 31668931),
  cohort = "pheno"
)

# the `anno` parameter can be specified to include variant annotations in the rowData of the genoMatrix
GT <- getGT(
  gdb,
  ranges = data.frame(CHROM = "chr21", start = 31659666, end = 31668931),
  cohort = "pheno",
  anno = "varInfo",
  annoFields = c(
    "VAR_id",
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "HighImpact",
    "ModerateImpact",
    "Synonymous"
  )
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
varinfo <- getAnno(
  gdb,
  table = "varinfo",
  where = "gene_name = 'UBQLN2' and ModerateImpact = 1"
)
GT <- getGT(
  gdb,
  VAR_id = varinfo$VAR_id,
  cohort = "pheno",
  includeVarInfo = TRUE,
  checkPloidy = "GRCh38"
)

# see ?genoMatrix for more details on the genoMatric class.
