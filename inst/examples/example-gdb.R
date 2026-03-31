library(rvatData)

# build a gdb directly from a vcf file
vcfpath <- rvat_example("rvatData.vcf.gz") # example vcf from rvatData package
gdbpath <- tempfile() # write gdb to temporary file
buildGdb(vcf = vcfpath, output = gdbpath, genomeBuild = "GRCh38")
gdb <- gdb(gdbpath)

# upload variant info
varinfo <- rvat_example("rvatData.varinfo") # example variant info from rvatData package
uploadAnno(gdb, name = "varInfo", value = varinfo)

# upload cohort info
pheno <- rvat_example("rvatData.pheno") # example cohort info from rvatData package
uploadCohort(gdb, name = "pheno", value = pheno)

# list annotations and cohorts present in gdb
listAnno(gdb)
listCohort(gdb)

# retrieve annotations for a genomic interval
varinfo <- getAnno(
  gdb,
  table = "varinfo",
  ranges = data.frame(CHROM = "chr1", start = 11013847, end = 11016874)
)
head(varinfo)
# see ?getAnno for more details

# retrieve cohort
pheno <- getCohort(gdb, cohort = "pheno")
head(pheno)
# see ?getCohort for more details

# delete table
uploadAnno(gdb, name = "varInfo2", value = varinfo)
dropTable(gdb, "varInfo2")

# retrieve gdb metadata
getGdbPath(gdb)
getGenomeBuild(gdb)

# retrieve genotypes
# first extract variant ids that we want to retrieve the genotypes for
# note that `where` here is an SQL compliant
varinfo <- getAnno(
  gdb,
  table = "varinfo",
  where = "gene_name = 'SOD1' and ModerateImpact = 1"
)

# retrieve genotypes for variants extracted above
GT <- getGT(gdb, VAR_id = varinfo$VAR_id, cohort = "pheno")

# retrieve genotypes for a given genomic interval
GT_fromranges <- getGT(
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
# Learn about the getGT method in ?getGT
# Learn about the genoMatrix format in ?genoMatrix

# Learn more about the downstream methods available for the gdb class in the relevant help pages
# e.g. ?mapVariants, ?subsetGdb, ?assocTest, ?writeVcf
close(gdb)
