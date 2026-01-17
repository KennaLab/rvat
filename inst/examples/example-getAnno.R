library(rvatData)
gdb <- create_example_gdb()

# retrieve full anno table
varinfo <- getAnno(gdb, table = "varInfo")
head(varinfo)

# extract a genomic range
varinfo <- getAnno(gdb, 
                   table = "varInfo", 
                   ranges = data.frame(CHROM = "chr1", start = 11013847, end = 11016874))
head(varinfo)

# keep only specified fields
varinfo <- getAnno(gdb, 
                   table = "varInfo", 
                   fields = c("VAR_id", "CHROM", "POS", "REF", "ALT", "ModerateImpact"),
                   ranges = data.frame(CHROM = "chr1", start = 11013847, end = 11016874))
head(varinfo)

# the `where` parameter can be used to to pass an SQL-compliant where clause t
varinfo <- getAnno(gdb, 
                   table = "varInfo", 
                   where = "gene_name = 'SOD1' and ModerateImpact = 1")
head(varinfo)


# the `inner` and `left` parameters can be used to perform inner and left join operations respectively
# e.g. we can use the `inner` parameter to filter e.g. based on a table containing QC-passing variants
# for example:
uploadAnno(gdb, name = "QCpass", value = data.frame(VAR_id = 1:100), skipRemap = TRUE, verbose = FALSE)
varinfo <- getAnno(gdb, 
                   inner = "QCpass",
                   table = "varInfo")
