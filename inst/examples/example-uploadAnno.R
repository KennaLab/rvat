library(rvatData)
gdb <- create_example_gdb()

# from data.frame
varinfo <- read.table(rvat_example("rvatData.varinfo"), header = TRUE)
uploadAnno(object = gdb, name = "varInfo", value = varinfo, overWrite = TRUE)

# similarly, an annotation table can be imported directly from file
filepath <- rvat_example("rvatData.varinfo")
uploadAnno(object = gdb, name = "varInfo2", value = filepath)

# if the annotation table already includes a 'VAR_id' field
# the `skipRemap` parameter can be set to TRUE to skip mapping based on
# CHROM,POS,REF,ALT.
anno <- mapVariants(gdb, 
                    ranges = data.frame(CHROM = "chr21", start = 31659666, end = 31668931, gene_name = "SOD1"),
                    verbose = FALSE)
uploadAnno(object = gdb, name = "gene", value = anno, skipRemap = TRUE)
