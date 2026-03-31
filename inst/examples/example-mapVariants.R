library(rvatData)
library(rtracklayer)
library(GenomicRanges)
gdb <- create_example_gdb()

# map variants to gene models
ranges <- GRanges(
  seqnames = c("chr21", "chr4"),
  ranges = IRanges(
    start = c(31659666, 169369704),
    end = c(31668931, 169612632)
  ),
  gene_name = c("SOD1", "NEK1")
)

mapVariants(gdb, ranges = ranges, uploadName = "gene", verbose = FALSE)

# similarly, ranges can be a data.frame
ranges <- data.frame(
  CHROM = c("chr21", "chr4"),
  start = c(31659666, 169369704),
  end = c(31668931, 169612632),
  gene_name = c("SOD1", "NEK1")
)

mapVariants(
  gdb,
  ranges = ranges,
  uploadName = "gene",
  verbose = FALSE,
  overWrite = TRUE
)

# often you'd want to map variants to a large set of ranges, such as ensembl models
# mapVariants supports several file formats, including gff/gtf, bed and ranges

# map variants using a gtf file
gtffile <- tempfile(fileext = ".gtf")
rtracklayer::export(
  makeGRangesFromDataFrame(ranges),
  con = gtffile,
  format = "gtf"
)

mapVariants(
  gdb,
  gff = gtffile,
  uploadName = "gene",
  verbose = FALSE,
  overWrite = TRUE
)

# map variants using a bed file

bedfile <- tempfile(fileext = ".bed")
rtracklayer::export(
  makeGRangesFromDataFrame(ranges),
  con = bedfile,
  format = "bed"
)
mapVariants(
  gdb,
  bed = bedfile,
  uploadName = "gene",
  verbose = FALSE,
  overWrite = TRUE
)

# see the variant annotation tutorial on the rvat website for more details
