library(rvatData)
vcf <- rvat_example("rvatData.vcf.gz")
output <- tempfile()
vcfInfo2Table(
  vcf = vcf,
  output = output
)
