library(rvatData)
vcfpath <- rvat_example("rvatData.vcf.gz")
gdbpath <- tempfile()

# build a gdb from vcf.
# the genomeBuild parameters stores the genome build in the gdb metadata
# this will be used to assign ploidies on sex chromosomes (diploid, XnonPAR, YnonPAR)
buildGdb(
  vcf = vcfpath,
  output = gdbpath,
  genomeBuild = "GRCh38"
)

# for large vcfs, the memlimit parameter can be lowered
buildGdb(
  vcf = vcfpath,
  output = gdbpath,
  genomeBuild = "GRCh38",
  memlimit = 100,
  overWrite = TRUE
)

# see ?gdb for more information on gdb-files,
# see ?concatGdb for concatenate gdb databases
