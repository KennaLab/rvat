library(dplyr)
library(GenomicRanges)
library(rvatData)

gdbpath <- rvat_example("rvatData.gdb")
testgdb <- withr::local_tempfile()
