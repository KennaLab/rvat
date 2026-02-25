library(dplyr)
library(GenomicRanges)
library(rvatData)

gdbpath <- rvat_example("rvatData.gdb")
testgdb <- withr::local_tempfile()

create_example_gdb <- function (gdbpath = NULL, filepath = NULL) 
{
    if (is.null(filepath)) {
        filepath <- tempfile()
    } else {
        if (file.exists(filepath)) {
            stop(sprintf("%s already exists. Please delete or use a different filepath.", 
                filepath))
        }
    }
    gdb <- test_path("data/rvatData.gdb")
    if (!file.copy(gdb, filepath)) {
        stop("Failed to copy example gdb")
    }
    gdb(filepath)
}
