library(rvatData)
varsetfile <- varSetFile(rvat_example("rvatData_varsetfile.txt.gz"))
gdb <- create_example_gdb()
ranges <- getRanges(varsetfile, gdb = gdb)
