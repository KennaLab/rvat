 library(rvatData)

 output <- tempfile()
 gdb <- create_example_gdb()
 writeVcf(gdb, VAR_id = 1:100, output = output)
