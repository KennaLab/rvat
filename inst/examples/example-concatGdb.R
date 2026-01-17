library(rvatData)
gdb <- create_example_gdb()

# to illustrate how concatGdb we'll first generate two small gdbs to concatenate
gdb1 <- tempfile()
subsetGdb(
  gdb,
  VAR_id = 1:100,
  output = gdb1
)

gdb2 <- tempfile()
subsetGdb(
  gdb,
  VAR_id = 101:200,
  output = gdb2
)

# write filepaths of gdbs to concatenate to a file
targets <- tempfile()
readr::write_lines(c(gdb1, gdb2), file = targets)
concatgdb <- tempfile()

# concatenate
concatGdb(
  targets = targets,
  output = concatgdb,
)
