library(rvatData)
gdb <- create_example_gdb()

# from data.frame
pheno <- read.table(rvat_example("rvatData.pheno"), header = TRUE)
uploadCohort(object = gdb, name = "cohortinfo", value = pheno)

# similarly, a cohort table can be imported directly from file
filepath <- rvat_example("rvatData.pheno")
uploadCohort(object = gdb, name = "cohortinfo2", value = filepath)
