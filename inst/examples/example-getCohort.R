library(rvatData)
gdb <- create_example_gdb()

# retrieve a cohort
cohort <- getCohort(gdb, cohort = "pheno")
head(cohort)

# retrieve a cohort, keep specified fields
cohort <- getCohort(gdb, cohort = "pheno", fields = c("IID", "sex", "pheno"))
head(cohort)
