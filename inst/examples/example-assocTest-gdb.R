library(rvatData)
gdb <- create_example_gdb()
varsetfile <- varSetFile(rvat_example("rvatData_varsetfile.txt.gz"))

# for example purposes, upload a small cohort
cohort <- getCohort(gdb, "pheno")
cohort <- cohort[cohort$IID %in% colnames(GTsmall),]
uploadCohort(gdb, name = "phenosmall", value = cohort)

# run a firth burden test on a binary phenotype
varsetlist <- getVarSet(varsetfile, unit = c("SOD1", "NEK1"), varSetName = "High")
rvb <- assocTest(gdb,
                 varSet = varsetlist,
                 cohort = "phenosmall",
                 pheno = "pheno",
                 covar = c("PC1", "PC2", "PC3", "PC4"),
                 test = "firth",
                 name = "example")


# run ACAT-v and SKAT tests
rvb <- assocTest(gdb,
                 varSet = varsetlist,
                 cohort = "phenosmall",
                 pheno = "pheno",
                 covar = c("PC1", "PC2", "PC3", "PC4"),
                 test = c("acatvfirth", "skat_burden_robust", "skato_robust"),
                 name = "example")

# run a burden test on a continuous phenotype
rvb <- assocTest(gdb,
                 varSet = varsetlist,
                 cohort = "phenosmall",
                 pheno = "age",
                 continuous = TRUE,
                 covar = c("PC1", "PC2", "PC3", "PC4"),
                 test = c("lm", "skat", "acatv"),
                 name = "example")

# run single variant tests on a binary phenotype 
sv <- assocTest(gdb,
                varSet = varsetlist,
                cohort = "phenosmall",
                pheno = "pheno",
                singlevar = TRUE,
                covar = c("PC1", "PC2", "PC3", "PC4"),
                test = c("firth", "glm", "scoreSPA"),
                name = "example",
                minCarriers = 1)

# similarly a list of VAR_ids can be specified instead of a varSetList/varSetFile
sv <- assocTest(gdb,
                VAR_id = 1:20,
                cohort = "phenosmall",
                pheno = "pheno",
                singlevar = TRUE,
                covar = c("PC1", "PC2", "PC3", "PC4"),
                test = c("firth", "glm", "scoreSPA"),
                name = "example",
                minCarriers = 1)

# apply variant filters
rvb <- assocTest(gdb,
                 varSet = varsetlist,
                 cohort = "phenosmall",
                 pheno = "pheno",
                 covar = c("PC1", "PC2", "PC3", "PC4"),
                 test = c("firth", "skat_robust", "acatv"),
                 name = "example", 
                 maxMAF = 0.05,
                 minCarriers = 1,
                 minCallrateVar = 0.9,
                 minCallrateSM = 0.95)

# Perform MAF-weighted burden tests (madsen-browning)
rvb <- assocTest(gdb,
                 varSet = varsetlist,
                 cohort = "phenosmall",
                 pheno = "pheno",
                 covar = c("PC1", "PC2", "PC3", "PC4"),
                 test = c("firth", "skat_robust", "acatv"),
                 MAFweights = "mb")

# CADD-weighted burden test
varsetlist_cadd <- getVarSet(varsetfile, unit = c("SOD1", "NEK1"), varSetName = "CADD")
rvb <- assocTest(gdb,
                 varSet = varsetlist_cadd,
                 cohort = "phenosmall",
                 pheno = "pheno",
                 covar = c("PC1", "PC2", "PC3", "PC4"),
                 test = c("firth", "skat_robust", "acatv"))

# Perform recessive burden test
rvb <- assocTest(gdb,
                 varSet = varsetlist,
                 cohort = "phenosmall",
                 pheno = "pheno",
                 covar = c("PC1", "PC2", "PC3", "PC4"),
                 test = c("firth", "skat_robust", "acatv"),
                 geneticModel = "recessive")

# Resampled burden test 
rvb <- assocTest(gdb,
                 varSet = varsetlist[1],
                 cohort = "phenosmall",
                 pheno = "pheno",
                 covar = c("PC1", "PC2", "PC3", "PC4"),
                 test = c("skat", "skat_burden", "acatv"),
                 name = "example",
                 methodResampling = "permutation",
                 nResampling = 100)

