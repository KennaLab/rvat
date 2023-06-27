## gdb -------------------------------------------------------------------------
### protected tables
gdb_protected_tables <- c("var","var_ranges", "SM","dosage","anno","cohort","meta")
## assocTest -------------------------------------------------------------------

### all implemented tests
assocTest_tests <- c("firth", "glm", "lm", "scoreSPA", "nbinom", 
                     "skat_burden","skat","skato", 
                     "skat_burden_robust","skat_robust", 
                     "skato_robust", "skat_burden_fwe","skat_fwe",
                     "skato_fwe", "acatv","acatvSPA","acatvfirth")

### singlevariant tests
assocTest_sv_tests <- c("lm", "firth", "glm", "nbinom", "scoreSPA")

#### continuous sv tests
assocTest_sv_cont_tests <- c("lm")

#### binary sv tests
assocTest_sv_bin_tests <- c("firth", "glm", "nbinom", "scoreSPA")

### rvb

#### cont. rvb tests
assocTest_rvb_cont_tests <- c("lm", "skat_burden", "skat", "skato",
                              "skat_burden_fwe", "skat_fwe", "skato_fwe", "acatv")

#### binary rvb tests
assocTest_rvb_bin_tests <- c("firth", "glm", "scoreSPA", "nbinom", "skat_robust", "skato_robust", "skat_burden_robust",
                             "skat_burden_fwe", "skat_fwe", "skato_fwe", "acatv", "acatvSPA","acatvfirth",
                             "skat_burden", "skat", "skato")

#### skat tests
assocTest_skat_tests <- c("skat_burden","skat","skato", "skat_burden_robust", "skat_robust", "skato_robust")

#### tests for which resampling is implemented
assocTest_resampling_tests <- c("skat", "skat_robust", "skat_burden", "skat_burden_robust", 
                                "skato_robust", "acatv")

#### tests for which an offset is implemented
assocTest_offset_tests <- c("firth", "glm", "lm", "scoreSPA", "nbinom", "skat_burden", "skat", "skato",
                            "skat_burden_robust", "skat_robust", "skato_robust", "skat_burden_fwe", "skat_fwe",
                            "skato_fwe", "acatvSPA","acatvfirth"
                            )

#### tests for which an aggregate should be calculated
assocTest_aggregate_tests <- c("firth", "glm", "lm", "nbinom")

## geneSetAssoc ----------------------------------------------------------------
geneSetAssoc_tests <- c("lm", "mlm", "fisher", "ttest", "ztest", "ACAT")
geneSetAssoc_tests_competitive <- c("lm", "mlm", "fisher")
geneSetAssoc_tests_competitive_condition <- c("lm")
geneSetAssoc_tests_competitive_threshold <- c("fisher")
geneSetAssoc_tests_competitive_nothreshold <- geneSetAssoc_tests_competitive[!geneSetAssoc_tests_competitive %in% geneSetAssoc_tests_competitive_threshold]
geneSetAssoc_tests_selfcontained <- c("ttest","ztest", "ACAT")
geneSetAssoc_tests_score <- c("lm", "mlm")

## rvatResult ------------------------------------------------------------------
columns_rvbResults <- list(cohort = c("character", "Rle"),
                           varSetName = c("character", "Rle"),
                           name = c("character", "Rle"),
                           unit = c("character", "Rle"),
                           pheno = c("character", "Rle"),
                           covar = c("character", "Rle"),
                           geneticModel = c("character", "Rle"),
                           MAFweight = c("character", "Rle"),
                           test = c("character", "Rle"),
                           nvar = c("integer", "numeric"),
                           caseCarriers = c("integer", "numeric"),
                           ctrlCarriers = c("integer","numeric"),
                           meanCaseScore = "numeric",
                           meanCtrlScore = "numeric",
                           caseN = c("integer", "numeric"),
                           ctrlN = c("integer", "numeric"),
                           caseCallRate =  c("numeric", "integer"),
                           ctrlCallRate = c("numeric", "integer"),
                           P = "numeric",
                           effect = "numeric",
                           effectSE = "numeric",
                           effectCIlower = "numeric",
                           effectCIupper = "numeric")
columns_rvbResults_rle <- c("cohort", "varSetName", "name", "pheno", "covar", "geneticModel", "MAFweight", "test")
columns_rvbResults_numeric <- c("nvar", "caseCarriers", "ctrlCarriers", "meanCaseScore", "meanCtrlScore", "caseN", "ctrlN", "caseCallRate",
                                "ctrlCallRate", "P", "effect", "effectSE", "effectCIlower", "effectCIupper")

columns_singlevarResults <- list(cohort = c("character", "Rle"),
                                varSetName = c("character", "Rle"),
                                name = c("character", "Rle"),
                                VAR_id = c("character", "Rle"),
                                pheno = c("character", "Rle"),
                                covar = c("character", "Rle"),
                                geneticModel = c("character", "Rle"),
                                test = c("character", "Rle"),
                                caseMAC = c("integer", "numeric"),
                                ctrlMAC = c("integer", "numeric"),
                                caseMAF = c("numeric", "integer"), 
                                ctrlMAF = c("numeric", "integer"), 
                                caseN = c("integer", "numeric"),
                                ctrlN = c("integer", "numeric"),
                                caseCallRate =  c("numeric", "integer"),
                                ctrlCallRate = c("numeric", "integer"),
                                P = "numeric",
                                effect = "numeric",
                                effectSE = "numeric",
                                effectCIlower = "numeric",
                                effectCIupper = "numeric")
columns_singlevarResults_rle <- c("cohort", "varSetName", "name", "pheno", "covar", "geneticModel", "test")
columns_singlevarResults_numeric <- c("caseMAC", "ctrlMAC", "caseMAF", "ctrlMAF", "caseN", "ctrlN", "caseCallRate",
                                      "ctrlCallRate", "P", "effect", "effectSE", "effectCIlower", "effectCIupper")

columns_gsaResults <- list(geneSetName = c("character", "Rle"),
                           test = c("character", "Rle"),
                           covar = c("character", "Rle"),
                           threshold = c("numeric","double"),
                           geneSetSize = c("integer","numeric"),
                           genesObs = c("integer","numeric"),
                           P = "numeric",
                           effect = "numeric",
                           effectSE = "numeric",
                           effectCIlower = "numeric",
                           effectCIupper = "numeric")

columns_gsaResults_rle <- c("geneSetName", "test", "covar")
columns_gsaResults_numeric <- c("threshold", "geneSetSize", "genesObs", "P", "effect", "effectSE", "effectCIlower", "effectCIupper")

columns_rvatResult <- list(
  singlevarResult = columns_singlevarResults,
  rvbResult = columns_rvbResults,
  gsaResult = columns_gsaResults
)

filtercolumns_rvatResult <- list(
  rvbResult = c("unit", "cohort", "varSetName","name", "pheno", "covar", "geneticModel", "MAFweight", "test"),
  singlevarResult = c("VAR_id", "cohort", "varSetName","name", "pheno", "covar", "geneticModel", "test"),
  gsaResult = c("geneSetName", "test","covar", "threshold")
)
