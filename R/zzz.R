.onLoad <- function(libname, pkgname){
  packageStartupMessage(sprintf("Welcome to RVAT v.%s! \nFor tutorials and updates see: https://kennalab.github.io/rvat/ \nPlease note that in version >= 0.3.0 effect sizes (including SE and CI) of firth and glm tests are on the log(OR) scale.", packageVersion("rvat")))
}
