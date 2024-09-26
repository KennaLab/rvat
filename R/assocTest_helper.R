# handle covariates ------------------------------------------------------------
.handle_covar <- function(coldata, covar, pheno) {
  
  # check if covar includes duplicates + check if covariates are available
  covar <- .handle_covar_keep_available(coldata = coldata, covar = covar)
  
  # create dummies if covariates include categorial data (factor/character fields)
  covar_types <- unlist(lapply(coldata[,covar, drop = FALSE], FUN = class))
  if(sum(covar_types %in% c("character", "factor")) > 0) {
    output <- .handle_covar_generate_dummies(coldata = coldata, covar = covar, pheno = pheno, covar_types = covar_types)
    covar <- output[["covar"]]
    coldata <- output[["coldata"]]
  }
  
  # drop covariates with zero covariance
  covar <- .handle_covar_zero_covariance(coldata = coldata, covar = covar)
  
  # return covar vector and coldata
  list(coldata = coldata, covar = covar)
}


.handle_covar_keep_available <- function(coldata, covar) {
  # check if duplicate covariates are specified, if so issue warning and keep unique
  if (length(unique(covar)) < length(covar)) {
    warning("Duplicate covariates are specified, unique covariates are kept.")
  }
  covar <- unique(covar)
  
  # check if all covariates are available, return error when not all are available
  if (mean(covar %in% colnames(coldata)) < 1) {
    stop(sprintf("The following covariate(s) are not available: %s",
                 paste(covar[!covar %in% colnames(coldata)], collapse=",")))
  }
  
  # return covar
  covar
}

.handle_covar_generate_dummies <- function(coldata, covar, pheno, covar_types) {
  modmatrix <- as.data.frame(coldata)
  modmatrix[names(covar_types[covar_types == "character"])] <- lapply(modmatrix[names(covar_types[covar_types == "character"])],
                                                                      factor)
  modmatrix <-  model.matrix.lm(as.formula(sprintf("%s ~ %s", pheno, paste(covar, collapse = "+"))), 
                                data = modmatrix, 
                                na.action = "na.pass")
  modmatrix <- modmatrix[,-1]
  colnames(modmatrix) <- make.names(colnames(modmatrix))
  covar <- colnames(modmatrix)
  coldata <- cbind(
    coldata[,setdiff(colnames(coldata), colnames(modmatrix))], modmatrix)
  list(coldata = coldata, covar = covar)
}

.handle_covar_zero_covariance <- function(coldata, covar, warning = TRUE) {
  # drop covariates with zero variance.
  covarKeep <- c()
  for (i2 in covar) { 
    if(var(coldata[,i2], na.rm = TRUE) > 0) {covarKeep <- c(covarKeep, i2)}
  }
  
  if (length(covar) > length(covarKeep) && warning) {
    warning(sprintf(
      "The following covariate(s) have zero covariance: %s",
      paste(covar[!covar %in% covarKeep], collapse = ",")
    ))
  }
  
  # return covariates to keep
  covarKeep
}
