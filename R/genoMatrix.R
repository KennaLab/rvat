#===============================================================================
# genoMatrix objects
#===============================================================================


setMethod("show", signature="genoMatrix",
          function(object)
            {
            cat("rvat genoMatrix Object\n")
            callNextMethod()
            })


# Constructor ------------------------------------------------------------------
#' @rdname genoMatrix
#' @usage NULL
genoMatrix=function(GT, SM, VAR_id, w=1, ploidy="diploid", varSetName="unnamed", unit="unnamed", cohortname="unnamed",verbose=TRUE)
{

  # Check if there are any duplicated VAR_ids
  if (sum(duplicated(VAR_id)) > 0) stop("Duplicated VAR_ids are provided.")
  # Validate provided ploidy values and ensure rowData will have correct nrow
  ploidyLevels=unique(ploidy)
  if (length(w)==1){w=rep(w,length(VAR_id))}

  # Label GT matrix
  colnames(GT)=SM$IID
  rownames(GT)=VAR_id

  # Remove samples with missing IID values
  missing=is.na(SM$IID)
  if (sum(missing)>0)
  {
    if(verbose) message(sprintf("%s/%s samples in the gdb are present in cohort '%s'",sum(!missing),nrow(SM), cohortname ))
    GT=GT[,!missing,drop=FALSE]
    SM=SM[!missing,,drop=FALSE]
    
  }

  # Construct SummarizedExperiment
  GT=SummarizedExperiment::SummarizedExperiment(assays=list(GT=GT),
                                                colData = S4Vectors::DataFrame(SM),
                                                rowData = S4Vectors::DataFrame(ploidy=ploidy, w=unname(w)))

  # set meta data
  S4Vectors::metadata(GT)$ploidyLevels=ploidyLevels
  S4Vectors::metadata(GT)$m=ncol(GT)
  S4Vectors::metadata(GT)$nvar=nrow(GT)
  S4Vectors::metadata(GT)$varSetName=varSetName
  S4Vectors::metadata(GT)$unit=unit
  S4Vectors::metadata(GT)$cohort=cohortname
  S4Vectors::metadata(GT)$geneticModel="allelic"
  S4Vectors::metadata(GT)$imputeMethod="none"

  # create genoMatrix object
  GT=new("genoMatrix",GT)

  # Reset genotype dosages in accordance with male / female / missing gender at sites with XnonPAR or YnonPAR ploidy
  GT=.resetSexChromDosage(GT)

  # Calculate allele frequencies
  SummarizedExperiment::rowData(GT)$AF=getAF(GT)

  # return
  GT
}

#' @importFrom S4Vectors setValidity2
S4Vectors::setValidity2("genoMatrix", function(object)
  {
  msg=NULL

  # Validate ploidyLevels
  if (mean(S4Vectors::metadata(object)$ploidyLevels %in% c("diploid","XnonPAR","YnonPAR"))<1)
  {msg=c(msg,"Ploidy must be set to one of diploid, XnonPAR or YnonPAR")}

  # Validate sample sex field
  if (! ("sex" %in% colnames(SummarizedExperiment::colData(object)))){msg=c(msg,"'sex' is a mandatory field for input SM tables. Accepted values for sex include: 0 (missing), 1 (male) or 2 (female). If needed you can encode all samples as '0' for missing.")}
  if (sum(is.na(SummarizedExperiment::colData(object)$sex))>0){msg=c(msg,"'sex' is a mandatory field for input SM tables. Accepted values for sex include: 0 (missing), 1 (male) or 2 (female). If needed you can encode all samples as '0' for missing.")}
  if (mean(SummarizedExperiment::colData(object)$sex %in% c(0,1,2))<1){msg=c(msg,"'sex' is a mandatory field for input SM tables. Accepted values for sex include: 0 (missing), 1 (male) or 2 (female). If needed you can encode all samples as '0' for missing.")}

  # Return
  if (length(msg)){msg} else{TRUE}
  })


setMethod(".resetSexChromDosage", signature="genoMatrix",
          definition=function(object)
            {
              if (mean(S4Vectors::metadata(object)$ploidyLevels=="diploid")==1){return(object)}
              
              ## XnonPAR
              assays(object)$GT[rowData(object)$ploidy == "XnonPAR", object$sex == 1] <- as.numeric(assays(object)$GT[rowData(object)$ploidy == "XnonPAR", object$sex == 1] >= 1)
              assays(object)$GT[rowData(object)$ploidy == "XnonPAR", object$sex == 0] <- NA_real_
              
              ## YnonPAR
              assays(object)$GT[rowData(object)$ploidy == "YnonPAR", object$sex == 1] <- as.numeric(assays(object)$GT[rowData(object)$ploidy == "YnonPAR", object$sex == 1] >= 1)
              assays(object)$GT[rowData(object)$ploidy == "YnonPAR", object$sex %in% c(0,2)] <- NA_real_
              
              ## return
              object
            })


setMethod("getAF", signature="genoMatrix",
          definition=function(object)
          {
            ## allele frequencies can only be calculated when geneticModel == "allelic"
            if(S4Vectors::metadata(object)$geneticModel == "allelic") {
              # Single ploidy level in genoMatrix
              if ( all(S4Vectors::metadata(object)$ploidyLevels == "diploid"))
              {
                AF <- Matrix::rowMeans(assays(object)$GT, na.rm=TRUE)/2
                return(ifelse(is.finite(AF),AF,0))
              } else {
                # Multiple ploidy levels
                AF <- rep(0, S4Vectors::metadata(object)$nvar)
                names(AF) <- rownames(object)
                
                if ("diploid" %in% S4Vectors::metadata(object)$ploidyLevels) {
                  AF[rowData(object)$ploidy == "diploid"] <- Matrix::rowMeans(assays(object)$GT[rowData(object)$ploidy == "diploid",,drop=FALSE], na.rm=TRUE)/2
                }
                
                if( "XnonPAR" %in% S4Vectors::metadata(object)$ploidyLevels) {
                  if(sum(object$sex == 0) > 0) warning(sprintf("GT contains variants with ploidy='XnonPAR', MAF is calculated in samples with non-missing sex info. Note that sex is missing for %s samples", sum(object$sex == 0)))
                  ac_males = Matrix::rowSums(assays(object)$GT[rowData(object)$ploidy == "XnonPAR", object$sex == 1,drop=FALSE] > 0, na.rm=TRUE)
                  ac_females <- Matrix::rowSums(assays(object)$GT[rowData(object)$ploidy == "XnonPAR", object$sex == 2,drop=FALSE], na.rm=TRUE)
                  non_missing_males = Matrix::rowSums(!is.na(assays(object)$GT[rowData(object)$ploidy == "XnonPAR",object$sex == 1,drop=FALSE]))
                  non_missing_females = Matrix::rowSums(!is.na(assays(object)$GT[rowData(object)$ploidy == "XnonPAR",object$sex == 2,drop=FALSE]))
                  AF[rowData(object)$ploidy == "XnonPAR"] <- (ac_males + ac_females)/(non_missing_males + (2*non_missing_females))
                } 
                
                if( "YnonPAR" %in% S4Vectors::metadata(object)$ploidyLevels) {
                  if(sum(object$sex == 0) > 0) warning(sprintf("GT contains variants with ploidy='YnonPAR', MAF is calculated in samples with non-missing sex info. Note that sex is missing for %s samples", sum(object$sex == 0)))
                  AF[rowData(object)$ploidy == "YnonPAR"] <- Matrix::rowMeans(assays(object)$GT[rowData(object)$ploidy == "YnonPAR", object$sex == 1,drop=FALSE],na.rm=TRUE)
                } 
                return(ifelse(is.finite(AF),AF,0))
              }
            } else {
              warning(sprintf("Allele frequencies are set to NA for geneticModel '%s'.", 
                              S4Vectors::metadata(object)$geneticModel
              ))
              return(rep(NA, nrow(object)))
            }
          })

setMethod("getMAF", signature="genoMatrix",
          definition=function(object)
          {
           object <- flipToMinor(object) 
           getAF(object)
          })

setMethod("getAC", signature="genoMatrix",
          definition=function(object)
          {
            ## calculate allele counts if geneticModel == 'allelic'
            if(S4Vectors::metadata(object)$geneticModel == "allelic") {
              
              if (S4Vectors::metadata(object)$imputeMethod != "none") {
                ac <- rowSums(floor(SummarizedExperiment::assays(object)$GT),na.rm = TRUE)
              } else {
                ac <- rowSums(SummarizedExperiment::assays(object)$GT,na.rm = TRUE)
              }
            } else {
              warning(sprintf("Allele counts can't be calculated for geneticModel '%s'. Returning NA vector.", 
                              S4Vectors::metadata(object)$geneticModel
              ))
              ac <- rep(NA, nrow(object))
            }
            ac
          })

setMethod("getMAC", signature="genoMatrix",
          definition=function(object)
          {
            object <- flipToMinor(object) 
            getAC(object)
          })

setMethod("getNCarriers", signature="genoMatrix",
          definition=function(object)
          {
            carriers <- rowSums(SummarizedExperiment::assays(object)$GT >= 1, na.rm=TRUE)
            carriers
          })

setMethod("getCR", signature="genoMatrix",
          definition=function(object, var = TRUE)
          {
            if (var) {
              Matrix::rowMeans(!is.na(SummarizedExperiment::assays(object)$GT))
            } else {
              Matrix::colMeans(!is.na(SummarizedExperiment::assays(object)$GT))
            }
          })


setMethod("summariseGeno", signature="genoMatrix",
          definition=function(object)
          {
            if (metadata(object)$imputeMethod != "none") stop("Please provide a non-imputed genoMatrix.")
            callRate=Matrix::rowSums(!is.na(SummarizedExperiment::assays(object)$GT))/S4Vectors::metadata(object)$m
            ref=Matrix::rowSums(SummarizedExperiment::assays(object)$GT==0,na.rm=TRUE)
            het=Matrix::rowSums(SummarizedExperiment::assays(object)$GT==1,na.rm=TRUE)
            hom=Matrix::rowSums(SummarizedExperiment::assays(object)$GT==2,na.rm=TRUE)
            hweP <- vector("numeric", nrow(object))
              
            if("XnonPAR" %in% S4Vectors::metadata(object)$ploidyLevels) {
              if(sum(object$sex == 0) > 0) warning(sprintf("GT contains variants with ploidy='XnonPAR', hweP is calculated within females, note that sex is missing for %s samples.", sum(object$sex == 0)))
              object_fem <- object[,SummarizedExperiment::colData(object)$sex==2]
              hweP[SummarizedExperiment::rowData(object)$ploidy == "XnonPAR"] <- 
                  hweTest(ref=Matrix::rowSums(SummarizedExperiment::assays(object_fem)$GT[SummarizedExperiment::rowData(object_fem)$ploidy == "XnonPAR",,drop=FALSE]==0, na.rm=TRUE),
                          het=Matrix::rowSums(SummarizedExperiment::assays(object_fem)$GT[SummarizedExperiment::rowData(object_fem)$ploidy == "XnonPAR",,drop=FALSE]==1, na.rm=TRUE),
                          hom=Matrix::rowSums(SummarizedExperiment::assays(object_fem)$GT[SummarizedExperiment::rowData(object_fem)$ploidy == "XnonPAR",,drop=FALSE]==2, na.rm=TRUE),
                          af=SummarizedExperiment::rowData(object_fem)$AF[SummarizedExperiment::rowData(object_fem)$ploidy == "XnonPAR"])
            }
            hweP[SummarizedExperiment::rowData(object)$ploidy == "diploid"] <- hweTest(ref[SummarizedExperiment::rowData(object)$ploidy == "diploid"],
                                                                                       het[SummarizedExperiment::rowData(object)$ploidy == "diploid"],
                                                                                       hom[SummarizedExperiment::rowData(object)$ploidy == "diploid"],
                                                                                       SummarizedExperiment::rowData(object)$AF[SummarizedExperiment::rowData(object)$ploidy == "diploid"])
            
            hweP[SummarizedExperiment::assays(object)$ploidy == "YnonPAR"]=1 
          
            output <- data.frame(
              VAR_id = rownames(object),
              AF = SummarizedExperiment::rowData(object)$AF,
              callRate = callRate,
              geno0 = ref,
              geno1 = het,
              geno2 = hom,
              hweP = hweP,
              stringsAsFactors=FALSE
            )
            return(output)
          })

hweTest=function(ref,het,hom,af)
{
  n=ref+het+hom
  ref0=round(n*(1-af)^2)
  het0=round(n*2*af*(1-af))
  hom0=round(n*(af)^2)
  X=((ref-ref0)^2)/ref0 + ((het-het0)^2)/het0 + ((hom-hom0)^2)/hom0
  X[!(is.finite(X))]=0
  return(1-pchisq(X,df=1))
}


# Setters ------------------------------------------------------------------------

#' @export
setMethod("[",signature=c("genoMatrix"),
          function(x, i, j, drop=TRUE)
          {
            out=callNextMethod()

            #Update variant counts and ploidyLevels in event of variant filtering
            if (!missing(i))
            {
              S4Vectors::metadata(out)$nvar=nrow(out)
              S4Vectors::metadata(out)$ploidyLevels=unique(rowData(out)$ploidy)
            }

            #Update sample counts and variant AF in event of sample filtering
            if (!missing(j))
            {
              S4Vectors::metadata(out)$m=ncol(out)
              if( S4Vectors::metadata(out)$nvar > 0) {
                SummarizedExperiment::rowData(out)$AF <- getAF(out)
              }
            }

            out
          })

#' @export
setMethod("updateGT", signature="genoMatrix",
          definition=function(object, SM = NULL, anno = NULL)
            {
            
            if (!is.null(SM) ) {
              if (!"IID" %in% colnames(SM)) stop("SM should contain an `IID` column")
              if (!(is.data.frame(SM) || is(SM, "DFrame"))) stop("`SM` should be of class `data.frame` or `DFrame`")
              rownames(SM) <- SM[["IID"]]
              
              # Validate input SM table
              ## check if genoData and SM contain same IIDs
              if (!all(SM[["IID"]] %in% colnames(object)) || !all(colnames(object) %in% SM[["IID"]])){stop("IID values in SM table do not match existing genoMatrix column names. 
                                                                                                            Renaming of genoMatrix columns can be performed by running: colnames(genoMatrixObject)=newNames")}
              
              ## match IIDs
              SM <- SM[match(colnames(object),SM$IID),,drop = FALSE]
              
              # Check for non-diploid ploidy
              if (mean(S4Vectors::metadata(object)$ploidyLevels=="diploid")<1){stop("Genotype dosage values at non-diploid sites were previously reset according to the original sample sex values. 
                                                                                     updateGT cannot be used in this setting, please instead rerun getGT with the corrected SM data.")}
              
              # Reset colData
              SummarizedExperiment::colData(object)=S4Vectors::DataFrame(SM)
            }
            
            if (!is.null(anno)) {
              if (!"VAR_id" %in% colnames(anno)) stop("anno should contain a `VAR_id` column")
              if (any(c("ploidy", "w", "AF") %in% colnames(anno))) {stop("`ploidy`, `w` and `AF` are protected rowData column names")}
              if (length(unique(anno$VAR_id)) < length(anno$VAR_id)) {stop ("`anno` shouldn't contain duplicated VAR_ids!")}
              if (!all(rownames(object) %in% anno$VAR_id)) {warning("Not all variants present in the genoMatrix are present in the anno table. Fields for missing variants will be filled with NAs.")}
              rowdata <- SummarizedExperiment::rowData(object)
              rowdata$VAR_id <- rownames(object)
              rowdata <- merge(rowdata[,c("VAR_id", "ploidy", "w", "AF")], anno, all.x = TRUE)
              rowdata$VAR_id = NULL
              SummarizedExperiment::rowData(object) <- rowdata
            }

            # Validate and return
            validObject(object)
            return(object)
            })

#' @export
setMethod("flipToMinor", signature="genoMatrix",
          definition=function(object)
            {
              if(S4Vectors::metadata(object)$geneticModel != "allelic") {
                warning("flipToMinor only applies when geneticModel == 'allelic', genoMatrix is returned unchanged.")
                return(object)
              }
              flip=SummarizedExperiment::rowData(object)$AF>0.5
              
              if(sum(flip) > 0) {
                if (mean(S4Vectors::metadata(object)$ploidyLevels=="diploid")==1)
                {
                  SummarizedExperiment::assays(object)$GT=abs(SummarizedExperiment::assays(object)$GT - 2*matrix(rep(flip,each=S4Vectors::metadata(object)$m),nrow=S4Vectors::metadata(object)$nvar, byrow=TRUE))
                } else
                {
                  GT <- SummarizedExperiment::assays(object)$GT
                  for (i in which(flip))
                  {

                    if (SummarizedExperiment::rowData(object)$ploidy[i]=="XnonPAR")
                    {
                      # Flips 0s and 1s for male
                      for (i2 in which(SummarizedExperiment::colData(object)$sex==1))
                      {
                        GT[i,i2]=abs(GT[i,i2]-1)
                      }
                      # Flip 0s and 2s for female
                      for (i2 in which(SummarizedExperiment::colData(object)$sex==2))
                      {
                        GT[i,i2]=abs(GT[i,i2]-2)
                      }
                    }
                    if (SummarizedExperiment::rowData(object)$ploidy[i]=="YnonPAR")
                    {
                      # Flip 0s and 1s
                      GT[i,]=abs(GT[i,]-1)
                    }
                    
                    if (SummarizedExperiment::rowData(object)$ploidy[i]=="diploid")
                    {
                      # Flip 0s and 2s
                      GT[i,]=abs(GT[i,]-2)
                    }
                    
                  }
                  object <- BiocGenerics:::replaceSlots(object, assays=Assays(SimpleList(GT=GT)))
                }
                SummarizedExperiment::rowData(object)$AF=getAF(object)
              }
              object
            })

#' @export
setMethod("recode", signature = "genoMatrix",
          definition = function(object, geneticModel, imputeMethod, weights, MAFweights)
          {
            if(nrow(object) == 0) {
              if(!missing(imputeMethod)) S4Vectors::metadata(object)$imputeMethod <- imputeMethod
              return(object)
            }
            if(!missing(geneticModel) && S4Vectors::metadata(object)$geneticModel != geneticModel) {
              if(S4Vectors::metadata(object)$geneticModel != "allelic")
                stop("Current geneticModel should be 'allelic' in order to apply dominant or recessive models.")
              if(S4Vectors::metadata(object)$imputeMethod != "none")
                stop("Provide a non-imputed genoMatrix to perform dominant/recessive recoding")
              object <- flipToMinor(object)
              if (geneticModel == "allelic")
              {SummarizedExperiment::assays(object)$GT=SummarizedExperiment::assays(object)$GT} else if (geneticModel == "dominant")
              {SummarizedExperiment::assays(object)$GT=(SummarizedExperiment::assays(object)$GT>0)*1} else if (geneticModel == "recessive")
              {SummarizedExperiment::assays(object)$GT=(SummarizedExperiment::assays(object)$GT==2)*1} else
              {stop(sprintf("%s does not represent a valid genetic model",geneticModel))}

              S4Vectors::metadata(object)$geneticModel <- geneticModel
              # if(geneticModel != "allelic") {
              #   warning(sprintf("Allele frequencies are set to NA for geneticModel '%s'.", 
              #                   S4Vectors::metadata(object)$geneticModel
              #   ))
              #   rowData(object)$AF <- NA
              # }
            }

            if(!missing(imputeMethod) && !is.null(imputeMethod)) {
              # Validate recode method and identify variants with missing genotype dosages
              if (!(imputeMethod %in% c("meanImpute","missingToRef"))){stop(sprintf("'%s' is not a recognized imputation method",imputeMethod))}
              if(S4Vectors::metadata(object)$imputeMethod != "none")
              {message(sprintf("GT is already imputed (method='%s'), skipping imputation", S4Vectors::metadata(object)$imputeMethod))
              } else {
                missingValues=which(Matrix::rowSums(is.na(SummarizedExperiment::assays(object)$GT)) > 0)

                # mean impute missing genotype dosages
                if (imputeMethod=="meanImpute")
                {
                  if (length(missingValues)>0)
                  {
                    GT <- SummarizedExperiment::assays(object)$GT
                    u <- Matrix::rowMeans(GT[missingValues,,drop=FALSE],na.rm=TRUE)
                    for (i in 1:length(missingValues))
                    {
                      gt=GT[missingValues[i],]
                      gt[is.na(gt)]=u[i]
                      GT[missingValues[i],]=gt
                    }
                    object <- BiocGenerics:::replaceSlots(object, assays=Assays(SimpleList(GT=GT)))
                  }
                }

                # Reset missing genotypes to 0
                if (imputeMethod=="missingToRef")
                {
                  SummarizedExperiment::assays(object)$GT[is.na(SummarizedExperiment::assays(object)$GT)]=0
                  SummarizedExperiment::rowData(object)$AF=getAF(object)
                }

                S4Vectors::metadata(object)$imputeMethod <- imputeMethod
              }
            }

            if(!missing(weights)) {
              if(length(weights) == 1) {
                rowData(object)$w <- as.numeric(rep(weights, nrow(object)))
              } else {
                if(length(weights) != nrow(object)) {stop("Length of `weights` should equal the number of variants in the genoMatrix.")}
                rowData(object)$w <- as.numeric(weights)
              }
            }

            if(!missing(MAFweights)) {
              if(!MAFweights %in% c("mb", "none")) stop("`MAFweights` parameter should be either 'none' or 'mb'.")
              if(MAFweights == "mb") {
                rowData(object)$w <- SummarizedExperiment::rowData(object)$w /
                  (sqrt(SummarizedExperiment::rowData(object)$AF*(1-SummarizedExperiment::rowData(object)$AF)))
              }
            }
            # Update 'aggregate' column (if present)
            if("aggregate" %in% colnames(colData(object))) colData(object)$aggregate <- NA_real_

            # Return genoMatrix object
            return(object)
          })


#' @export
setMethod("aggregate", signature = "genoMatrix",
          definition = function(x, returnGT = TRUE)
          {
            
            callRate <- Matrix::rowSums(!is.na(SummarizedExperiment::assays(x)$GT))/S4Vectors::metadata(x)$m
            if(sum(callRate < 1) > 0) {
              stop("genoMatrix contains missing values, impute missing values using the `recode` method.")
            }
            
            if(any(is.na(rowData(x)$w) | is.infinite(rowData(x)$w))) {
              x <- x[!is.na(rowData(x)$w) & !is.infinite(rowData(x)$w),]
            }
            # Generate aggregate counts
            colData(x)$aggregate <- Matrix::rowSums(
              t(SummarizedExperiment::assays(x)$GT) %*%
                diag(SummarizedExperiment::rowData(x)$w,
                     ncol=S4Vectors::metadata(x)$nvar,
                     nrow=S4Vectors::metadata(x)$nvar)
            )
            
            if(returnGT) return(x) else return(colData(x)$aggregate)
          })



#' @export
setMethod("aggregateGT", signature = "genoMatrix",
          definition = function(object, returnGT = TRUE)
          {
            aggregate(x = object, returnGT = returnGT)
          })

#' @export
setMethod("getCarriers", signature = "genoMatrix",
          definition = function(object, 
                                VAR_id = NULL, 
                                colDataFields = NULL, 
                                rowDataFields = NULL, 
                                groupBy = NULL, 
                                aggregate = FALSE,
                                imputeMethod = "meanImpute"
                                ) {
            
            ## check if specified VAR_ids, colDataFields and rowDatafields are present in the data
            if (is.null(VAR_id))  {
              VAR_id = rownames(object)
            } else {
              VAR_id <- as.character(VAR_id)
              if (!all(VAR_id %in% rownames(object))) {stop("Not all specified VAR_ids are present in the genotype matrix!")}
            }
            object <- object[VAR_id,]
            
            if ( is.null(groupBy) ) {
              
              if (!is.null(colDataFields)) {
                if (!all(colDataFields %in% colnames(colData(object)))) {
                  stop ("Not all specified colDataFields are present in colData(object)!")
                }
              }
              
              if (!is.null(rowDataFields)) {
                if (!all(rowDataFields %in% colnames(rowData(object)))) {
                  stop ("Not all specified rowDataFields are present in rowData(object)!")
                }
              }
              
              carriers <- lapply(
                1:nrow(object),
                function(i)  {
                  x <- assays(object)$GT[i,]
                  carriers <- colnames(object)[!is.na(x) & x >= 1]
                  if(length(carriers) > 0) {
                    x <- x[carriers]
                    data.frame(VAR_id = rownames(object)[i], 
                               IID = carriers,
                               genotype = x,
                               stringsAsFactors=FALSE
                    )
                  } else {
                    NULL
                  }
                })
              carriers <- do.call(rbind, carriers)
              
              if(!is.null(carriers) && nrow(carriers) > 0) {
                rownames(carriers)=NULL
                if (!is.null(colDataFields)) {
                  carriers <- dplyr::left_join(
                    carriers,
                    cbind(IID = colnames(object), as.data.frame(colData(object))[,setdiff(colDataFields, "IID"),drop=FALSE]),
                    by="IID"
                  )
                }
                
                if (!is.null(rowDataFields)) {
                  carriers <- dplyr::left_join(
                    carriers,
                    cbind(VAR_id = rownames(object), as.data.frame(rowData(object))[,setdiff(rowDataFields, "VAR_id"),drop=FALSE]),
                    by="VAR_id"
                  )
                }
              } else {
                carriers <- data.frame()
              }
            } else {
              
              if (aggregate) {
               
                ## freqs per cohort
                object <- flipToMinor(object)
                object <- aggregate(recode(object, imputeMethod = imputeMethod))
                
                carriers <- as.data.frame(colData(object))[,c("aggregate", groupBy),drop=FALSE] %>%
                  dplyr::group_by(dplyr::across(dplyr::all_of(groupBy))) %>%
                  dplyr::summarize(
                    meanBurdenScore  = mean(aggregate),
                    .groups="drop"
                  )
                
              } else {
                groupBy <- unique(c("VAR_id", groupBy))
                carriers <- 
                  cbind(tibble::as_tibble(colData(object)[,setdiff(groupBy, "VAR_id"),drop=FALSE]),
                        tibble::as_tibble(t(assays(object)$GT[VAR_id,,drop=FALSE]))) %>% 
                  tidyr::pivot_longer(dplyr::all_of(VAR_id),
                                      names_to = "VAR_id", 
                                      values_to = "genotype") %>%
                  dplyr::group_by(dplyr::across(dplyr::all_of(groupBy))) %>%
                  dplyr::summarize(
                    carrierN = sum(genotype >= 1, na.rm=TRUE),
                    carrierFreq = mean(genotype >= 1,na.rm=TRUE),
                    carrierFreqSE = sd(genotype >= 1, na.rm=TRUE)/sqrt(dplyr::n()),
                    .groups="drop"
                  )
              }
            }
            
            carriers
          })
