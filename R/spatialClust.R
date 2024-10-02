#' @rdname spatialClust
#' @usage NULL
#' @export
setMethod("spatialClust", signature="gdb",
          definition = function(object,output,varSetName,unitTable,unitName,windowSize,overlap,intersection=c(),where=c(),weightName=1, posField="POS",minTry=5,warning=TRUE)
          {
            if (length(windowSize)!=length(overlap)){stop(sprintf("Number of provided window sizes (%s) does not match number of provided overlaps (%s)",
                                                                  paste(windowSize,collapse=","),paste(overlap,collapse=",")))}
            output=gzfile(output,"w")
            query=sprintf("select distinct %s as unit, VAR_id, %s as weight, %s as POS from %s", unitName, weightName, posField, unitTable)
            if (length(intersection)>0){intersection=unlist(strsplit(intersection,split=","))}
            for (i in intersection)
            {
              query=sprintf("%s inner join %s using (VAR_id)",query,i)
            }
            if (length(where)>0){query=sprintf("%s where %s",query, where)}
            query=sprintf("select unit, group_concat(VAR_id) as VAR_id, group_concat(weight) as weight, '%s' as varSetName, group_concat(POS) as POS from (%s) x group by unit",varSetName, query)
            
            metadata <- list(
              rvatVersion = as.character(packageVersion("rvat")),
              gdbId = getGdbId(object),
              genomeBuild = getGenomeBuild(object),
              creationDate = as.character(round(Sys.time(), units = "secs"))
            )
            .write_rvat_header(filetype = "varSetFile", 
                               metadata = metadata, 
                               con = output)
            
            handle=RSQLite::dbSendQuery(object,query)
            while (!RSQLite::dbHasCompleted(handle))
            {
              varSet=RSQLite::dbFetch(handle,n=1)
              runDistanceCluster(varSet=varSet,posField=posField,windowSize=windowSize,overlap=overlap,output=output,minTry=minTry,warning=warning)
            }
            RSQLite::dbClearResult(handle)
            close(output)
          }
)


# runDistanceCluster
# Parser to handle input/ output of variant data for distanceCluster
runDistanceCluster=function(varSet,posField,windowSize,overlap,output,minTry=5,warning=TRUE)
{
  unit=varSet$unit
  varSetName=varSet$varSetName
  varSet=data.frame(VAR_id=unlist(strsplit(varSet$VAR_id,",")),
                    weight=unlist(strsplit(varSet$weight,",")),
                    POS=unlist(strsplit(as.character(varSet$POS),",")),
                    stringsAsFactors=FALSE
  )
  
  if (nrow(varSet) < minTry ) {
    if (warning) warning(sprintf("varSet contains fewer than %s variants ('minTry' parameter), all variants are returned as a single cluster.", minTry))
    for (wi in 1:length(windowSize)) {
      varSeti=paste(paste(unit,'0',sep="_"),paste(varSet$VAR_id,collapse=","),paste(varSet$weight,collapse=","),paste(varSetName,windowSize[wi],overlap[wi],sep="_"),sep="|")
      write(varSeti,output,append=TRUE)
    }
    return(NULL)
  }
  
  if (posField=="HGVSc") {
    varSet$POS=as.integer(stringr::str_replace(varSet$POS,"c.([0-9]*).*","\\1"))
  } else {
    varSet$POS=as.integer(varSet$POS)
  }
  
  for (wi in 1:length(windowSize))
  {
    parts=distanceCluster(pos=varSet$POS,names=varSet$VAR_id,windowSize=windowSize[wi],overlap=overlap[wi],return_what="names", write_output=FALSE, path_output=NA,warning=warning)
    for (part in names(parts))
    {
      varSeti=apply(varSet[varSet$VAR_id %in% parts[[part]],],2,paste,collapse=",")
      varSeti=paste(paste(unit,part,sep="_"),varSeti["VAR_id"],varSeti["weight"],paste(varSetName,windowSize[wi],overlap[wi],sep="_"),sep="|")
      write(varSeti,output,append=TRUE)
    }
  }
}

# distanceCluster
# Adapted from https://github.com/heidefier/cluster_wgs_data (Fier, GenetEpidemiol, 2017). 
# Partitions a supplied set of variant until the mean distance between variants equals the median (as would be expected for a homogeneous poisson process).

distanceCluster<-function(pos,names,windowSize=100,overlap=50,return_what="names", write_output=TRUE, path_output=NA,warning=TRUE){
  
  if (!(return_what %in% c("names","positions"))){
    stop('Please enter "names" OR "positions" as argument for return_what')
  }
  
  if (write_output & is.na(path_output)){
    stop('Please specify a valid output name in path_output','\n')
  }
  
  b<-as.numeric(pos)    
  comb<-data.frame(b = b, names = names,stringsAsFactors=FALSE)
  comb<-comb[order(b),]
  b<-as.numeric(comb[,1])
  names<-comb[,2]
  bdiff<-diff(b)
  
  # calculate lambda under the assumption that distances follow exp distribution
  
  lambda<-log(2)/median(bdiff)
  
  #observed mean values in the windows
  
  obsMeanWindow <- zoo::rollapply(bdiff,windowSize,mean,by=overlap,align="left")
  if (length(obsMeanWindow) == 0) {
    if(warning) warning("No valid partitioning for specified parameters, all variants are returned as a single cluster.")
    return(list("0" = names))
  }
  
  #theoretical mean value given medians in the windows
  ws<-as.numeric(windowSize)
  sh<-as.numeric(overlap)
  
  theoMeanWindow <-zoo::rollapply(bdiff,ws,function(x) 1/(log(2)/median(x)),by=sh,align="left")
  
  #actual positions
  raPositions <- zoo::rollapply(b,ws+1,function(x) x,by=sh,align="left")
  
  maximas<-vector("list",length(obsMeanWindow))
  for (i in 1:length(obsMeanWindow)){
    if (obsMeanWindow[i] <= theoMeanWindow[i]) {
      maximas[[i]]<-NA
    } else {
      ord<-order(diff(raPositions[i,]),decreasing=TRUE)
      j<-1
      excl<-mean(diff(raPositions[i,])[-c(ord[1])])
      while (theoMeanWindow[i]<excl){
        j<-j+1
        excl<-mean(diff(raPositions[i,])[-c(ord[1:j])]) #The largest distances are removed until the median-informed theoretical means are larger than the mean of the edited window
      }
      j<-j-1
      maximas[[i]]<-raPositions[i,][c(ord[1:j])+1] 
    }
  }
  maxima <- unlist(maximas)
  maxima <- maxima[!is.na(maxima)]
  maximaf <- sort(unique(maxima))
  
  if (length(maximaf) == 0) {
    if(warning) warning("No valid partitioning for specified parameters, all variants are returned as a single cluster.")
    return(list("0" = names))
  }
  
  if (!is.na(path_output)) file.create(path_output)
  if (return_what=="positions") {
    lo_final<-split(b,findInterval(b,maximaf))
    ll<-as.numeric(lapply(lo_final,length))
    if (return_what=="positions" && write_output) lapply(lo_final, write, path_output, append=TRUE, ncolumns=max(ll) )
  } else {
    lo_final<-split(names,findInterval(b,maximaf))
    ll<-as.numeric(lapply(lo_final,length))
    if (return_what=="names" && write_output)  lapply(lo_final, write, path_output, append=TRUE, ncolumns=max(ll) )
  }
  
  return(lo_final)
  
}
