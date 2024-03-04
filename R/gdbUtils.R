#-------------------------------------------------------------------------------
# Creation of gdb files from vcf
#-------------------------------------------------------------------------------

#' Create a gdb file.
#'
#' Creates a new [`gdb`] file and returns a connection object of type gdb-class. 
#' The gdb can be structured and populated using a provided vcf file. 
#' If no input variant file is provided then only an empty gdb is created.
#'
#' @param output Path for output [`gdb`] file
#' @param vcf Input vcf file used to structure and populate gdb. Warning this function makes the following of assumptions: 1) strict adherence to vcf format (GT subfield first element in genotype firelds), 2) multiallelic records have been split, 3) desired genotype QC has already been applied (DP,GQ filters), 4) GT values conform to the set {0/0,0/1,1/0,1/1,./.,0|0,0|1,1|0,1|1,.|.}. Multiallelic parsing and genotype QC can be performed using vcftools and/or accompanying parser scripts included on the rvat github.
#' @param skipIndexes Flag to skip generation of indexes for var and dosage table (VAR_id;CHROM, POS,REF,ALT). Typically only required if you plan to use gdbConcat to concatenate a series of separately generated gdb files before use
#' @param skipVarRanges Flag to skip generation of ranged var table. Typically only required if you plan to use gdbConcat to concatenate a series of separately generated gdb files before use
#' @param overWrite overwrite if `output` already exists? Defaults to `FALSE`, in which case an error is raised.
#' @export
buildGdb=function(output,vcf=c(),skipIndexes=FALSE, skipVarRanges=FALSE,overWrite=FALSE,memlimit=1000)
{

  # Create gdb file
  if (file.exists(output))
  {
    if (overWrite){file.remove(output)} else {stop(sprintf("output gdb already exists '%s'. Must set overWrite=TRUE to replace existing files.",output))}
  }
  mygdb=gdb(output)

  # Create tables
  message(sprintf("%s\tCreating gdb tables",Sys.time()))
  DBI::dbExecute(mygdb,"create table var (VAR_id integer primary key, CHROM text, POS int, ID text, REF text, ALT text, QUAL text, FILTER text, INFO text, FORMAT text);")
  DBI::dbExecute(mygdb,"create table SM (IID text, sex int)")
  DBI::dbExecute(mygdb,"create table dosage (VAR_id integer primary key, GT BLOB);")
  DBI::dbExecute(mygdb,"create table anno (name text,value text,date text)")
  DBI::dbExecute(mygdb,"create table cohort (name text,value text,date text)")
  DBI::dbExecute(mygdb,"create table meta (name text,value text)")
  
  # Import variant records
  if (length(c(vcf))>1){stop("Can only build gdb based on a single input file.")}
  if (length(vcf)>0){populateGdbFromVcf(mygdb,vcf,memlimit=memlimit)}

  # Generate Indexes
  if (! skipIndexes)
  {
    message(sprintf("%s\tCreating var table indexes",Sys.time()))
    DBI::dbExecute(mygdb,"create index var_idx on var (VAR_id)")
    DBI::dbExecute(mygdb,"create index var_idx2 on var (CHROM,POS,REF,ALT)")
    message(sprintf("%s\tCreating SM table index",Sys.time()))
    DBI::dbExecute(mygdb,"create index SM_idx on SM (IID)")
    message(sprintf("%s\tCreating dosage table index",Sys.time()))
    DBI::dbExecute(mygdb,"create index dosage_idx on dosage (VAR_id)")
  }
  
  # Generate var_ranges table
  if(!skipVarRanges) {
    message(sprintf("%s\tCreating ranged var table",Sys.time()))
    addRangedVarinfo(mygdb, overwrite=TRUE)
  }
  
  # Add rvat version to meta table
  DBI::dbExecute(mygdb,"insert into meta values (:name, :value)",
                 params = list(name = "rvatVersion", 
                               value = as.character(packageVersion("rvat"))))
  
  # Complete
  message(sprintf("%s\tComplete",Sys.time()))

  # Return gdb-class object
  return(mygdb)
}


setMethod("populateGdbFromVcf", signature="gdb",
          definition=function(object,vcf,memlimit=1000)
          {
            # Open vcf connection
            if (vcf=="-"){vcf="stdin"} else if (!file.exists(vcf)){stop(sprintf("Input vcf %s does not exist",vcf))}
            if (substr(vcf,nchar(vcf)-2,nchar(vcf))==".gz"){con=gzcon(file(vcf,open='r'))} else {con=file(vcf,open="r")}

            # Skip over vcf meta-data
            while (length(i <- readLines(con,n=1)) > 0)
            {
              if (substr(i,1,2)=="##"){next}
              break
            }

            # Parse vcf header line
            if (!grepl("^#CHROM",i)){stop("Invalid vcf header")}
            header=unlist(strsplit(i,split="\t"))
            width=length(header)
            m=width-9
            if(m<=0){stop("Invalid parsing of vcf header line. No samples detected")}
            message(sprintf("%s sample IDs detected",m))
            DBI::dbWriteTable(object,name="SM",value=data.frame("IID"=header[10:(m+9)],"sex"=0),overwrite=TRUE)

            # Parse vcf records
            message(sprintf("%s\tParsing vcf records",Sys.time()))
            counter=0
            while (length(records <- readLines(con,n=memlimit)) > 0)
            {

              # Increment row counter
              counter=counter+length(records)

              # Parse records
              records=matrix(unlist(stringi::stri_split_fixed(records,"\t")),ncol=width,byrow=TRUE)
              for (i in 1:nrow(records))
              {
                if(grepl(",",records[i,5]))
                {
                  alleles=unlist(strsplit(records[i,5],split=","))
                  gt=substr(records[i,-(1:9),drop=FALSE],1,3)
                  carrierIndex=which(!(gt %in% c(".:0","./.","0/0",".|.","0|0")))
                  carrierGt=strsplit(gt[carrierIndex],split="\\/|\\|")
                  carrierN=length(carrierIndex)
                  for (ai in 1:length(alleles))
                  {
                    # Write variant data
                    insertVarRecord(object,record=c(records[i,1:4],alleles[ai],records[i,6:9]))

                    # Write genotype data
                    gtia=gt
                    for (carrier in 1:carrierN)
                    {
                      gtia[carrierIndex[carrier]]=paste(ifelse(carrierGt[[carrier]]==ai,"1","0"),collapse="/")
                    }
                    insertDosageRecord(object,record=gtia)
                  }
                } else {
                  insertVarRecord(object,record=records[i,1:9])
                  insertDosageRecord(object,record=substr(records[i,-(1:9),drop=FALSE],1,3))
                }

              }
              # Commit
              message(sprintf("%s\tProcessing completed for %s records. Committing to db.",Sys.time(),counter))
            }
          })


setMethod("insertVarRecord", signature="gdb",
          definition=function(object, record)
          {
            DBI::dbExecute(object,"insert into var(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT) values (:chrom, :pos, :id, :ref, :alt, :qual, :flt, :info, :form)",
                           params=list(chrom=record[1],pos=record[2],id=record[3],ref=record[4],alt=record[5],qual=record[6],flt=record[7],info=record[8],form=record[9]))
          })

setMethod("insertDosageRecord", signature="gdb",
          definition=function(object,record)
          {
            record[record=="0/0"]="0"
            record[record=="./0"]="0"
            record[record=="0/1"]="1"
            record[record=="1/0"]="1"
            record[record=="1/1"]="2"
            record[record=="./."]="N"
            record[record=="0|0"]="0"
            record[record=="0|1"]="1"
            record[record=="1|0"]="1"
            record[record=="./1"]="1"
            record[record==".|1"]="1"
            record[record=="1|1"]="2"
            record[record==".|."]="N"
            record[record==".:0"]="N"
            obs=sort(unique(c(record)))
            if (sum(! obs %in% c("0","1","2","N"))>0)
            {
              stop(sprintf("Invalid genotype codes after parsing. Expected (0,1,2,N), observed (%s) ",
                           paste(obs,collapse=",")))
            }
            
            record=I(list(memCompress(paste(record,collapse=""),type="gzip")))
            DBI::dbExecute(object,"insert into dosage(GT) values (:record)",params=list(record=record))
          })

## add GenomicRanges var
setMethod("addRangedVarinfo", 
          signature="gdb",
          definition=
            function(object, overwrite = FALSE){
              
              ## check if 'var_ranges' is already present in the gdb
              if("var_ranges" %in% DBI::dbListTables(object)) {
                if (!overwrite) stop("'var_ranges' is already present in the gdb. Set `overwrite = TRUE` to overwrite.")
                if (overwrite) {
                  message("'var_ranges' is already present in the gdb, it will be overwritten.")
                  dropTable(object, "var_ranges")
                }
              }
              
              chroms <- unique(getAnno(object, "var", fields = "CHROM")$CHROM)
              DBI::dbExecute(object,"create table var_ranges (CHROM text, ranges blob)")
              
              ## generate GRanges for each chromosome anad insert into var_ranges
              for(chrom in chroms) {
                ## get varinfo for current chromosome
                vars <- getAnno(object,
                                "var",
                                where = sprintf("CHROM = '%s'", chrom))
                
                ## convert into genomicRanges 
                gr <- GenomicRanges::GRanges(
                  seqnames = vars$CHROM,
                  ranges = IRanges::IRanges(
                    start = vars$POS,
                    end = vars$POS+stringr::str_length(vars$REF)-1
                  ),
                  VAR_id = vars$VAR_id
                )
                names(gr) <- vars$VAR_id
                
                ## insert GRanges into gdb as blobs per chromosome
                DBI::dbExecute(object, 
                               statement = "insert into var_ranges(CHROM,ranges) values (:CHROM, :ranges)",                                params = list(
                                 CHROM = chrom,
                                 ranges = list(serialize(gr,connection=NULL))))
                
              }
            })

#-------------------------------------------------------------------------------
# Concatenate existing gdb


#' Concatenate gdb databases
#'
#' Function to concatenate [`gdb`] databases. Only retains content of base tables (SM, var, dosage).
#'
#' @param targets File listing full paths of gdbs to concatenate
#' @param output Output gdb file path.
#' @param skipRemap Flag to skip reseting of VAR_id to row id after concatenation. Defaults to `FALSE`.
#' @param skipIndexes Flag to skip generation of standard var and dosage table indexes (VAR_id;CHROM, POS,REF,ALT).
#' Defaults to `FALSE`.
#' @export
concatGdb=function(targets,output,skipRemap=FALSE,skipIndexes=FALSE)
{
  gdb=scan(targets,what="character")
  
  ## check if more than two files are included
  if (length(gdb) < 2) {stop("Require at least 2 valid gdb files for merging")}
  
  ## check if no duplicate files are included
  if (sum(duplicated(gdb)) > 0) {stop("Duplicate files are included in the targets file")}
  
  ## check if filepaths are valid 
  for (i in gdb){if (!file.exists(i)){stop(sprintf("Invalid file path '%s'",i))}}
  
  ## connect to output and check if writable
  if (file.exists(output)){file.remove(output)}
  tryCatch({db=DBI::dbConnect(DBI::dbDriver("SQLite"),output)}, error=function(e){stop(sprintf("Could not write to output path '%s'",output))})

  # Create tables
  message(sprintf("%s\tCreating db tables",Sys.time()))
  DBI::dbExecute(db,"drop table if exists var;")
  DBI::dbExecute(db,"create table var (VAR_id integer, CHROM text, POS int, ID text, REF text, ALT text, QUAL text, FILTER text, INFO text, FORMAT text);")
  DBI::dbExecute(db,"drop table if exists dosage;")
  DBI::dbExecute(db,"create table dosage (VAR_id integer, GT BLOB);")

  # Copy SM table and create new anno and cohort metadata tables
  message("Creating SM, anno and cohort tables")
  DBI::dbExecute(db,"attach :src as src", params=list(src=gdb[1]))
  DBI::dbExecute(db,"create table SM as select * from src.SM")
  DBI::dbExecute(db,"detach src")
  DBI::dbExecute(db,"create table anno (name text,value text,date text)")
  DBI::dbExecute(db,"create table cohort (name text,value text,date text)")

  # Copy var and dosage data from source db
  versions <- vector(mode = "character", length = length(gdb))
  names(versions) <- gdb
  for (i in gdb)
  {
    message(sprintf("Merging '%s'",i))
    DBI::dbExecute(db,"attach :src as src", params=list(src=i))
    DBI::dbExecute(db,"insert into var select * from src.var")
    DBI::dbExecute(db,"insert into dosage select * from src.dosage")
    versions[i] <- dbGetQuery(db, "select * from src.meta where name = 'rvatVersion'")$value
    DBI::dbExecute(db,"detach src")
  }
  version <- unique(versions)
  if ( length(version) > 1 ) {
    warning("Not all input gdbs were created with the same RVAT version, 
we recommend generating the input gdbs with the same RVAT version.")
    version <- version[1]
  }
  

  # Reset VAR_id to row id
  if (!(skipRemap))
  {
    message(sprintf("%s\tReseting VAR_id to rowid",Sys.time()))
    DBI::dbExecute(db,"update var set VAR_id=rowid")
    DBI::dbExecute(db,"update dosage set VAR_id=rowid")
  }

  # Create indexes
  if (!(skipIndexes))
  {
    message(sprintf("%s\tCreating var table indexes",Sys.time()))
    DBI::dbExecute(db,"create index var_idx on var (VAR_id)")
    DBI::dbExecute(db,"create index var_idx2 on var (CHROM,POS,REF,ALT)")
    message(sprintf("%s\tCreating dosage table index",Sys.time()))
    DBI::dbExecute(db,"create index dosage_idx on dosage (VAR_id)")
  }
  
  message(sprintf("%s\tCreating ranged var table",Sys.time()))
  addRangedVarinfo(gdb(output), overwrite=TRUE)
  
  DBI::dbExecute(db,"create table meta (name text,value text)")
  DBI::dbExecute(db,"insert into meta values (:name, :value)",
                 params = list(name = "rvatVersion", 
                               value = version))
  
  message(sprintf("%s\tComplete",Sys.time()))
}

