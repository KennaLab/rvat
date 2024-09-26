#===============================================================================
# Convenience functions to help extract variant annotations from vcf files
#===============================================================================


#' Convert vcf info field to table format
#'
#' Convert vcf info field to table format. Requires valid vcf where INFO fields are specified in header.
#' @param vcf vcf file path.
#' @param output output path
#' @param splitMultiallelic Returns one row per alternative allele instead of one row per variant. Default=TRUE.
#' @export
vcfInfo2Table=function(vcf,output,splitMultiallelic=TRUE)
{
  # Validate input and open connections
  if (vcf=="-"){vcf="stdin"} else if (!file.exists(vcf)){stop(sprintf("Input vcf %s does not exist",vcf))}
  if (substr(vcf,nchar(vcf)-2,nchar(vcf))==".gz"){con=gzcon(file(vcf,open='r'))} else {con=file(vcf,open="r")}
  if (output=="-"){output=stdout()} else {if (file.exists(output)){file.remove(output)}; output=file(output,open="w")}

  # Load INFO fields from vcf meta-data
  fields=list()
  while (length(i <- readLines(con,n=1)) > 0)
  {
    if (substr(i,1,1)!="#"){stop("Invalid vcf header.")}
    if (substr(i,1,6)=="#CHROM"){break}
    if (grepl("^##INFO",i)){fields[[sub("##INFO=<ID=([^,]*),.*","\\1",i)]]="NA"}
  }
  nfields=length(names(fields))
  for (i in names(fields)){message(sprintf("Detected INFO field %s ",i))}
  write(paste(c("CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER",names(fields)),collapse="\t"),file=output,append=TRUE)

  # Parse records (skipping in memory header line)
  sitePattern=sprintf("(^%s).*",paste(rep("[^\t]*",8),collapse="\t"))
  while (length(i <- readLines(con,n=1)) > 0)
  {
    out=fields
    i=unlist(strsplit(sub(sitePattern,"\\1",i),"\t"))
    if (i[[8]]=="."){write(paste(c(i[1:7],unlist(out)),collapse="\t"),file=output,append=TRUE); next}
    for (field in strsplit(unlist(strsplit(i[[8]],";")),split="="))
    {
      out[[field[1]]]=field[2]
    }
    if (length(out) != nfields){warning(sprintf("Observed field not described in vcf meta information - %s\n",paste(i,collapse="|")))}
    if (splitMultiallelic)
    {
      alleles=unlist(strsplit(i[5],split=","))
      for (ai in alleles)
      {
        write(paste(c(i[1:4],ai,i[6:7],unlist(out)),collapse="\t"), file=output, append=TRUE)
      }
    } else
    {
      write(paste(c(i[1:7],unlist(out)),collapse="\t"), file=output, append=TRUE)
    }
  }
}