# General use command line parser
getArgs=function(help="",required=c())
        {
         args=commandArgs(TRUE)
         args=do.call("rbind",lapply(args,parseArgs))
         labels=args[,1]
         args=as.list(args[,2])
         names(args)=labels
         if ("help" %in% names(args)){stop(help)}
         if (length(required)>0){if (mean(required %in% names(args))!=1){stop(help)}}
         message("RVAT arguments detected:")
         for (i in names(args))
	        {message(sprintf("--%s=%s",i,args[[i]]))}
        return(args)
        }

parseArgs=function(arg){return(c(gsub("^[-]*([^=]*)(.*)","\\1",arg),gsub("^[-]*([^=]*)=(.*)","\\2",arg)))}
