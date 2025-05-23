% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allClasses.R, R/varSet.R
\docType{class}
\name{varSetFile}
\alias{varSetFile}
\alias{varSetFile-class}
\alias{show,varSetFile-method}
\alias{metadata,varSetFile-method}
\alias{getGdbId,varSetFile-method}
\alias{getRvatVersion,varSetFile-method}
\alias{length,varSetFile-method}
\alias{listUnits,varSetFile-method}
\alias{listVarSets,varSetFile-method}
\alias{getVarSet,varSetFile-method}
\title{Class to manage interactions with a varSetFile}
\description{
An S4 class to manage interactions with varSetFiles. varSets can be generated from annotations using
the \code{\link[=buildVarSet]{buildVarSet()}} method. Specific units and/or annotations can be loaded using
the \code{getVarSet} method. A varSetFile can be used as input in \code{\link[=assocTest]{assocTest()}} to perform
burden/single variant association tests on the varSets included in the varSetFile.
}
\section{Build a varSetFile}{

\itemize{
\item \code{buildVarSet(x, ...)}: Generate a varSetList or \code{\link[=varSetFile]{varSetFile()}} that stores weighted variant sets
for use in association testing. This can be based on 1) annotations uploaded to the gdb or 2) a data.frame including annotations.
See \code{\link[=buildVarSet]{buildVarSet()}} for details.
}
}

\section{Connect to a varSetFile}{

\itemize{
\item \code{varSetFile(path, memlimit = 5000)}: Connect to a varSetFile object.
}
}

\section{Getters}{

In the following code snippets, x is a varSetFile object.
\itemize{
\item \code{getVarSet(x, unit, varSetName)}: Retrieve varSets for specified units and/or varSetNames.
Output will be of class \code{\link{varSetList}}.
\item \code{listUnits(x)}: Return a vector of all units included in the varSetFile
\item \code{listVarSets(x)}: Return a vector of all varSetNames included in the varSetFile
}
}

\section{Association testing}{

A varSetFile can be directly supplied to the \code{\link[=assocTest]{assocTest()}} gdb method, using the
\code{varSet} parameter. Association tests will then be performed for each varSet included
in the varSetFile.
}

\examples{
library(rvatData)

# Build a varSetFile including variants with a moderate predicted impact
gdb <- create_example_gdb()
varsetfile_moderate <- tempfile()
buildVarSet(object = gdb, 
            output = varsetfile_moderate,
            varSetName = "Moderate", 
            unitTable = "varInfo", 
            unitName = "gene_name",
            where = "ModerateImpact = 1")

# connect to the varSetFile
varsetfile <- varSetFile(varsetfile_moderate)

# list included units and varSets 
units <- listUnits(varsetfile)
head(units)
varsets <- listVarSets(varsetfile)
head(varsets)

# basic operations
length(varsetfile)

# metadata
metadata(varsetfile)
getRvatVersion(varsetfile)
getGdbId(varsetfile)

# retrieve varSets 
varsets <- getVarSet(varsetfile, unit = c("SOD1", "FUS"))
# this returns a varSetList, which is an in-memory representation of varSets
# most of the methods that work on a varSetFile also work on a varSetList (see ?varSetLit for details)
getVarSet(varsets, unit = "SOD1")

# see e.g., ?assocTest and ?aggregate for downstream methods that can loop through varsetfiles and varsetlists.


}
\seealso{
\code{\link{varSetList}}

\code{\link{buildVarSet}}

\code{\link{assocTest}}
}
\keyword{varSetFile}
