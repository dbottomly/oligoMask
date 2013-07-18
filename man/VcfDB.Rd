\name{VcfDB-class}
\Rdversion{1.1}
\docType{class}
\alias{VcfDB-class}
\alias{VcfDB}
\alias{vcfDb,VcfDB-method}
\alias{vcfDb}
\alias{tbsl,VcfDB-method}
\alias{tbsl}

\title{Class \code{"VcfDB"}}
\description{
A simple class representing an SQLite database representation of a VCF file for use with the oligoMask package
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("VcfDB", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{db.path}:}{Object of class \code{"character"} containing the file path to the database}
    \item{\code{tbsl}:}{Object of class \code{"TableSchemaList"} containing the schema of the database}
     \item{\code{start.var.table}:}{Object of class \code{"character"} when determining how to join the tables in order to perform the masking procedure, the table to start from.   }
    \item{\code{end.var.table}:}{Object of class \code{"character"} when determining how to join the tables in order to perform the masking procedure, the table to end at.}
	  \item{\code{var.mask.probe.id}:}{Object of class \code{"character"} column that the variant mask database uses to denote the probe ID}
    \item{\code{var.mask.var.id}:}{Object of class \code{"character"} column that the variant mask database uses to denote each variant}
  }
}
\section{Methods}{
  \describe{
    \item{vcfDb}{\code{signature(obj = "VcfDB")}: Returns the filename}
    \item{tbsl}{\code{signature(obj = "VcfDB")}: Returns the \code{TableSchemaList} object}
	 }
}

\author{
Daniel Bottomly
}

\seealso{
\code{\linkS4class{TableSchemaList}}
}
\examples{

temp.vcf.db <- new("VcfDB", db.path=om.db.file(), tbsl=SangerTableSchemaList())

vcfDb(temp.vcf.db)
tbsl(temp.vcf.db)

}
\keyword{classes}
