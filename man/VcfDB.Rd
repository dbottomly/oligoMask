\name{VcfDB-class}
\Rdversion{1.1}
\docType{class}
\alias{VcfDB-class}
\alias{VcfDB}
\alias{vcfDb,VcfDB-method}
\alias{vcfDb}
\alias{tbsl,VcfDB-method}
\alias{tbsl}
\alias{searchCols,VcfDB-method}
\alias{searchDict,VcfDB-method}
\alias{searchTables,VcfDB-method}
\alias{searchCols}
\alias{searchDict}
\alias{searchTables}
\alias{getProbeVars}
\alias{getProbeVars,VcfDB-method}

\title{Class \code{"VcfDB"}}
\description{
A simple class representing an SQLite database representation of a VCF file for use with the oligoMask package
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("VcfDB", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{db.file}:}{Object of class \code{"character"} containing the file path to the database}
    \item{\code{tbsl}:}{Object of class \code{"TableSchemaList"} containing the schema of the database}
     \item{\code{start.var.table}:}{Object of class \code{"character"} when determining how to join the tables in order to perform the masking procedure, the table to start from.   }
    \item{\code{end.var.table}:}{Object of class \code{"character"} when determining how to join the tables in order to perform the masking procedure, the table to end at.}
	  \item{\code{var.mask.probe.id}:}{Object of class \code{"character"} column that the variant mask database uses to denote the probe ID}
    \item{\code{var.mask.var.id}:}{Object of class \code{"character"} column that the variant mask database uses to denote each variant}
    \item{\code{search.cols}:}{Object of class \code{"list"} A list of lists defining labels and a dictionary for common searches.  Each list should contain
    	a name and have the following elements:
    	\code{table} the table the column should be found in
    	\code{column} the column to be searched against
    	\code{dict} a named character vector in the form key=value where the key is a consistent label that can be database agnostic whereas the label is how
    	that label translates to the current database.}
  }
}
\section{Methods}{
  \describe{
    \item{dbFile}{\code{signature(obj = "VcfDB")}: Returns the filename}
    \item{schema}{\code{signature(obj = "VcfDB")}: Returns the \code{TableSchemaList} object}
    \item{searchCols}{\code{signature(obj = "VcfDB")}, name: Return the column(s) to be searched for a given search label}
    \item{searchDict}{\code{signature(obj = "VcfDB")}, name, value=NULL: Return the database values for a given search label and key}
    \item{searchTables}{\code{signature(obj = "VcfDB")}, name:Return the table name associated with a given search label}
    \item{getProbeVars}{\code{signature(obj = "VcfDB")}, name:Returns a \code{data.frame} containing information on the masked probes.}
    }
}

\author{
Daniel Bottomly
}

\seealso{
\code{\linkS4class{TableSchemaList}}
}
\examples{

temp.vcf.db <- new("VcfDB", db.file=om.db.file(), tbsl=SangerTableSchemaList())

dbFile(temp.vcf.db)
schema(temp.vcf.db)

searchCols(temp.vcf.db)

searchDict(temp.vcf.db, "mapping.status", "unique")

searchTables(temp.vcf.db)

}
\keyword{classes}
