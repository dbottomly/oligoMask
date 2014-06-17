\name{VariantMaskParams-class}
\Rdversion{1.1}
\docType{class}
\alias{VariantMaskParams-class}
\alias{VariantMaskParams}
\alias{getProbeDf,VariantMaskParams-method}
\alias{getProbeDf}
\alias{maskDb,VariantMaskParams-method}
\alias{maskDb}
\alias{validProbeQuery,VariantMaskParams-method}
\alias{validProbeQuery}

\title{Class \code{"VariantMaskParams"}}
\description{
A class containing the parameter values used for the filtering procedure.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("VariantMaskParams", ...)} or more conveniently through the 
VariantMaskParams constructor function as VariantMaskParams(var.db, geno.filter=FALSE, rm.unmap=TRUE, rm.mult=TRUE, mask.type="static").
}
\section{Slots}{
  \describe{
    \item{\code{mask.type}:}{Object of class \code{"character"} Type of mask, currently only the 'static' mask is available}
    \item{\code{geno.filter}:}{Object of class \code{"logical"} indicating whether only the probes overlapping the high quality genotypes be used (default is FALSE)}
    \item{\code{rm.unmap}:}{Object of class \code{"logical"} indicating whether the probes that are unmapped should be removed as well (default is TRUE)}
    \item{\code{rm.mult}:}{Object of class \code{"logical"} indicating whether the probes that map to multiple locations should be removed as well (default is TRUE)}
     \item{\code{oligo.probe.id}:}{Object of class \code{"character"} indicating the database column that the oligo package uses to denote the probe ID}
  \item{\code{var.db}:}{Object of class \code{"VcfDB"} corresponding to the VCF database to be used for masking}
  }
}
\section{Methods}{
  \describe{
    \item{getProbeDf}{\code{signature(object = "VariantMaskParams"), gene.fs, target, sortBy="fsetid"}: Retrieve a \code{data.frame} describing the probe-(meta)probeset relationships of the probes not masked out based off the database connection in \code{gene.fs} and the summary type specified in \code{target} currently either 'core' or 'probeset'.}
    \item{maskDb}{\code{signature(obj = "VariantMaskParams")}: Retrieve the filename of the variant database}
    \item{validProbeQuery}{\code{signature(object = "VariantMaskParams"), target, should.add=TRUE}: Returns an SQLite query to be used for selecting the probes to be masked out and excluded from the basic query performed by oligo.}
	 }
}
\author{
Daniel Bottomly
}
\seealso{
	\code{\link{create.sanger.mouse.vcf.db}}, \code{\linkS4class{VcfDB}}
}
\examples{
data(SunGeneFS)
temp.vcf.db <- new("VcfDB", db.path=om.db.file(), tbsl=SangerTableSchemaList())
var.mask.pars <- new("VariantMaskParams", var.db=temp.vcf.db)

validProbeQuery(var.mask.pars, target="core")

probe.df <- getProbeDf(var.mask.pars, gene.fs=SunGeneFS, target="core")
head(probe.df)

maskDb(var.mask.pars)

}
\keyword{classes}
