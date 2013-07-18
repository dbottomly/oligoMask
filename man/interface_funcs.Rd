\name{Interface functions/methods}
\alias{maskRMA}
\alias{maskRMA,GeneFeatureSet-method}
\title{
Functions/methods interfacing with other packages
}
\description{
Provides the main functionality behing oligoMask.  These functions/methods allow the user to apply variant masks to oligo classes such as \code{GemeFeatureSet}.
}

\usage{
	\S4method{maskRMA}{GeneFeatureSet}(object, background=TRUE, normalize=TRUE, subset=NULL, target="core", apply.mask=FALSE, mask.params=NULL)
	}

\arguments{
  \item{object}{
	A valid object of class \code{GeneFeatureSet} (for now)}
  \item{background}{Logical-perform RMA background correction?}
  \item{normalize}{Logical-perform quantile normalization?}
  \item{subset}{To be implemented}
  \item{target}{Type of summarization to be performed (e.g. 'core' or 'probeset')}
  \item{apply.mask}{Boolean indicating whether the variant mask should be applied}
  \item{mask.params}{An object of class \code{VariantMaskParams}}
}

\section{Methods}{
	\describe{
		\item{maskRMA}{\code{signature(object = "GeneFeatureSet")} Return a \code{ExpressionSet} object with the probes summarized at the level indicated by \code{target} using the RMA method.  Probes determined to be unreliable are removed prior to the RMA calculations.}	}
}



\author{
Daniel Bottomly
}
\seealso{
\code{\link{rma}}, \code{\linkS4class{GeneFeatureSet}}
}
\examples{
	data(SunGeneFS)
	temp.vcf.db <- new("VcfDB", db.path=om.db.file(), tbsl=SangerTableSchemaList())
	var.mask.pars <- new("VariantMaskParams", var.db=temp.vcf.db)
	
	masked.exprs <- maskRMA(SunGeneFS, apply.mask=TRUE, mask.params=var.mask.pars)
}
\keyword{utilities}

