\name{oligoMask-package}
\alias{oligoMask-package}
\alias{oligoMask}
\alias{om.db.file}
\alias{om.tab.file}
\alias{om.vcf.file}
\alias{om.lo.file}
\alias{test.pd.1.0.st.v1}
\alias{SunGeneFS}

\docType{package}
\title{
Tools for masking out problematic Affymetrix probes using VCF files and the oligo package
}
\description{
It has been established that the presence of single nucleotide variants in transcript regions that are interograted by microarray probes can cause spurious decreases in signal.  In the case of Affymetrix arrays, more than one probe is present in a given probeset/transcript-cluster and this provides the opportunity to remove probes that could potentially be affected prior to background correction, normalization and summarization as in the RMA method. This package is intended to provide the user a convient interface between several packages to faciliate analysis of Affymetrix arrays especially the Gene and Exon varieties. Although, originally designed for the analysis of expression data from mouse recombinant inbred intercrosses use with human expression arrays is also possible.
}
\details{
\tabular{ll}{
Package: \tab oligoMask\cr
Version: \tab 0.99.01\cr
Date: \tab 2013-5-13\cr
Imports: \tab methods, DBI, RSQLite,VariantAnnotation,batch,reshape2,graph,RBGL,oligo,oligoClasses,BSgenome,Biostrings,seqnames.db\cr
Depends: \tab R (>= 3.0.0)\cr
Suggests: \tab RUnit, BiocGenerics, BSgenome.Mmusculs.UCSC.mm9\cr
License: \tab GPL-3\cr
biocViews: \tab DataImport\cr
Collate: \tab helper_functions.R tableSchemaList.R vcf_parsing.R VariantMaskParams.R\cr
Built: \tab R 3.0.0; ; 2013-06-10 21:05:11 UTC; unix\cr
}

Index:
\preformatted{
TableSchemaList-class   Class '"TableSchemaList"'
VariantMaskParams-class
                        Class '"VariantMaskParams"'
VcfDB-class             Class '"VcfDB"'
maskRMA                 Method implementing variant masking along with the RMA algorithm
create.sanger.mouse.vcf.db	A function for the conversion of Sanger VCF files to an SQLite database
om.db.file	Function returning a path to a sample oligoMask SQLite database
om.tab.file	Function returning a path to a sample tab-delimited probe sequence file from Affymetrix
om.vcf.file	Function returning a path to a subset of a Sanger VCF file containing SNVs for several mouse inbred strains.
om.lo.file	Function returning a path to the lifted over coordinates from the om.tab.file
test.pd.1.0.st.v1	A small sample 'AffyGenePDInfo' class for testing purposes
SunGeneFS	A subset of the Affy Gene ST expression data from Sun et al. 2012
}

Further information is available in the following vignettes:
\tabular{ll}{
\code{oligoMask} \tab  (source)\cr
}
}
\author{
Daniel Bottomly

Maintainer: Daniel Bottomly <bottomly@ohsu.edu>
}

\references{
Sun, Wei, et al. "Transcriptome atlases of mouse brain reveals differential expression across brain regions and genetic backgrounds." G3: Genes| Genomes| Genetics 2.2 (2012): 203-211.
}

\keyword{ package }
