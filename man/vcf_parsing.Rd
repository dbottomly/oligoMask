\name{VCF to DB conversions}
\alias{make.vcf.table}
\alias{add.vcf.to.db}
\alias{create.sanger.mouse.vcf.db}
\title{
Functions for conversion of VCF files to SQLite databases
}
\description{
Functionality to create SQLite databases from several types of VCF including those that represent the genetic diversity of populations (e.g. 1000 genomes) or inbred strains (Sanger mouse genomes project).
}
\usage{
make.vcf.table(db.obj, window.size, vcf.name, probe.grange, vcf.type="SNV", use.tables=NULL, limit=NULL, should.debug=FALSE, vcf.param=NULL, filter.func=NULL, filter.params=list())
create.sanger.mouse.vcf.db(sanger.db, vcf.files, vcf.labels, probe.tab.file, strain.names, bs.genome, keep.category="main", window.size=1000, max.mismatch=1, limit.chr=NULL, should.debug=FALSE, package.info=NULL)
}

\arguments{
  
  \item{db.obj}{
    A \code{VcfDB} object
  }
  
  \item{sanger.db}{
    A \code{VcfDB} object containing a SangerTableSchemaList() and path to the database or package to be created.
  }
  
  \item{use.tables}{
	The subset of tables the procedure should be limited to or NULL.
}
  \item{should.debug}{
	Logical indicating whether additional messages should be displayed to the user
}
\item{window.size}{
	The number of aligned probes to be processed at a given time
}
\item{vcf.name}{
	Name of the VCF file
}
\item{strain.names}{
	character vector containing the genotype columns to be used
}
\item{probe.grange}{
	GRanges object containing the probe alignments information
}

\item{vcf.type}{
	Label associate with the VCF file, usually 'SNV' or 'INDEL'
}
\item{limit}{
	Number of iterations to be limited to for testing purposes
}
\item{vcf.param}{
	An object of class \code{ScanVcfParam}
}
\item{filter.func}{
	A function to be used to filter a list as returned by \code{scanVcf}  For example see oligoMask:::filter.sanger.vcf.  If missing or NULL, no
	filtering will be performed.
}
\item{filter.params}{
	A list containing named elements with any additional values to be passed to \code{filter.func} if applicable.
}
\item{vcf.files}{
	A character vector containing the path to one or more VCF files.
}
\item{vcf.labels}{
	A character vector containing the desired label for each VCF file (e.g. 'SNV' or 'INDEL').
}
\item{probe.tab.file}{
	A character vector containing the path to the tab delimited probe sequence file distributed by Affymetrix.
}
\item{bs.genome}{
	An object of class \code{BSgenome}
}

\item{keep.category}{
	Categories of probes to keep from the \code{probe.tab.file} e.g. 'main'.  These will be the probes that go through the SNP masking process and will be 
	candidates for removal due to being impacted by variants.
}
\item{max.mismatch}{
	maximum number of mismatches to allow in the realignments.
}
\item{limit.chr}{
	A character vector containing the chromosomes to consider or NULL in which case all available chromosomes from the specified genome will be used.
}
\item{package.info}{
	Either NULL or a named list containing the following elements: "AUTHOR", "AUTHOREMAIL", "BOWTIE_PATH", "GENOME_PATH", "VCF_QUERY_CMD", "VCF_TYPE".  If NULL is supplied, a database will
	be generated at the path specified in dbFile of the Database object.  If the list is provided, a package will be generated containing the database along with other relevant metadata supplemented by the values supplied in
	package.info.  For this list, AUTHOR and AUTHOREMAIL indicates the author of the package and their email.  BOWTIE_PATH, GENOME_PATH and VCF_QUERY_CMD indicates the path to a bowtie executable, the genome
	used for creation of the database and the path to a vcf-query tool from vcftools or htslib, these are not necessary unless the package is to be tested afterwards and can be set to an arbitrary string.  They can always be added in later.  VCF_TYPE is simply a desription of the type of VCF database created for instance NODvB6 or CC.
}

}
\author{
Daniel Bottomly
}
\seealso{
\code{\link{scanVcf}}, \code{\linkS4class{TableSchemaList}}, \code{\linkS4class{Database}}, \code{\linkS4class{ScanVcfParam}}, \code{\linkS4class{BSgenome}}
}
\examples{
    if (require(BSgenome.Mmusculus.UCSC.mm9))
    {
    	vcf.files <- om.vcf.file()
    	vcf.labels <- "SNV"
    	probe.tab.file <- om.tab.file()
    	#seven strains besides the reference strain that are used in the CC
    	strain.names <- c("CASTEiJ", "AJ", "PWKPhJ", "129S1", "NZO", "NODShiLtJ", "WSBEiJ")
    	bs.genome <- BSgenome.Mmusculus.UCSC.mm9

    	#probably shouldn't do this in real analyses unless you have to
    	seqlevels(bs.genome) <- sub("chr", "", seqlevels(bs.genome))

    	db.schema <- SangerTableSchemaList()
    	db.name=tempfile()

    	keep.category <- "main"
    	window.size <- 1000
    	max.mismatch <- 0
    	should.debug <- TRUE

    	#duplications in this file, remove and recreate...
    	tab.aln <- read.delim(om.lo.file(), sep="\t", header=TRUE, stringsAsFactors=FALSE)
    	limit.chr <- GRanges(seqnames="19", ranges=IRanges(start=min(tab.aln$start), end=max(tab.aln$end)/4), strand="*")
	
	sanger.db <- new("VcfDB", db.file=db.name, tbsl=db.schema)
	
    	create.sanger.mouse.vcf.db(sanger.db, vcf.files, vcf.labels, probe.tab.file, strain.names, bs.genome, keep.category, window.size, max.mismatch, limit.chr, should.debug)
	
    	db.con <- dbConnect(SQLite(), db.name)

    	var.ovls <- dbGetQuery(db.con, "SELECT * FROM probe_to_snp NATURAL JOIN reference NATURAL JOIN probe_align NATURAL JOIN probe_info limit 5")

    	var.ovls
      }

}
\keyword{utilities}

