\name{VCF to DB conversions}
\alias{populate.db.tbl.schema.list}
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
populate.db.tbl.schema.list(db.con, db.schema, ins.vals, use.tables = NULL, should.debug = FALSE)
make.vcf.table(db.schema, window.size, vcf.name, db.con, probe.grange, vcf.type="SNV", use.tables=NULL, limit=NULL, should.debug=FALSE, vcf.param=NULL, filter.func=NULL, filter.params=list())
create.sanger.mouse.vcf.db(vcf.files, vcf.labels, probe.tab.file, strain.names, bs.genome, db.schema, db.name="test.db", keep.category="main", window.size=1000, max.mismatch=1, limit.chr=NULL, should.debug=FALSE)
}

\arguments{
  \item{db.con}{
	A connection to an SQLite database.
}
  \item{db.schema}{
	An object of class \code{TableSchemaList}.
}
  \item{ins.vals}{
	List or other type of data to be inserted into the database using the function in the \code{dta.func} element of the TableSchemaList.
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
\item{db.name}{
	A character vector containing the name of the output database name.
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

}
\author{
Daniel Bottomly
}
\seealso{
\code{\link{scanVcf}}, \code{\linkS4class{TableSchemaList}}, \code{\linkS4class{ScanVcfParam}}, \code{\linkS4class{BSgenome}}
}
\examples{
	print("hello")
}
\keyword{utilities}

