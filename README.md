The `oligoMask` package
=========

The `oligoMask` package provides tools for masking out problematic Affymetrix probes 
using VCF files and the oligo package

It has been established that the presence of single nucleotide variants in 
transcript regions that are interrogated by microarray probes can cause spurious 
decreases in signal.  In the case of Affymetrix arrays, more than one probe is present 
in a given probeset/transcript-cluster and this provides the opportunity to remove probes 
that could potentially be affected prior to background correction, normalization and 
summarization as in the RMA method. This package is intended to provide the user a convenient 
interface between several packages to facilitate analysis of Affymetrix arrays especially the 
Gene and Exon varieties. Although, originally designed for the analysis of expression data 
from mouse recombinant inbred inter-crosses use with human expression arrays is also possible.

Installation
----------

The simplest approach to installing the package is to download and install the files in the latest release either using R CMD INSTALL package\_name.tar.gz or via install.packages("package\_name.tar.gz", repos=NULL, type="source") from an R session.

The package can also be installed manually by first clicking 'Download ZIP' and unzipping the resulting 
'oligoMask-master.zip' file into a convenient directory.  From within R in the same directory as 
'oligoMask-master' type:

install.packages("oligoMask-master", repos=NULL, type="source")

or for instance using the `devtools` package:

install_github(username="dbottomly", repo="oligoMask", ref="master")

How to contribute
---------

Contributions are encouraged through the standard fork/pull procedures.  Feel free to send me an email with any 
questions.

Contact
---------

Dan Bottomly
bottomly@ohsu.edu
