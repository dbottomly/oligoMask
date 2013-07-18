The `oligoMask` package
=========

The `oligoMask` package provides tools for masking out problematic Affymetrix probes 
using VCF files and the oligo package

It has been established that the presence of single nucleotide variants in 
transcript regions that are interograted by microarray probes can cause spurious 
decreases in signal.  In the case of Affymetrix arrays, more than one probe is present 
in a given probeset/transcript-cluster and this provides the opportunity to remove probes 
that could potentially be affected prior to background correction, normalization and 
summarization as in the RMA method. This package is intended to provide the user a convenient 
interface between several packages to faciliate analysis of Affymetrix arrays especially the 
Gene and Exon varieties. Although, originally designed for the analysis of expression data 
from mouse recombinant inbred intercrosses use with human expression arrays is also possible.
