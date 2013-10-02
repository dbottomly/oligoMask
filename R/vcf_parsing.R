#TODO: if no strains/sample names are specified, parse based on the REF and ALT and maybe another field that would indicate populations

setClass(Class="VcfDB", representation=list(db.path="character", tbsl="TableSchemaList", start.var.table="character", end.var.table="character", var.mask.probe.id="character", var.mask.var.id="character"),
         prototype=prototype(start.var.table="reference", end.var.table="probe_info", var.mask.probe.id="probe_id", var.mask.var.id="ref_id"))

setGeneric("vcfDb", def=function(obj, ...) standardGeneric("vcfDb"))
setMethod("vcfDb", signature("VcfDB"), function(obj)
          {
                return(obj@db.path)
          })

setGeneric("tbsl", def=function(obj, ...) standardGeneric("tbsl"))
setMethod("tbsl", signature("VcfDB"), function(obj)
          {
                return(obj@tbsl)
          })

setMethod("show", signature("VcfDB"), function(object)
        {
			#a hack for now to determine whether the object is a package
			split.path <- strsplit(object@db.path, .Platform$file.sep)[[1]]
			if (split.path[length(split.path)] == "package.db" && split.path[length(split.path)-1] == "extdata")
			{
				pack.name <- split.path[length(split.path)-2]
				pkg.des <- packageDescription(pack.name)
				message(pkg.des$Description)
				message(paste("Built on:", gsub(";", "", pkg.des$Built)))
				message(paste("Package version:", pkg.des$Version))
			}
			else
			{
				#probably a standalone database
				message("An object of class VcfDB pointing to a standalone DB")
			}
			
			db.con <- dbConnect(SQLite(), object@db.path)
			
			if (isTRUE(all.equal(om.CC.mogene.2.1.st@tbsl, SangerTableSchemaList())))
			{
				probe.count <- dbGetQuery(db.con, "SELECT COUNT(DISTINCT(probe_id)) FROM probe_info")[,1]
				message(paste("Containing alignments for", probe.count, "probes"))
				
				var.count <- dbGetQuery(db.con, "SELECT COUNT(DISTINCT(ref_id)) FROM reference")[,1]
				strain.count <- dbGetQuery(db.con, "SELECT COUNT(DISTINCT(strain)) FROM genotype")[,1]
				message(paste("Containing", var.count , "variants from", strain.count, "inbred strains"))
			}
			else
			{
				stop("ERROR: Unknown type of TableSchemaList")
			}
			
			invisible(dbDisconnect(db.con))
        })

#summary method
#probe.cat.counts <- dbGetQuery(db.con, "SELECT align_status AS Alignment_Status, COUNT(DISTINCT(probe_ind)) AS Count from probe_info GROUP BY align_status")
#				format(probe.cat.counts)

filter.sanger.vcf <- function(snp.vcf.list, param.list)
{
	if (class(snp.vcf.list) != "list")
	{
		stop("ERROR: snp.vcf.list needs to be a list object")
	}
	if ("strain.names" %in% names(param.list) == FALSE)
	{
		stop("ERROR: filter function filter.sanger.vcf needs param.list to have a named element 'strain.names'")
	}
	
    strain.names <- param.list$strain.names
    
	if (is.character(strain.names) == FALSE || length(strain.names) == 0)
	{
		stop("ERROR: strain.names needs to be a character vector with at least one element")
	}
	
	#also need to check whether all strain.names are in the GENO$GT matrix
	
	#iterate over the list and only keep those elements which have variants for one of the seven strain
    keep.snps.list <- lapply(snp.vcf.list, function(x)
                           {
                                if (nrow(x$GENO$GT) > 0)
                                {
				    diff.names <- setdiff(strain.names, colnames(x$GENO$GT))
				    if (length(diff.names) > 0)
				    {
					    stop(paste("ERROR:in strain.names", paste(diff.names, collapse=","), "differs from colnames of GENO$GT"))
				    }
                                    return(apply(x$GENO$GT[,strain.names, drop=FALSE], 1, function(x) any(x != "0/0")))
                                }
                                else
                                {
                                    return(FALSE)
                                }
                           })
    
    #first filter for those with at least one variant

    keep.records <- sapply(keep.snps.list, any)
    
    snp.vcf.list <- snp.vcf.list[keep.records]
    keep.snps.list <- keep.snps.list[keep.records]
    
    #next remove entries within a record that do not contain variants in one of the strains
    
    stopifnot(all(names(snp.vcf.list) == names(keep.snps.list)))
    
    sub.snp.vcf.list <- lapply(names(snp.vcf.list), function(x)
                           {
                                keep.elems <- keep.snps.list[[x]]
                                
                                #added due to issue concerning duplicate names for indels
                                temp.snp.vcf <- snp.vcf.list[[x]]
                                names(temp.snp.vcf$rowData) <- NULL
                                
                                return(with(temp.snp.vcf, list(rowData=rowData[keep.elems,], REF=REF[keep.elems,], ALT=ALT[keep.elems],
                                                                    GENO=list(GT=GENO$GT[keep.elems,strain.names, drop=FALSE], FI=GENO$FI[keep.elems,strain.names, drop=FALSE]))))
                           })
	
	names(sub.snp.vcf.list) <- names(snp.vcf.list)
    
    return(sub.snp.vcf.list)
    
}

populate.db.tbl.schema.list <- function(db.con, db.schema, ins.vals=NULL, use.tables=NULL, should.debug=FALSE)
{
	if (class(db.con) != "SQLiteConnection")
	{
		stop("ERROR: db.con needs to be of class SQLiteConnection")
	}
	
	if (class(db.schema) != "TableSchemaList")
	{
		stop("ERROR: db.schema needs to be of class TableSchemaList")
	}
	
	if (missing(ins.vals) || is.null(ins.vals))
	{
		stop("ERROR: ins.vals cannot be missing or NULL")
	}
	
    if (missing(use.tables) || is.null(use.tables) || is.na(use.tables))
    {
        use.tables <- schemaNames(db.schema)
    }
    else if (all(use.tables %in% schemaNames(db.schema)) == FALSE)
    {
        stop("ERROR: Invalid values for use.tables")
    }
    
    #schemaNames should be arranged in the order of population
    for(i in use.tables)
    {
        message(paste("Starting", i))
        #if table doesn't exist, then create it
        if (dbExistsTable(db.con, tableName(db.schema, i, mode="normal")) == FALSE)
        {
            if (should.debug) message("Creating database table")
            if (should.debug) message(createTable(db.schema, i, mode="normal"))
            dbGetQuery(db.con, createTable(db.schema, i, mode="normal"))
        }
        
        #then merge with existing databases as necessary

        if (shouldMerge(db.schema, i))
        {
            if (should.debug) message("Creating temporary table for merging")
            
            if (dbExistsTable(db.con, tableName(db.schema, i, mode="merge")))
            {
                stop("ERROR: Temporary tables should not exist prior to this loop")
            }
            
            if (should.debug) message(createTable(db.schema, i, mode="merge"))
            dbGetQuery(db.con, createTable(db.schema, i, mode="merge"))
            
            if (should.debug) message("Adding to temporary table")
            if (should.debug) message(insertStatement(db.schema, i, mode="merge"))
            #first add the data to temporary database
            dbBeginTransaction(db.con)
            dbGetPreparedQuery(db.con, insertStatement(db.schema, i, mode="merge"), bind.data = bindDataFunction(db.schema, i, ins.vals, mode="merge"))
            dbCommit(db.con)
            
            #merge from temporary into main table
            if (should.debug) message("Merging with existing table(s)")
            if (should.debug) message(mergeStatement(db.schema, i))
            dbGetQuery(db.con, mergeStatement(db.schema, i))
            
            #then also drop intermediate tables
            if (should.debug) message("Removing temporary table")
            if (should.debug) message(paste("DROP TABLE", tableName(db.schema, i, mode="merge")))
            dbGetQuery(db.con, paste("DROP TABLE", tableName(db.schema, i, mode="merge")))
        }else
        {
            if (should.debug) message("Adding to database table")
            if (should.debug) message(insertStatement(db.schema, i))
            #add the data to database
            dbBeginTransaction(db.con)
            dbGetPreparedQuery(db.con, insertStatement(db.schema, i), bind.data = bindDataFunction(db.schema, i, ins.vals))
            dbCommit(db.con)
        }
        
    }
}

#based off of the 'GenomeSearching vignette in BSgenome'
run.analysis.strand <- function(probe.dss, bs.genome, seqnames, strands, max.mismatch, debug)
{
        probe.grange <- GRanges()
		
		if (debug == TRUE) message("Aligning to genome")
        
		if (is.character(seqnames) == TRUE)
		{
			seq.lens <- seqlengths(bs.genome)[seqnames]
			seq.grange <- GRanges(seqnames=names(seq.lens), ranges=IRanges(start=1, end=as.numeric(seq.lens)), strand="*")
		}
		else if (class(seqnames) == "GRanges")
		{
			seq.grange <- seqnames
		}
		else
		{
			stop("ERROR: Unexpected class for seqnames, should be either GRanges or character")
		}
		
        for (strand in strands)
        {
                message("Preprocessing probe sequences")
                if (strand == "-")
                {
                        dict0 <- PDict(x=reverseComplement(probe.dss), max.mismatch=max.mismatch, tb.start=NA, tb.end=NA, tb.width=NA, algorithm="ACtree2", skip.invalid.patterns=FALSE)
                }
                else
                {
                        dict0 <- PDict(x=probe.dss, max.mismatch=max.mismatch, tb.start=NA, tb.end=NA, tb.width=NA, algorithm="ACtree2", skip.invalid.patterns=FALSE)
                }
                
				split.grange <- split(seq.grange, as.character(seqnames(seq.grange)))
				
				for (seqname in names(split.grange))
				{	
						cur.grange <- split.grange[[seqname]]
						
						for (i in 1:length(cur.grange))
						{
							message("Retrieving chromosome sequence")
							subject <- subseq(bs.genome[[seqname]], start=start(cur.grange)[i], end=end(cur.grange)[i])
							message(paste("Aligning to", strand, "of chromosome", seqname, paste(start(cur.grange)[i], end(cur.grange)[i], sep="-"),"..."))
							mindex <- matchPDict(dict0, subject, max.mismatch=max.mismatch)
							
							chr.matches <- unlist(mindex)
							match.names <- names(chr.matches)
							
							chr.matches <- shift(chr.matches, shift=as.integer(start(cur.grange)[i]-1))
							
							if (length(chr.matches) > 0)
							{
								probe.grange <- suppressWarnings(append(probe.grange, GRanges(seqnames=Rle(seqname), ranges=chr.matches, strand=Rle(strand), probe.name=match.names)))
							}
						}
						 
				} 
        }
        
        return(probe.grange)
}

#probe.tab.file <- "/Users/bottomly/Desktop/resources/microarray/MoGene-2_1-st-v1.mm10.probe.tab/MoGene-2_1-st-v1.mm10.probe.tab"
#library(BSgenome.Mmusculus.UCSC.mm10)
#bs.genome <- BSgenome.Mmusculus.NCBI.Build37
#max.mismatch <- 1
#debug <- TRUE
#save(probe.aligns, file="temp_probe_aligns_ncbi37.RData")
#keep.category="main"
#vcf.labels <- c("SNV", "INDEL")
#vcf.files <- c("/Users/bottomly/Desktop/resources/vcfs/mgp.v2.snps.annot.reformat.vcf.gz", "/Users/bottomly/Desktop/resources/vcfs/mgp.v2.indels.annot.reformat.vcf.gz")
#db.schema <- new("TableSchemaList")
#db.name <- "test.db"
#window.size=1000

create.sanger.mouse.vcf.db <- function(vcf.files, vcf.labels, probe.tab.file, strain.names, bs.genome, db.schema, db.name="test.db", keep.category="main", window.size=1000, max.mismatch=1, limit.chr=NULL, should.debug=FALSE, package.info=NULL)
{
        if (is.null(vcf.files) == TRUE || is.null(vcf.labels) == TRUE || is.character(vcf.files) == FALSE || is.character(vcf.labels) == FALSE || length(vcf.files) != length(vcf.labels))
        {
            stop("ERROR: vcf.labels and vcf.files need to be character vectors of equal length")
        }
        
        if (! all(file.exists(vcf.files)))
        {
            stop("ERROR: all vcf.files do not exist")
        }
        
        if (!file.exists(probe.tab.file))
        {
            stop("ERROR: probe.tab.file does not exist")
        }
        
        if (is.character(strain.names) == FALSE || length(strain.names) == 0)
        {
            stop("ERROR: strain.names needs to be a positive length character vector")
        }
        
        if (file.exists(db.name))
        {
            unlink(db.name, recursive=TRUE)
        }
		
	 if (missing(limit.chr) || is.null(limit.chr))
	 {
		 limit.chr <- seqnames(bs.genome)
	 }
	 else if (is.character(limit.chr) == TRUE && all(limit.chr %in% seqnames(bs.genome)) == FALSE)
	 {
		 stop("ERROR: If specified as a character, limit.chr needs to be an element of seqnames(bs.genome)")
	 }
	 else if (class(limit.chr) == "GRanges" && (length(limit.chr) != 1 || all(as.character(seqnames(limit.chr)) %in% seqnames(bs.genome)) == FALSE))
	 {
		 stop("ERROR: If specified as a GRanges object, limit.chr needs to be of length 1 and have its seqnames be in seqnames(bs.genome)")
	 }
	 else if (class(limit.chr) %in% c("character", "GRanges") == FALSE)
	 {
		 stop("ERROR: Unexpected class for limit.chr needs to be either character or GRanges or NULL")
	 }
	 
	 if (class(bs.genome) != "BSgenome")
	 {
		 stop("ERROR: bs.genome needs to be of class BSgenome")
	 }
        
	 if (is.null(package.info) == FALSE && class(package.info) == "list" && all(names(package.info) %in% c("AUTHOR", "AUTHOREMAIL", "BOWTIE_PATH", "GENOME_PATH", "VCF_QUERY_CMD", "VCF_TYPE")))
	 {
		 # browser()
		  actual.db.name <- file.path(db.name, "inst", "extdata", "package.db")
		  actual.tbsl.name <- file.path(db.name, "inst", "extdata", "tbsl.RData")
		  
		  package.desc <- packageDescription("oligoMask")
		  
		  var.type.str <- paste("list(", paste(paste(vcf.labels, "=", paste0("'", vcf.files, "'")) , collapse=","), ")")
		  
		  #for now the chip manufacturer is always Affy, maybe have the ability to change at some point...
		  
		  if (class(limit.chr) == "GRanges")
		  {
			   use.chr.vec <- unique(as.character(seqnames(limit.chr)))
		  }
		  else
		  {
			  use.chr.vec <- limit.chr 
		  }
		  
		  syms <- list(VERSION=package.desc$Version, MANUF="Affymetrix", CHIPNAME=gsub("-", "_", strsplit(basename(probe.tab.file), "\\.")[[1]][1]), LIC=package.desc$License, TBSLDATA=basename(actual.tbsl.name), DB_NAME=basename(actual.db.name),
			  GENOME_PACKAGE=bs.genome@seqs_pkgname, LIMIT_CHR=paste0("c(", paste(paste0("'",use.chr.vec, "'"), collapse=","),")"), VAR_TYPE=var.type.str,
			  NUM_MISMATCH=as.character(max.mismatch))
		  
		  syms <- append(syms, package.info)
	 
		  Biobase::createPackage(pkgname=db.name, destinationDir=".", originDir=system.file(file.path("extdata", "oligoMask.template"), package = "oligoMask"), symbolValues=syms, unlink=TRUE)
		  
		  #the extdata directory doesn't get transferred over probably as it is empty...
		  if (file.exists(dirname(actual.tbsl.name))==FALSE)
		  {
			   dir.create(dirname(actual.tbsl.name), recursive=TRUE)
		  }
		  
		  save(db.schema, file=actual.tbsl.name)
		       
		  db.name <- actual.db.name
	 }
	 else if (is.null(package.info) == FALSE && class(package.info) == "list" && all(names(package.info) %in% c("AUTHOR", "AUTHOREMAIL", "BOWTIE_PATH", "GENOME_PATH", "VCF_QUERY_CMD", "VCF_TYPE")) == FALSE)
	 {
		  stop("ERROR: if a list is supplied to package.info, it needs to have the 'AUTHOR', 'AUTHOREMAIL', 'BOWTIE_PATH', 'GENOME_PATH', 'VCF_QUERY_CMD' and 'VCF_TYPE' names supplied to it")
	 }
	 else if (is.null(package.info) == FALSE && class(package.info) != "list")
	 {
		  stop("ERROR: need to either set package.info as NULL or provide a list")
	 }
	
        db.con <- dbConnect(SQLite(), db.name)
        
        sanger.vcf.param <- ScanVcfParam(trimEmpty=FALSE, fixed="ALT", geno=c("GT", "FI"), info=NA, strand="*")
        
        probe.tab <- read.delim(probe.tab.file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
        
        chr.tab <- table(unlist(lapply(vcf.files, function(x) rownames(header(scanVcfHeader(x))[["contig"]]))))
        common.chrs <- names(chr.tab)[chr.tab == length(vcf.files)]
        
	 if (length(common.chrs) == 0)
	 {
			 message("NOTICE: chromosome annotations not found in VCF, using those from the supplied genome")
			 common.chrs <- as.character(seqnames(bs.genome))
			 
	 }else if(seqnameStyle(Seqinfo(seqnames=common.chrs)) != seqnameStyle(bs.genome))
	 {
                stop("ERROR: Chromosome names are not in the same style")
	 }
        
        if (all(names(probe.tab) %in% c("Probe.ID", "Transcript.Cluster.ID", "probe.x", "probe.y", "assembly", "seqname", "start", "stop", "strand", "probe.sequence", "target.strandedness", "category")))
        {
                keep.tab <- probe.tab[probe.tab$category %in% keep.category,]
                
                dup.probes <- keep.tab[keep.tab$Probe.ID %in% keep.tab$Probe.ID[duplicated(keep.tab$Probe.ID)],]
                
                #there can be duplicate probes assigned to different probesets
                #so first check to see whether probes with the same IDs have the same probe sequence and then keep the unique probes
                if (nrow(dup.probes) > 0)
                {
                        same.seq <- all(tapply(dup.probes$probe.sequence, dup.probes$Probe.ID, function(x) length(unique(x))) == 1)
                
                        if (! same.seq)
                        {
                                stop("ERROR: Expected duplicated probe IDs to have the same sequences")
                        }  
                }
                
                keep.tab <- keep.tab[!duplicated(keep.tab$Probe.ID),]
                
                probe.info <- data.frame(probe_id=keep.tab$Probe.ID, fasta_name=with(keep.tab, paste0(seqname, ":", start,  "-", stop, ";", strand, ";", probe.sequence)), align_status="NonMapped", stringsAsFactors=FALSE)
                
                #as keep.tab and probe.info are aligned
                probe.seqs <- as.character(keep.tab$probe.sequence)
                names(probe.seqs) <- probe.info$fasta_name
                
                stopifnot(sum(duplicated(names(probe.seqs))) == 0)
                
                probe.dss <- DNAStringSet(probe.seqs, use.names=TRUE)
				
                probe.aligns <- run.analysis.strand(probe.dss, bs.genome, limit.chr, c("+", "-"), max.mismatch, should.debug)
                
                #from here create the probe_info and probe_align tables
                
                probe.info <- add.status.probe.info(probe.info, probe.aligns, common.chrs)
				
                if (should.debug) message("Loading probe_info and probe_align tables")
                
                populate.db.tbl.schema.list(db.con, db.schema, list(probe_info=probe.info, probe_align=probe.aligns), use.tables=c("probe_info", "probe_align"), should.debug=should.debug)
                
                if (should.debug) message("Loading data from VCF files")
					
				#keep only the unique probes for mapping
				
				probe.aligns <- probe.aligns[values(probe.aligns)$probe.name %in% probe.info$fasta_name[probe.info$align_status == "UniqueMapped"]]
                
                for(i in 1:length(vcf.files))
                {
                    make.vcf.table(db.schema=db.schema, window.size=window.size, vcf.name=vcf.files[i], db.con=db.con, probe.grange=probe.aligns, vcf.type=vcf.labels[i], use.tables=c("vcf_annot", "reference", "allele", "genotype", "probe_to_snp"), limit=NULL, should.debug=should.debug, vcf.param=sanger.vcf.param,
									   filter.func=filter.sanger.vcf, filter.params=list(strain.names=strain.names))
                }
                
        }
        else
        {
                stop("ERROR: Unsupported tab file specified")
        }
		
	 invisible(dbDisconnect(db.con))
}

#now divide the probes into the different alignment categories
#probes that uniquely mapped to a genotyped chromosome
#probes that uniquely mapped to a non-genotyped contig
#probes that mapped to multiple locations
#probes that did not map

add.status.probe.info <- function(probe.info, probe.aligns, common.chrs)
{
	dup.probes <- unique(values(probe.aligns)$probe.name[duplicated(values(probe.aligns)$probe.name)])
	unique.probes <- setdiff(values(probe.aligns)$probe.name, dup.probes)
	unique.not.geno <- setdiff(values(probe.aligns[as.character(seqnames(probe.aligns)) %in% common.chrs == FALSE])$probe.name, dup.probes)
	
	probe.info$align_status <- "UnMapped"
	probe.info$align_status[probe.info$fasta_name %in% unique.probes] <- "UniqueMapped"
	probe.info$align_status[probe.info$fasta_name %in% unique.not.geno] <- "UniqueMappedNotGeno"
	probe.info$align_status[probe.info$fasta_name %in% dup.probes] <- "MultiMapped"
	
	non.mapped <- setdiff(probe.info$fasta_name, values(probe.aligns)$probe.name)
	#just a sanity check
	if (length(non.mapped) != sum(probe.info$align_status == "UnMapped"))
	{
		stop("ERROR: number of non-mapped is not the expected quantity")
	}
	
	return(probe.info)
}

make.vcf.table <- function(db.schema, window.size, vcf.name, db.con, probe.grange, vcf.type="SNV", use.tables=NULL, limit=NULL, should.debug=FALSE, vcf.param=NULL, filter.func=NULL, filter.params=list())
{
    
    if (missing(limit) || is.null(limit))
    {
        loop.lims <- seq(1, length(probe.grange), window.size)
    }
    else if (is.numeric(limit))
    {	
        loop.lims <- seq(1, length(probe.grange), window.size)
		loop.lims <- loop.lims[1:min(limit, length(loop.lims))]
    }
    else
    {
        stop("ERROR: Invalid limit supplied")
    }
    
    if (missing(vcf.param) || is.null(vcf.param))
    {
        vcf.param <- ScanVcfParam()
    }
	else if (class(vcf.param) != "ScanVcfParam")
	{
		stop("ERROR: vcf.param needs to be an object of class 'ScanVcfParam'")
	}
	
	#make sure the seqlevels are concordant with the current VCF file if seqnames are defined...
	
	temp.head <- rownames(header(scanVcfHeader(vcf.name))[["contig"]])
	
	if (! is.null(temp.head))
	{
		probe.grange <- keepSeqlevels(probe.grange, temp.head)
	}
	
    for (i in loop.lims)
    {
        if (should.debug) message(i)
        if (((i + window.size)-1) > length(probe.grange))
        {
            use.probes <- window(probe.grange, start=i, end=length(probe.grange))
        }else
        {
            use.probes <- window(probe.grange, start=i, width=window.size)
        }
        
        vcfWhich(vcf.param) <- use.probes
        
        snp.vcf.list <- scanVcf(file=vcf.name, param=vcf.param)
        
		if (missing(filter.func) == FALSE && is.null(filter.func) == FALSE)
		{
			sub.snp.vcf.list <- filter.func(snp.vcf.list, filter.params)
		}
		else
		{
			sub.snp.vcf.list <- snp.vcf.list	
		}
        
        if (length(sub.snp.vcf.list) == 0)
        {
           if (should.debug) message("Skipping iteration due to lack of data")
            
        }
        else
        {
            if (should.debug) message("Adding to database")
			
            snp.type.vcf.list <- list(vcf_list=sub.snp.vcf.list, vcf_annot=c(vcf_name=vcf.name, type=vcf.type))
            
            populate.db.tbl.schema.list(db.con, db.schema, snp.type.vcf.list, use.tables, should.debug)
        }
    
    }
}
