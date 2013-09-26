test.alignments <- function()
{
    db.name="@DB_NAME@"
    limit.chr <- @LIMIT_CHR@
    bs.genome <- @GENOME_PACKAGE@
    
    db.con <- dbConnect(SQLite(), db.name)
    
    probe.info <- dbGetQuery(db.con, "SELECT * FROM probe_info")
    
    #write to a fastq file
    temp <- strsplit(probe.info$fasta_name, ";")
    fasta.vec <- sapply(temp, "[[", 3)
    names(fasta.vec) <- probe.info$probe_ind
    probe.set <- DNAStringSet(fasta.vec, use.names=TRUE)
    
    use.probe.fa <- tempfile()
    
    writeXStringSet(probe.set, filepath=use.probe.fa)
    
    use.file <- tempfile()
    
    system(paste("@BOWTIE_PATH@ -f -v @NUM_MISMATCH@ -M 1 --best -S -p 10 @GENOME_PATH@", use.probe.fa, paste0(use.file, ".sam")))
    
    asBam(file=paste0(use.file, ".sam"), destination=use.file, overwrite=TRUE)
    sortBam(file=paste0(use.file, ".bam"), destination=use.file, overwrite=TRUE)
    indexBam(paste0(use.file, ".bam"), overwrite=TRUE)
    
    use.param <- ScanBamParam(flag = scanBamFlag(), simpleCigar = FALSE,
         reverseComplement = FALSE, tag = "XM",
         what = scanBamWhat())
    
    scanned.bam <- scanBam(paste0(use.file, ".bam"), param=use.param)
    
    unique.mapped <- which(is.na(scanned.bam[[1]]$tag$XM))
    multi.mapped <- which(is.na(scanned.bam[[1]]$tag$XM) == FALSE & scanned.bam[[1]]$tag$XM == 2)
    non.mapped <- which(is.na(scanned.bam[[1]]$tag$XM) == FALSE & scanned.bam[[1]]$tag$XM == 0)
    
    #check to make sure those unique alignments found in the bowtie alignments only are due to unique mapping to alternative haplotypes/contigs
    bwt.diff.uniqs <- setdiff(scanned.bam[[1]]$qname[unique.mapped], probe.info$probe_ind[probe.info$align_status %in% c("UniqueMapped", "UniqueMappedNotGeno")])
    bwt.diff.chrs <- as.character(scanned.bam[[1]]$rname)[scanned.bam[[1]]$qname %in% bwt.diff.uniqs]
    checkTrue(any(bwt.diff.chrs %in% limit.chr) == FALSE)
    
    #check to make sure those unique alignments found in the snp mask only are due to alignments to multimappings with non-included chromosomes (really just a sanity check by examining whether they are multimapping)
    om.diff.uniqs <- setdiff(probe.info$probe_ind[probe.info$align_status %in% c("UniqueMapped", "UniqueMappedNotGeno")], scanned.bam[[1]]$qname[unique.mapped])
    om.diff.tags <- scanned.bam[[1]]$tag$XM[scanned.bam[[1]]$qname %in% om.diff.uniqs]
    checkTrue(all(om.diff.tags == 2))
    
    #now for those that are unique and in common, check the sanity of the coordinates
    
    probe.aln <- dbGetQuery(db.con, "SELECT probe_ind, probe_chr, probe_start, probe_end FROM probe_info NATURAL JOIN probe_align")
    
    bwt.chr.pos <- with(scanned.bam[[1]], data.frame(probe_ind=qname, probe_chr=as.character(rname), probe_start=pos, probe_end=pos+24, stringsAsFactors=FALSE))
    bwt.chr.pos <- bwt.chr.pos[unique.mapped,]
    bwt.chr.pos <- bwt.chr.pos[bwt.chr.pos$probe_ind %in% probe.aln$probe_ind,]
    probe.aln <- probe.aln[probe.aln$probe_ind %in% bwt.chr.pos$probe_ind,]
    
    bwt.chr.pos$probe_ind <- as.integer(bwt.chr.pos$probe_ind)
    
    #remove those with chromosomes not used in the original probe alignments for additional analysis
    
    bwt.oob.chrs <- bwt.chr.pos[bwt.chr.pos$probe_chr %in% limit.chr == FALSE,]
    probe.aln.oob.chrs <- probe.aln[probe.aln$probe_ind %in% bwt.oob.chrs$probe_ind,]
    
    bwt.chr.pos <- bwt.chr.pos[bwt.chr.pos$probe_chr %in% limit.chr,]
    probe.aln <- probe.aln[probe.aln$probe_ind %in% bwt.oob.chrs$probe_ind == FALSE,]
    
    #this should take care of most
    bwt.chr.pos <- bwt.chr.pos[do.call("order", bwt.chr.pos),]
    probe.aln <- probe.aln[do.call("order", probe.aln),]
    rownames(bwt.chr.pos) <- NULL
    rownames(probe.aln) <- NULL
    
    checkEquals(probe.aln, bwt.chr.pos)
    
    #can actually implement checks if this is the case, but isn't for this example
    checkTrue(nrow(bwt.oob.chrs) == 0)
    
    dbDisconnect(db.con)
}

test.sanger.genotypes <- function()
{
    base.cmd <- "@VCF_QUERY_CMD@ -f '%CHROM:%POS\t%REF[\t%SAMPLE=%TGT]\n' -r $$REGION -s $$SAMPLES $$VCF_FILE"
    
    #in practice, the vcf file and db.name will be placed in the template upon creation..
    
    vcf.files <- list('SNV'="@SNV_PATH@", 'INDEL'="@INDEL_PATH@")
    db.name="@DB_NAME@"
    
    db.con <- dbConnect(SQLite(), db.name)
    
    probe.aligns <- dbGetQuery(db.con, "SELECT * FROM probe_align NATURAL LEFT OUTER JOIN probe_to_snp")
    
    probe.no.genos <- probe.aligns[is.na(probe.aligns$ref_id),]
    
    use.cols <- dbGetQuery(db.con, "SELECT DISTINCT(strain) FROM genotype")[,1] 
    
    pipe.cmd <- sub("$$SAMPLES", paste(use.cols, collapse=","), base.cmd, fixed=TRUE)
    
    #first retrieve those probes not expected to overlap anything
    
    is.homozyg.lines <- function(res)
    {
        temp.res <- data.frame(t(sapply(1:nrow(res), function(x)
                   {
                        paste(use.cols, "=", paste(rep(res[x,2],2), collapse="/"), sep="")
                   })), stringsAsFactors=FALSE)
            
            res <- res[,3:ncol(res)]
            names(res) <- names(temp.res)
            return(sapply(1:nrow(temp.res), function(x) all(t(temp.res)[,x] == t(res)[,x])))
    }
    
    for (i in 1:nrow(probe.no.genos))
    {
        print(i)
        use.interval <- paste0(probe.no.genos[i,2], ":", paste(probe.no.genos[i,3:4], collapse="-"))
        #if no interval exists, it will throw an error about not finding any lines, catch this and test it appropriately
        
        for (j in vcf.files)
        {
            print(j)
            cur.com <- sub("$$VCF_FILE",j, pipe.cmd, fixed=TRUE)
            print (cur.com)
            
            checkTrue(tryCatch(all(is.homozyg.lines(read.delim(pipe(sub("$$REGION", use.interval , cur.com, fixed=TRUE)), colClasses = "character",sep="\t", header=FALSE, stringsAsFactors=FALSE))), error = function(e) any(grepl("no lines available in input", e$message) )))
        }
       
    }
    
    #then do the rest
    
    probe.with.genos <- probe.aligns[is.na(probe.aligns$ref_id) == FALSE,]
    
    #seperate those probes that overlap multiple variants from the rest
    
    dup.probes <- probe.with.genos$probe_ind[duplicated(probe.with.genos$probe_ind)]
    
    dup.probes.with.genos <- probe.with.genos[probe.with.genos$probe_ind %in% dup.probes,]
    
    nd.probes.with.genos <- probe.with.genos[probe.with.genos$probe_ind %in% dup.probes == FALSE,]
    
    all.geno.refs <- dbGetQuery(db.con, "SELECT * FROM reference NATURAL JOIN genotype NATURAL JOIN allele NATURAL JOIN vcf_annot")
    
    split.geno.refs <- split(all.geno.refs, as.character(all.geno.refs$ref_id))
    
    #deal with the non-duplicated probes first
    
    .check.single.probes.geno <- function(probe.to.ref, db.genos)
    {
        use.interval <- paste0(probe.to.ref[1,2], ":", paste(probe.to.ref[1,3:4], collapse="-"))
        
        split.db.genos <- split(db.genos, db.genos$type)
        
        for (z in names(split.db.genos))
        {
            cur.com <- sub("$$VCF_FILE",vcf.files[[z]], pipe.cmd, fixed=TRUE)
        
            genos <- read.delim(pipe(sub("$$REGION", use.interval , cur.com, fixed=TRUE)), colClasses = "character",sep="\t", header=FALSE, stringsAsFactors=FALSE)
            
            keep.rows <- genos[is.homozyg.lines(genos)==FALSE,]
            
            temp.melt.vcf <- melt(keep.rows, id.vars=c("V1", "V2"))
            
            split.samp.genos <- strsplit(temp.melt.vcf$value, "=|\\/")
            
            temp.melt.vcf$strain <- sapply(split.samp.genos, "[", 1)
            temp.geno.ind.list <- unlist(lapply(1:length(split.samp.genos), function(x) rep(x, length(split.samp.genos[[x]])-1)))
            temp.geno.list <- unlist(lapply(1:length(split.samp.genos), function(x) split.samp.genos[[x]][2:length(split.samp.genos[[x]])]))
            
            temp.melt.vcf <- cbind(temp.melt.vcf[temp.geno.ind.list,], genos=temp.geno.list)
            temp.melt.vcf$genos <- as.character(temp.melt.vcf$genos)
            temp.melt.vcf <- temp.melt.vcf[,c("V1", "strain", "genos")]
            names(temp.melt.vcf) <- c("ref", "strain", "alleles")
            temp.melt.vcf <- temp.melt.vcf[!duplicated(temp.melt.vcf),]
            temp.melt.vcf <- temp.melt.vcf[do.call("order", temp.melt.vcf),]
            rownames(temp.melt.vcf) <- NULL
            
            temp.db.genos <- split.db.genos[[z]][,c("seqnames", "start", "strain", "alleles")]
            temp.db.genos$ref <- paste(temp.db.genos$seqnames, temp.db.genos$start, sep=":")
            temp.db.genos <- temp.db.genos[,c("ref", "strain", "alleles")]
            temp.db.genos <- temp.db.genos[!duplicated(temp.db.genos),]
            temp.db.genos <- temp.db.genos[do.call("order", temp.db.genos),]
            rownames(temp.db.genos) <- NULL
            
            checkEquals(temp.melt.vcf, temp.db.genos)
        }
    }
    
    for (i in 1:nrow(nd.probes.with.genos))
    {
        .check.single.probes.geno(nd.probes.with.genos[i,], split.geno.refs[[as.character(nd.probes.with.genos$ref_id[i])]])
    }
    
    split.dups <- split(dup.probes.with.genos, dup.probes.with.genos$probe_ind)
    
    for (i in 1:length(split.dups))
    {
        cur.dups <- split.dups[[i]]
        .check.single.probes.geno(probe.to.ref=cur.dups[1,], do.call("rbind", split.geno.refs[as.character(cur.dups$ref_id)]))
    }
}
