
full.tests <- FALSE

stopifnot(require(RSQLite))
stopifnot(require(BSgenome.Mmusculus.UCSC.mm9))
stopifnot(require(VariantAnnotation))
stopifnot(require(reshape2))

data("SunGeneFS")

#currently not testing the following trivial methods that are exported

#SangerTableSchemaList, tbsl, vcfDb, om.vcf.file, om.tab.file, VariantMaskParams, maskDb
#searchTables, searchCols, searchDict


test.createTable <- function()
{
    tbsl <- new("TableSchemaList")
    
    valid.tables <- names(tbsl@tab.list)
    
    db.con <- dbConnect(SQLite(), "temp.db")
    
    for(i in valid.tables)
    {
        print(i)
        for (j in c("normal", "merge"))
        {
            print(j)
            
            f.keys <- tbsl@tab.list[[i]]$foreign.keys
            
            #if there are no foreign keys available, don't allow create table statements to be generated
            if (j == "merge" && is.null(f.keys))
            {
                checkException(createTable(tbsl, table.name=i, mode="merge"))
            }
            else
            {
                if (is.null(f.keys) || j == "normal")
                {
                    add.cols <- character(0)
                    add.type <- character(0)
                    add.pk <- integer(0)
                }
                else
                {
                    #the basic table should already exist so can retrieve the previous info on coltypes
                    
                    temp.prag <- do.call("rbind", lapply(names(f.keys), function(x)
                           {
                                dbGetQuery(db.con, paste("pragma table_info(",x,")"))
                           }))
                    
                    key.vals <- sapply(f.keys, "[[", "ext.keys")
                    key.prag <- temp.prag[temp.prag$name %in% key.vals,]
                    add.cols <- key.prag$name
                    add.type <- key.prag$type
                    add.pk <- rep(0L, length(add.type))
                }
                
                prag.tab.name <- ifelse(j=="merge", paste0(i, "_temp"), i)
                
                checkTrue(is.null(dbGetQuery(db.con, createTable(tbsl, table.name=i, mode=j))))
                tab.prag <- dbGetQuery(db.con, paste("pragma table_info(",prag.tab.name,")"))
                sub.prag <- tab.prag[,c("name", "type", "pk")]
                
                col.types <- sapply(strsplit(tbsl@tab.list[[i]]$db.schema, "\\s+"), "[", 1)
                col.names <- tbsl@tab.list[[i]]$db.cols
                is.pk <- as.integer(grepl("PRIMARY KEY", tbsl@tab.list[[i]]$db.schema))
                
                col.names <- append(col.names, add.cols)
                col.types <- append(col.types, add.type)
                is.pk <- append(is.pk, add.pk)
                
                query.dta <- data.frame(name=col.names, type=col.types, pk=is.pk, stringsAsFactors=FALSE)
                
                if (j == "merge")
                {
                    query.dta <- query.dta[query.dta$pk == 0 & query.dta$name %in% sapply(f.keys, "[[", "local.keys") == FALSE,]
                }
                
                ord.prag <- sub.prag[do.call("order", sub.prag),]
                ord.query <- query.dta[do.call("order", query.dta),]
                
                rownames(ord.prag) <- NULL
                rownames(ord.query) <- NULL
                
                checkEquals(ord.prag, ord.query)
            }
        }
        
    }
    
    dbDisconnect(db.con)
    file.remove("temp.db")
}

test.insertStatement <- function()
{
    set.seed(123)
    
    tbsl <- new("TableSchemaList")
    
    valid.tables <- names(tbsl@tab.list)
    
    db.con <- dbConnect(SQLite(), "temp.db")
    
    for(i in valid.tables)
    {
        print(i)
        for(j in c("normal", "merge"))
        {
            print(j)
            f.keys <- tbsl@tab.list[[i]]$foreign.keys
            
            if (j == "merge" && is.null(f.keys))
            {
                checkException(insertStatement(tbsl, i, j))
            }
            else
            {
                #first create the tables
                checkTrue(is.null(dbGetQuery(db.con, createTable(tbsl, table.name=i, mode=j))))
                
                prag.tab.name <- ifelse(j=="merge", paste0(i, "_temp"), i)
                tab.prag <- dbGetQuery(db.con, paste("pragma table_info(",prag.tab.name,")"))
                
                #create a couple lines of fake data to insert into the database
                
                ins.dta <- as.data.frame(matrix(sample.int(10000, 10*nrow(tab.prag)), ncol=nrow(tab.prag), nrow=10, dimnames=list(NULL, tab.prag$name)), stringsAsFactors=fALSE)
                
                for(p in colnames(ins.dta))
                {
                    if (tab.prag$type[tab.prag$name == p] == "TEXT")
                    {
                        ins.dta[,p]  <- as.character(ins.dta[,p])
                    }
                }
                
                #load into the database
                
                dbBeginTransaction(db.con)
                checkTrue(is.null(dbGetPreparedQuery(db.con, insertStatement(tbsl, i, mode=j), bind.data = ins.dta)))
                dbCommit(db.con)
                
                #check whether it respects should.ignore
                
                ignore.match <- regexpr(pattern="INSERT\\s+OR\\s+IGNORE", text=insertStatement(tbsl, i, mode=j), perl=TRUE)
                
                if (tbsl@tab.list[[i]]$should.ignore)
                {
                    checkTrue(ignore.match != -1)
                }
                else
                {
                    checkTrue(ignore.match == -1)
                }
            }
        }
    }
    
    dbDisconnect(db.con)
    file.remove("temp.db")
}

test.mergeStatement <- function()
{  
    tbsl <- new("TableSchemaList")
    
    valid.tables <- names(tbsl@tab.list)
    
    for(i in valid.tables)
    {
        #again if there are no foreign keys make sure the query dies
        f.keys <- tbsl@tab.list[[i]]$foreign.keys
        print(i)
        if (is.null(f.keys))
        {
            checkException(mergeStatement(tbsl, i))
        }
        else
        {
            cur.stat <- mergeStatement(tbsl, i)
            
            #is the table definition consistent
            
            tab.match <- regexpr(pattern=paste0(i, "\\s+\\(\\s+([\\w+_,]+)\\s+\\)"), text=cur.stat, perl=TRUE)
            tab.str <- substr(cur.stat, start=attr(tab.match, "capture.start"), stop=attr(tab.match, "capture.start")+attr(tab.match, "capture.length")-1)
            split.tab <- strsplit(tab.str, ",")[[1]]
            
            tab.cols <- tbsl@tab.list[[i]]$db.cols
            tab.cols <- tab.cols[tbsl@tab.list[[i]]$db.schema != "INTEGER PRIMARY KEY AUTOINCREMENT"]
            
            checkTrue(length(intersect(split.tab, tab.cols)) == length(union(split.tab, tab.cols)))
            
            #is the select statement consistent with the table definition
            
            select.match <- regexpr(pattern=paste0("SELECT\\s+", tab.str), text=cur.stat, perl=TRUE)
            
            checkTrue(select.match != -1)
            
            #are the joins sane
            
            join.base <- sapply(f.keys, function(x) paste0("\\(", paste(x$ext.keys, collapse=","), "\\)"))
            
            join.str <- paste0(i, "_temp\\s+", paste(paste("JOIN", names(join.base), "USING", join.base, sep="\\s+"), collapse="\\s+"))
            
            join.match <- regexpr(pattern=join.str, text=cur.stat, perl=TRUE)
            
            checkTrue(join.match != -1)
            
            #is it respecting should.ignore
            
            ignore.match <- regexpr(pattern="INSERT\\s+OR\\s+IGNORE", text=cur.stat, perl=TRUE)
            
            if (tbsl@tab.list[[i]]$should.ignore)
            {
                checkTrue(ignore.match != -1)
            }
            else
            {
                checkTrue(ignore.match == -1)
            }
        }
    }
}

test.populate.db.tbl.schema.list <- function()
{
    #populate a subset of the tbsl using the example data  
    
    .tc.func <- function(x)
    {
        temp.x <- x[!duplicated(x$Transcript.Cluster.ID),"Transcript.Cluster.ID", drop=FALSE]
        names(temp.x) <- "TC_ID"
        return(temp.x)
    }
    
    .probe.func <- function(x)
    {
        temp.x <- x[,c("Probe.ID", "probe.x", "probe.y", "Transcript.Cluster.ID")]
        names(temp.x) <- c("Probe_ID", "Probe_X", "Probe_Y", "TC_ID")
        return(temp.x)
    }
    
    tab.list <- list(transcript_cluster=list(db.cols=c("TC_PK", "TC_ID"),
                                db.schema=c("INTEGER PRIMARY KEY AUTOINCREMENT", "INTEGER"),
                                db.constr="",
                                dta.func=.tc.func, should.ignore=FALSE, foreign.keys=NULL),
        probe=list(db.cols=c("Probe_PK", "Probe_ID", "Probe_X", "Probe_Y", "TC_PK"),
                   db.schema=c("INTEGER PRIMARY KEY AUTOINCREMENT", "INTEGER", "INTEGER", "INTEGER", "INTEGER"),
                   db.constr="",
                   dta.func=.probe.func, should.ignore=FALSE, foreign.keys=list(transcript_cluster=list(ext.keys="TC_ID", local.keys="TC_PK"))))
    
    tbsl <- new("TableSchemaList", tab.list=tab.list)
    
    db.con <- dbConnect(SQLite(), "temp.db")
    
    probe.tab.file <- om.tab.file()
    
    probe.tab <- read.delim(probe.tab.file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    
    populate.db.tbl.schema.list(db.con, db.schema=tbsl, ins.vals=probe.tab, use.tables=NULL, should.debug=TRUE)

    test.query <- dbGetQuery(db.con, "SELECT Probe_ID, Probe_X, Probe_Y, TC_ID FROM probe JOIN transcript_cluster USING (TC_PK)")
    names(test.query) <- c("Probe.ID", "probe.x", "probe.y", "Transcript.Cluster.ID")
    
    test.query <- test.query[do.call("order", test.query),]
    rownames(test.query) <- NULL
    
    sub.probe.tab <- probe.tab[,c("Probe.ID", "probe.x", "probe.y", "Transcript.Cluster.ID")]
    sub.probe.tab <- sub.probe.tab[do.call("order", sub.probe.tab),]
    rownames(sub.probe.tab) <- NULL
    
    checkEquals(test.query, sub.probe.tab)
    
    dbDisconnect(db.con)
    file.remove("temp.db")
    
    #should be able to easily populate a single table without dependencies
    
    db.con <- dbConnect(SQLite(), "temp.db")
    
    populate.db.tbl.schema.list(db.con, db.schema=tbsl, ins.vals=probe.tab, use.tables="transcript_cluster", should.debug=TRUE)
    
    comp.tab <- .tc.func(probe.tab)
    comp.tab <- cbind(TC_PK=1:nrow(comp.tab), comp.tab)
    rownames(comp.tab) <- NULL
    
    test.query.2 <- dbGetQuery(db.con, "SELECT * FROM transcript_cluster")
    
    checkEquals(comp.tab, test.query.2)
    
    dbDisconnect(db.con)
    file.remove("temp.db")
    
    #should throw an error if attempting to populate a table with dependencies
    
    db.con <- dbConnect(SQLite(), "temp.db")
    checkException(populate.db.tbl.schema.list(db.con, db.schema=tbsl, ins.vals=probe.tab, use.tables="probe", should.debug=TRUE))
    
    dbDisconnect(db.con)
    file.remove("temp.db")
}

test.make.vcf.table <- function()
{
    
    db.schema <- new("TableSchemaList")
    #ust a relatively strange number < the length of probe.grange
    window.size <- 109
    vcf.name <- om.vcf.file()
    db.con <- dbConnect(SQLite(), "temp.db")
    
    lo.probe.dta <- read.delim(om.lo.file(), sep="\t", header=TRUE, stringsAsFactors=FALSE)
    #duplications exist here, so remove them ahead of time
    lo.probe.dta <- lo.probe.dta[!duplicated(lo.probe.dta),]
    
    probe.grange <- with(lo.probe.dta, GRanges(seqnames=Rle(sub("chr", "", seqnames)), ranges=IRanges(start=start, end=end), strand=Rle("*")))
    
    vcf.type <- "SNV"
    use.tables <- c("vcf_annot", "reference", "allele", "genotype")
    
    if (full.tests == TRUE)
    {
        limit <- NULL
    }else
    {
        limit <- 5
    }
    
    
    should.debug <- TRUE
    vcf.param <- ScanVcfParam(trimEmpty=FALSE, fixed="ALT", geno=c("GT", "FI"), info=NA, strand="*")
    filter.func <- oligoMask:::filter.sanger.vcf 
    filter.params <- list(strain.names=c("CASTEiJ", "AJ", "PWKPhJ", "129S1", "NZO", "NODShiLtJ", "WSBEiJ"))
    
    make.vcf.table(db.schema, window.size, vcf.name, db.con,probe.grange, vcf.type, use.tables, limit, should.debug, vcf.param, filter.func, filter.params)
    
    all.dta <- dbGetQuery(db.con, "SELECT seqnames, start, end, filter, alleles, allele_num, geno_chr, strain FROM reference NATURAL JOIN allele NATURAL JOIN genotype NATURAL JOIN vcf_annot")
    
    if (is.null(limit) == FALSE)
    {
        vcfWhich(vcf.param) <- probe.grange[1:(window.size*limit)]
    }else
    {
        vcfWhich(vcf.param) <- probe.grange
    }
    
    vcf.list <- scanVcf(vcf.name, param=vcf.param)
    filt.vcf <- filter.func(vcf.list, filter.params)
    
    filt.vcf.dta <- do.call("rbind", lapply(filt.vcf, function(x)
           {
                temp.dta <- cbind(as.data.frame(x$rowData), ref=as.character(x$REF), alt=sapply(x$ALT, paste, collapse=","), gt=x$GENO$GT, fi=x$GENO$FI)
           }))
    
    #the same variant can impact several probes, so remove the duplicates as the conversion to the database should do so
    
    filt.vcf.dta <- filt.vcf.dta[!duplicated(filt.vcf.dta),]
    
    #filter of TRUE implies that all of the strains passed the FI filter for each genotype
    
    filt.vcf.dta$filter <- apply(filt.vcf.dta[,grep("fi\\.", colnames(filt.vcf.dta))], 1, function(x) all(is.na(x) == FALSE & x==1))
    
    filt.vcf.dta <- filt.vcf.dta[,-grep("fi\\.", colnames(filt.vcf.dta))]
    
    melt.vcf.dta <- melt(filt.vcf.dta, measure.vars=colnames(filt.vcf.dta)[grep("gt\\.", colnames(filt.vcf.dta))], variable.name="strain")
    melt.vcf.dta$strain <- sub("gt\\.", "", melt.vcf.dta$strain)
    split.gt <- do.call("rbind", strsplit(melt.vcf.dta$value, "/"))
    melt.vcf.dta <- cbind(melt.vcf.dta, chr=split.gt)
    sec.melt.vcf <- melt(melt.vcf.dta, measure.vars=colnames(melt.vcf.dta)[grep("chr", colnames(melt.vcf.dta))], value.name="allele_num", variable.name="geno_chr")
    sec.melt.vcf$geno_chr <- sub("chr\\.", "", sec.melt.vcf$geno_chr)
    
    unk.ref.alts <- strsplit(paste(".", sec.melt.vcf$ref, sec.melt.vcf$alt, sep=","), ",")
    suppressWarnings(sec.melt.vcf$allele_num <- as.numeric(sec.melt.vcf$allele_num))
    sec.melt.vcf$allele_num <- ifelse(is.na(sec.melt.vcf$allele_num), -1, sec.melt.vcf$allele_num)
    
    sec.melt.vcf$alleles <- sapply(1:length(unk.ref.alts), function(x) unk.ref.alts[[x]][sec.melt.vcf$allele_num[x]+2])
    
    #go from here checking these results against the populated database
    
    #convert everyone to compatable datatypes prior to comparison
    
    for (i in 1:ncol(sec.melt.vcf))
    {
        sec.melt.vcf[,i] <- as.character(sec.melt.vcf[,i])
    }
    
    for (i in 1:ncol(all.dta))
    {
        all.dta[,i] <- as.character(all.dta[,i])
    }
    
    sec.melt.vcf <- sec.melt.vcf[,colnames(all.dta)]
    
    #order and compare
    
    all.dta <- all.dta[do.call("order", all.dta),]
    rownames(all.dta) <- NULL
    sec.melt.vcf <- sec.melt.vcf[do.call("order", sec.melt.vcf),]
    rownames(sec.melt.vcf) <- NULL
    
    checkEquals(all.dta, sec.melt.vcf)
    
    dbDisconnect(db.con)
    file.remove("temp.db")
}

#need additional test code to ensure that only uniquely mapping probes are assigned to snps etc..
test.create.sanger.mouse.vcf.db <- function()
{
    #note that the example snps here are actually build 37
    #require(BSgenome.Mmusculus.UCSC.mm10)
    vcf.files <- om.vcf.file()
    vcf.labels <- "SNV"
    probe.tab.file <- om.tab.file()
    #seven strains besides the reference strain that are used in the CC
    strain.names <- c("CASTEiJ", "AJ", "PWKPhJ", "129S1", "NZO", "NODShiLtJ", "WSBEiJ")
    bs.genome <- BSgenome.Mmusculus.UCSC.mm9
    
    #probably shouldn't do this in real analyses unless you have to
    seqlevels(bs.genome) <- sub("chr", "", seqlevels(bs.genome))
    
    db.schema <- SangerTableSchemaList()
    db.name="test.db"
    
    keep.category <- "main"
    window.size <- 1000
    max.mismatch <- 0
    should.debug <- TRUE
    
    #duplications in this file, remove and recreate...
    tab.aln <- read.delim(om.lo.file(), sep="\t", header=TRUE, stringsAsFactors=FALSE)
    
    if (full.tests == TRUE)
    {
        limit.chr <- "19"
    }
    else
    {
        limit.chr <- GRanges(seqnames="19", ranges=IRanges(start=min(tab.aln$start), end=max(tab.aln$end)/4), strand="*")
    }
    
    create.sanger.mouse.vcf.db(vcf.files, vcf.labels, probe.tab.file, strain.names, bs.genome, db.schema, db.name, keep.category, window.size, max.mismatch, limit.chr, should.debug)
    
    db.con <- dbConnect(SQLite(), db.name)
    
    checkTrue(length(intersect(names(db.schema@tab.list), dbListTables(db.con))) == length(names(db.schema@tab.list)))
    #are the unique alignments consistent with Affy's file?
    
    probe.aln <- dbGetQuery(db.con, "SELECT * FROM probe_align NATURAL JOIN probe_info")
    #there are duplicates in this file that are discarded by the default db population approach, so will do so here...
    tab.aln <- tab.aln[!duplicated(tab.aln),]
    sub.tab.aln <- tab.aln[tab.aln$Probe.ID %in% probe.aln$probe_id,]
    
    checkTrue(nrow(probe.aln) == nrow(sub.tab.aln))
    
    merged.annot <- merge(probe.aln, sub.tab.aln, by.x="probe_id", by.y="Probe.ID", all=TRUE, incomparables=NA, sort=FALSE)
    
    checkTrue(nrow(merged.annot) == nrow(probe.aln))
    
    checkTrue(sum(merged.annot$probe_start == merged.annot$start) == nrow(merged.annot) && sum(merged.annot$probe_end == merged.annot$end) == nrow(merged.annot))
    
    #not sure I can say anything useful about the multi or non-mapped probes without using external tools...
    
    #are the genotype alignments consistent with the vcf file
    
    use.range <- with(probe.aln, GRanges(seqnames=Rle(probe_chr), ranges=IRanges(start=probe_start, end=probe_end), strand="*", fasta_name=fasta_name))
    
    sanger.vcf.param <- ScanVcfParam(trimEmpty=FALSE, fixed="ALT", geno=c("GT", "FI"), info=NA, strand="*")
    
    vcf.list <- scanVcf(vcf.files, param=sanger.vcf.param)
    
    ovls <- findOverlaps(query=vcf.list[[1]]$rowData, subject=use.range)
    
    informative.hits <- apply(vcf.list[[1]]$GENO$GT[queryHits(ovls),strain.names], 1, function(x) all(x %in% c("0/0", "./."))== FALSE)
    use.hits <- which(informative.hits)
    
    #A fix for R-2.13 as there were duplicate names...
    names(vcf.list[[1]]$rowData) <- NULL
    
    probe.snp <- cbind(as.data.frame(vcf.list[[1]]$rowData)[queryHits(ovls)[use.hits],], fasta_name=values(use.range)$fasta_name[subjectHits(ovls)[use.hits]])
    
    var.ovls <- dbGetQuery(db.con, "SELECT * FROM probe_to_snp NATURAL JOIN reference NATURAL JOIN probe_align NATURAL JOIN probe_info")
    
    #there will be some of these which are not variants in one of the 8 strains, check that here...
    temp <- merge(probe.snp, var.ovls, all=TRUE, incomparables=NA, sort=FALSE)
    
    checkTrue(nrow(temp) == nrow(probe.snp))
    checkTrue(sum(is.na(temp)) == 0)
    
    dbDisconnect(db.con)
    file.remove(db.name)
}

test.validProbeQuery<- function()
{   
    var.db <- new("VcfDB", db.path=om.db.file(), tbsl=SangerTableSchemaList(), start.var.table="reference", end.var.table="probe_info", var.mask.probe.id="probe_id", var.mask.var.id="ref_id")
    var.mask.par <- VariantMaskParams(var.db=var.db, geno.filter=FALSE, rm.unmap=TRUE, rm.mult=TRUE, mask.type="static")
    
    db.con <- dbConnect(SQLite(), var.mask.par@var.db@db.path)
    
    number.multi.un <- dbGetQuery(db.con, "SELECT COUNT(*) FROM probe_info WHERE align_status in ('MultiMapped', 'UnMapped')")
    number.un <- dbGetQuery(db.con, "SELECT COUNT(*) FROM probe_info WHERE align_status in ('UnMapped')")
    number.multi <- dbGetQuery(db.con, "SELECT COUNT(*) FROM probe_info WHERE align_status in ('MultiMapped')")
    number.unique.var <- dbGetQuery(db.con, "SELECT COUNT(DISTINCT(probe_align_id)) FROM probe_to_snp")
    number.unique.filt.var <- dbGetQuery(db.con, "SELECT COUNT(DISTINCT(probe_align_id)) FROM probe_to_snp NATURAL JOIN reference WHERE filter = 'TRUE'")
    
    target <- "probeset"
    
    #default setting
    var.mask.par@geno.filter <- FALSE
    var.mask.par@rm.mult <- TRUE
    var.mask.par@rm.unmap <- TRUE
    
    test.1 <- dbGetQuery(db.con, validProbeQuery(var.mask.par, target, should.add=FALSE))
    
    checkTrue(nrow(test.1) == number.unique.var[,1] + number.multi.un[,1])
    
    var.mask.par@geno.filter <- TRUE
    var.mask.par@rm.mult <- TRUE
    var.mask.par@rm.unmap <- TRUE
    test.2 <- dbGetQuery(db.con, validProbeQuery(var.mask.par, target, should.add=FALSE))
    
    #this should be a smaller version of test.1 as it is filtered
    checkTrue(nrow(test.2) <= nrow(test.1))
    rm(test.1)
    checkTrue(nrow(test.2) == number.unique.filt.var[,1] + number.multi.un[,1])
    rm(test.2)
    
    var.mask.par@geno.filter <- FALSE
    var.mask.par@rm.mult <- FALSE
    var.mask.par@rm.unmap <- TRUE
    
    test.3 <- dbGetQuery(db.con, validProbeQuery(var.mask.par, target, should.add=FALSE))
    
    checkTrue(nrow(test.3) == number.unique.var[,1] + number.un[,1])
    rm(test.3)
    
    var.mask.par@geno.filter <- TRUE
    var.mask.par@rm.mult <- FALSE
    var.mask.par@rm.unmap <- TRUE
    
    test.4 <- dbGetQuery(db.con, validProbeQuery(var.mask.par, target, should.add=FALSE))
    
    checkTrue(nrow(test.4) == number.unique.filt.var[,1] + number.un[,1])
    rm(test.4)
    
    var.mask.par@geno.filter <- FALSE
    var.mask.par@rm.mult <- TRUE
    var.mask.par@rm.unmap <- FALSE
    
    test.5 <- dbGetQuery(db.con, validProbeQuery(var.mask.par, target, should.add=FALSE))
    
    checkTrue(nrow(test.5) == number.unique.var[,1] + number.multi[,1])
    
    var.mask.par@geno.filter <- TRUE
    var.mask.par@rm.mult <- TRUE
    var.mask.par@rm.unmap <- FALSE
    
    test.6 <- dbGetQuery(db.con, validProbeQuery(var.mask.par, target, should.add=FALSE))
    
    checkTrue(nrow(test.6) == number.unique.filt.var[,1] + number.multi[,1])
    
    var.mask.par@geno.filter <- FALSE
    var.mask.par@rm.mult <- FALSE
    var.mask.par@rm.unmap <- FALSE
    
    test.7 <- dbGetQuery(db.con, validProbeQuery(var.mask.par, target, should.add=FALSE))
    
    checkTrue(nrow(test.7) == number.unique.var[,1])
    
    var.mask.par@geno.filter <- TRUE
    var.mask.par@rm.mult <- FALSE
    var.mask.par@rm.unmap <- FALSE
    
    test.8 <- dbGetQuery(db.con, validProbeQuery(var.mask.par, target, should.add=FALSE))
    
    checkTrue(nrow(test.8) == number.unique.filt.var[,1])
    
    dbDisconnect(db.con)
}

test.getProbeDf <- function()
{   
    var.db <- new("VcfDB", db.path=om.db.file(), tbsl=SangerTableSchemaList(), start.var.table="reference", end.var.table="probe_info", var.mask.probe.id="probe_id", var.mask.var.id="ref_id")
    var.mask.par <- VariantMaskParams(var.db=var.db)
    
    #for testing convenience dbGetQuery(db(sun.gene.fs), "detach database var_mask")
    
    db.con <- dbConnect(SQLite(), var.mask.par@var.db@db.path)
    
    #only go through all combinations if full.tests is TRUE as it gets time intenstive
    #otherwise just run the default parameters using 'core' and 'probeset' as the targets
    if (full.tests == TRUE)
    {
        possib.vals <- expand.grid(list(filter=c(TRUE, FALSE), rm.unmap=c(TRUE,FALSE), rm.multi=c(TRUE, FALSE), target=c("core", "probeset")))
        possib.vals$target <- as.character(possib.vals$target)
    }
    else
    {
        possib.vals <- data.frame(filter=rep(var.mask.par@geno.filter, 2), rm.unmap=rep(var.mask.par@rm.unmap, 2), rm.multi=rep(var.mask.par@rm.mult, 2), target=c("core", "probeset"), stringsAsFactors=FALSE)
    }
   
    
    for (i in 1:nrow(possib.vals))
    {
        print(i)
        var.mask.par@geno.filter <- possib.vals$filter[i]
        var.mask.par@rm.unmap <- possib.vals$rm.unmap[i]
        var.mask.par@rm.mult <- possib.vals$rm.multi[i]
        
        featureInfo.old <- oligo:::stArrayPmInfo(object=SunGeneFS, target=possib.vals$target[i])
        
        featureInfo.mask <- getProbeDf(object=var.mask.par, gene.fs=SunGeneFS, target=possib.vals$target[i])
        
        masked.probes <- dbGetQuery(db.con, validProbeQuery(var.mask.par, possib.vals$target[i], should.add=FALSE))
    
        featureInfo.old.rm <- featureInfo.old[featureInfo.old$fid %in% masked.probes[,var.mask.par@var.db@var.mask.probe.id] == FALSE,]
        
        featureInfo.old.rm <- featureInfo.old.rm[do.call("order", featureInfo.old.rm),]
        rownames(featureInfo.old.rm) <- NULL
        
        featureInfo.mask <- featureInfo.mask[do.call("order", featureInfo.mask),]
        rownames(featureInfo.mask) <- NULL
        
        checkEquals(featureInfo.old.rm, featureInfo.mask)
    }
    
    dbDisconnect(db.con)
}

test.maskRMA <- function()
{
    
    var.db <- new("VcfDB", db.path=om.db.file(), tbsl=SangerTableSchemaList(), start.var.table="reference", end.var.table="probe_info", var.mask.probe.id="probe_id", var.mask.var.id="ref_id")
    var.mask.par <- VariantMaskParams(var.db=var.db)
    
    #test.con <- dbConnect(SQLite(), var.db@db.path)
    #dbGetQuery(test.con, "explain query plan SELECT probe_id FROM (SELECT probe_id , COUNT( ref_id ) > 0 AS contains_var , SUM(IFNULL( filter = 'TRUE' ,0)) > 0 AS passes_filter , align_status FROM probe_info NATURAL LEFT OUTER JOIN probe_align NATURAL LEFT OUTER JOIN probe_to_snp NATURAL LEFT OUTER JOIN reference GROUP BY probe_id ) WHERE ( contains_var = 1 AND align_status = 'UniqueMapped' )")
    
    #dbDisconnect(test.con)
    
    #is the result you get from maskRMA in the default manner the same as applying medpolish in a similar manner.
    #will only test the default here as it is time intensive and other unitTests cover the retrieval and filtering procedures...
    
    if (full.tests == TRUE)
    {
        basic <- rma(SunGeneFS, background=FALSE, normalize=FALSE, target="core")
        same.as.basic <- maskRMA(SunGeneFS, background=FALSE, normalize=FALSE, target="core", apply.mask=FALSE, mask.params=var.mask.par)
        checkEquals(exprs(basic), exprs(same.as.basic))
    }
    
    
    #now with applying the mask.params
    
    applied.mask <- maskRMA(SunGeneFS, background=FALSE, normalize=FALSE, target="core", apply.mask=TRUE, mask.params=var.mask.par)
    
    raw.exprs <- exprs(SunGeneFS)
    probe.annot <- getProbeDf(object=var.mask.par, gene.fs=SunGeneFS, target="core")
    
    missing.probes <- setdiff(rownames(exprs(SunGeneFS)), probe.annot$fid)
    
    featureInfo.old <- oligo:::stArrayPmInfo(object=SunGeneFS, target="core")
    
    missing.probes.probesets <- unique(featureInfo.old$fsetid[featureInfo.old$fid %in% missing.probes])
    
    sub.probe.annot <- probe.annot[probe.annot$fsetid %in% missing.probes.probesets,]
    
    split.probes <- split(as.data.frame(raw.exprs[sub.probe.annot$fid,]), sub.probe.annot$fsetid)
    
    #from gene_st_test.R
    for (x in names(split.probes))
    {
        temp <- medpolish(log2(as.matrix(split.probes[[x]])), trace.iter=FALSE)
        
        final.exprs <- temp$overall + temp$col
        
        checkEquals(final.exprs, exprs(applied.mask)[x,])
    }
    
}