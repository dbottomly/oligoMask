
full.tests <- FALSE

stopifnot(require(RSQLite))
stopifnot(require(BSgenome.Mmusculus.UCSC.mm9))
stopifnot(require(VariantAnnotation))
stopifnot(require(reshape2))
stopifnot(require(oligo))

data("SunGeneFS")

#currently not testing the following trivial methods that are exported

#tbsl, vcfDb, om.vcf.file, om.tab.file, VariantMaskParams, maskDb
#searchTables, searchCols, searchDict


test.SangerTableSchemaList <- function()
{
    ##encountered an issue with the reference table not being unique with regards to vcf_type, so adjusted the queries and added a unit test
    
    #test data taken from the real vcf files
    all.strain.names <- c("129P2","129S1","129S5","AJ","AKRJ","BALBcJ","C3HHeJ","C57BL6NJ","CASTEiJ","CBAJ","DBA2J","FVBNJ","LPJ","NODShiLtJ","NZOHlLtJ","PWKPhJ","SPRETEiJ","WSBEiJ")
    
    snp.vcf.list <- list(`1:16511945-16511969`=list(rowData=GRanges(seqnames=Rle("1"), ranges=IRanges(start=16511962, end=16511962), strand=Rle("*")),
                                                    REF=DNAStringSet("G"),
                                                    ALT="A",
                                                    GENO=list(GT=matrix(c(rep("0/0", 16), "1/1", "."), ncol=18, dimnames=list(NULL, all.strain.names)),
                                                              FI=matrix(c(rep(1, 17), NA),ncol=18, dimnames=list(NULL, all.strain.names)))),
                         `1:72635669-72635693`=list(rowData=GRanges(seqnames=Rle("1"), ranges=IRanges(start=72635685, end=72635685), strand=Rle("*")),
                                                    REF=DNAStringSet("G"),
                                                     ALT="A",
                                                     GENO=list(GT=matrix(c(rep("0/0", 3), rep("1/1", 2), "0/0", "1/1", rep("0/0", 3), "1/1", rep("0/0", 5), rep("1/1", 2)), ncol=18, dimnames=list(NULL, all.strain.names)),
                                                               FI=matrix(rep(1, 18), ncol=18, dimnames=list(NULL, all.strain.names)))))
    
    indel.vcf.list <- list(`1:5070129-5070153`=list(rowData=GRanges(seqnames=Rle("1"), ranges=IRanges(start=5070140, end=5070175), strand=Rle("*")),
                                   REF=DNAStringSet("TCTAGGTACCCTTGGTTTTCTTGGGAAAGGCTGGGC"),
                                   ALT="T",
                                   GENO=list(GT=matrix(c(rep(".", 10), "1/1", rep(".", 7)), ncol=18, dimnames=list(NULL, all.strain.names)),
                                             FI=matrix(c(rep(NA, 10), 1, rep(NA, 7)),ncol=18, dimnames=list(NULL, all.strain.names)))),
                           `1:72635669-72635693`=list(rowData=GRanges(seqnames=Rle("1"), ranges=IRanges(start=72635685, end=72635685), strand=Rle("*")),
                                   REF=DNAStringSet("G"),
                                   ALT="GTC",
                                   GENO=list(GT=matrix(c(rep(".", 8), "1/1", rep(".", 9)), ncol=18, dimnames=list(NULL, all.strain.names)),
                                             FI=matrix(c(rep(NA, 8), 1, rep(NA, 9)),ncol=18, dimnames=list(NULL, all.strain.names)))))
    
    
    db.file <- tempfile()
    tbsl <- SangerTableSchemaList()
    
    test.db <- Database(tbsl, db.file)
    
    #pre-add the probe alignment table
    probe.align.dta <- data.frame(probe_align_id=1:3, probe_chr="1", probe_start=c(5070129, 16511945,72635669), probe_end=c(5070153, 16511969, 72635693), probe_ind=1:3, stringsAsFactors=FALSE)
    
    db.con <- dbConnect(SQLite(), db.file)
    
    dbWriteTable(db.con, "probe_align", probe.align.dta, row.names=FALSE)
    
    dbDisconnect(db.con)
    
     populate(test.db, vcf_list=snp.vcf.list, vcf_annot=c(vcf_name="test_1.txt", type="SNV"),
          use.tables=c("vcf_annot", "reference", "allele", "genotype", "probe_to_snp"), should.debug=TRUE)
    
    #snvs
    #all.snp.list <- list(vcf_list=snp.vcf.list, vcf_annot=c(vcf_name="test_1.txt", type="SNV"))
    #populate.db.tbl.schema.list(db.con, db.schema=tbsl, ins.vals=all.snp.list, use.tables=c("vcf_annot", "reference", "allele", "genotype", "probe_to_snp"), should.debug=TRUE)
    populate(test.db, vcf_list=snp.vcf.list, vcf_annot=c(vcf_name="test_1.txt", type="SNV"),
          use.tables=c("vcf_annot", "reference", "allele", "genotype", "probe_to_snp"), should.debug=TRUE)
    
    #indels
    
    #all.indel.list <- list(vcf_list=indel.vcf.list, vcf_annot=c(vcf_name="test_2.txt", type="INDEL"))
    #populate.db.tbl.schema.list(db.con, db.schema=tbsl, ins.vals=all.indel.list, use.tables=c("vcf_annot", "reference", "allele", "genotype", "probe_to_snp"), should.debug=TRUE)
    populate(test.db,vcf_list=indel.vcf.list, vcf_annot=c(vcf_name="test_2.txt", type="INDEL"),use.tables=c("vcf_annot", "reference", "allele", "genotype", "probe_to_snp"), should.debug=TRUE)
    
    db.con <- dbConnect(SQLite(), dbFile(test.db))
    
    test <- dbGetQuery(db.con, "SELECT * FROM reference NATURAL JOIN genotype NATURAL JOIN allele NATURAL JOIN vcf_annot")
    
    dbDisconnect(db.con)
    
    all.vcf <- c(snp.vcf.list, indel.vcf.list)
    
    test.target <- do.call("rbind", lapply(1:length(all.vcf), function(x)
           {
                cur.vcf <- all.vcf[[x]]
                
                if (x %in% 1:2)
                {
                    vcf_name = "test_1.txt"
                    type = "SNV"
                }
                else
                {
                    vcf_name = "test_2.txt"
                    type = "INDEL"
                }
                
                filter <- all(is.na(cur.vcf$GENO$FI) == FALSE & cur.vcf$GENO$FI == 1)
                geno_chr <- rep(1:2, length(all.strain.names))
                strain <- as.character(sapply(all.strain.names, function(y) rep(y,2)))
                allele_num <- as.character(sapply(as.character(cur.vcf$GENO$GT), function(y)
                                     {
                                        if (y != ".")
                                        {
                                            return(strsplit(y, "\\/")[[1]])
                                        }
                                        else
                                        {
                                            return(c(-1,-1))
                                        }
                                     }))
                geno.vec <- c(".", as.character(cur.vcf$REF), cur.vcf$ALT)
                alleles <- geno.vec[as.numeric(allele_num)+2]
                
                return(data.frame(seqnames=as.character(seqnames(cur.vcf$rowData)), start=start(cur.vcf$rowData), end=end(cur.vcf$rowData), filter=filter, geno_chr=geno_chr, allele_num=allele_num, strain=strain, alleles=alleles,  vcf_name=vcf_name, type=type, stringsAsFactors=FALSE))
           }))
    
    test.target.sort <- test.target[do.call("order", test.target),]
    rownames(test.target.sort) <- NULL
    
    test.target.sort$filter <- as.character(test.target.sort$filter)
    test.target.sort$allele_num <- as.numeric(test.target.sort$allele_num)
    
    test <- test[,names(test.target.sort)]
    test.sort <- test[do.call("order", test),]
    rownames(test.sort) <- NULL
    
    checkEquals(test.sort, test.target.sort)
}

###still need to fix me:  attributes are not identical across measure variables; they will be dropped

test.make.vcf.table <- function()
{
    
    #db.schema <- new("TableSchemaList")
    db.schema <- SangerTableSchemaList()
    #ust a relatively strange number < the length of probe.grange
    window.size <- 109
    vcf.name <- om.vcf.file()
    #db.con <- dbConnect(SQLite(), tempfile())
    
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
    
    db.obj <- Database(db.schema, tempfile())
    
    make.vcf.table(db.obj, window.size, vcf.name, probe.grange, vcf.type, use.tables, limit, should.debug, vcf.param, filter.func, filter.params)
    
    db.con <- dbConnect(SQLite(), dbFile(db.obj))
    
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
                temp.dta <- cbind(as.data.frame(x$rowData), ref=as.character(x$REF), alt=sapply(x$ALT, paste, collapse=","), gt=x$GENO$GT, fi=x$GENO$FI,
                                  stringsAsFactors=F)
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
    }else
    {
        limit.chr <- GRanges(seqnames="19", ranges=IRanges(start=min(tab.aln$start), end=max(tab.aln$end)/4), strand="*")
    }
    
    #first, creating the database without a package
    
    sanger.db <- new("VcfDB", db.file=db.name, tbsl=db.schema)
    
    create.sanger.mouse.vcf.db(sanger.db, vcf.files, vcf.labels, probe.tab.file, strain.names, bs.genome, keep.category, window.size, max.mismatch, limit.chr, should.debug, package.info=NULL)
    examine.vcf.db(db.name, db.schema, tab.aln, vcf.files, strain.names)
    
    #then after creating the database with a package
    
    db.name <- "test.package"
    
    sanger.db@db.file <- db.name
    
    if (file.exists(db.name))
    {
        unlink(db.name,recursive=T)
    }
    
    package.info <- list(AUTHOR="test", AUTHOREMAIL="test@test.com", BOWTIE_PATH="/Users/bottomly/Desktop/github_projects/bowtie_build/bowtie",
                         GENOME_PATH="/Users/bottomly/Desktop/resources/sequences/GRCm38_68", VCF_QUERY_CMD="htscmd vcfquery", VCF_TYPE="CC")
    
    create.sanger.mouse.vcf.db(sanger.db, vcf.files, vcf.labels, probe.tab.file, strain.names, bs.genome, keep.category, window.size, max.mismatch, limit.chr, should.debug, package.info=package.info)
    
    examine.vcf.db(file.path(db.name, "inst", "extdata", "package.db"), db.schema, tab.aln, vcf.files, strain.names)
    
    #then do an Rcheck to make sure the package builds etc...
    
    #only do this test if these files exist, which is only likely to be true on my local machine
    
    #bowtie.work <- system(package.info$BOWTIE_PATH) != 127
    #vcf.query.cmd.work <- system(package.info$VCF_QUERY_CMD) != 127
    #
    #if (bowtie.work && file.exists(paste0(package.info$GENOME_PATH, ".1.ebwt")) && vcf.query.cmd.work)
    #{
    #    #not running the tests as this is only a partial alignment example and so they will not work...
    #    did.work <- system(paste("R-3.0.1 CMD check --no-tests",db.name))
    #    checkTrue(did.work == 0)
    #}
    
    unlink(db.name, recursive=T)
}

examine.vcf.db <- function(db.name, db.schema, tab.aln, vcf.files, strain.names)
{
    db.con <- dbConnect(SQLite(), db.name)
    
    checkTrue(length(intersect(names(db.schema@tab.list), dbListTables(db.con))) == length(names(db.schema@tab.list)))
    #are the unique alignments consistent with Affy's file?
    
    probe.aln <- dbGetQuery(db.con, "SELECT * FROM probe_align NATURAL JOIN probe_info")
    #there are duplicates in this file that are discarded by the default db population approach, so will do so here...
    tab.aln <- tab.aln[!duplicated(tab.aln),]
    sub.tab.aln <- tab.aln[tab.aln$Probe.ID %in% probe.aln$probe_id,]
    
    checkTrue(nrow(probe.aln) == nrow(sub.tab.aln))
    
    merged.annot <- merge(probe.aln, sub.tab.aln, by.x="probe_id", by.y="Probe.ID", all=TRUE, sort=FALSE)
    
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
    temp <- merge(probe.snp, var.ovls, all=TRUE,  sort=FALSE)
    
    checkTrue(nrow(temp) == nrow(probe.snp))
    checkTrue(sum(is.na(temp)) == 0)
    
    dbDisconnect(db.con)
    file.remove(db.name)
}

test.getProbeDf <- function()
{
    DEACTIVATED("Need to revise logic behind this test")
    
    var.db <- new("VcfDB", db.file=om.db.file())
    var.mask.par <- VariantMaskParams(var.db=var.db)
    
    #for testing convenience dbGetQuery(db(sun.gene.fs), "detach database var_mask")
    
    db.con <- dbConnect(SQLite(), var.mask.par@var.db@db.file)
    
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
    
    var.db <- new("VcfDB", db.file=om.db.file())
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
    
    #try again using the "before.summary" option in a basic sense
    
    applied.mask.2 <- maskRMA(SunGeneFS, background=TRUE, normalize=TRUE, target="core", mask.type="before.summary", apply.mask=TRUE, mask.params=var.mask.par)
    
    init.exprs <- exprs(SunGeneFS)
    
    init.exprs.bg <- preprocessCore:::rma.background.correct(init.exprs)
    
    init.exprs.norm <- preprocessCore:::normalize.quantiles(init.exprs.bg)
    
    checkTrue(all(dim(init.exprs.norm) == dim(init.exprs)))
    
    rownames(init.exprs.norm) <- rownames(init.exprs)
    colnames(init.exprs.norm) <- colnames(init.exprs)
    
    split.probes.2 <- split(as.data.frame(init.exprs.norm[sub.probe.annot$fid,]), sub.probe.annot$fsetid)
    
    #from gene_st_test.R
    for (x in names(split.probes.2))
    {
        temp <- medpolish(log2(as.matrix(split.probes.2[[x]])), trace.iter=FALSE)
        
        final.exprs <- temp$overall + temp$col
        
        checkEquals(final.exprs, exprs(applied.mask.2)[x,])
    }
    
}