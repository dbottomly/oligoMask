###taken from the PD package definitions to support the test data

##need to uncomment me when refactoring is done...
globals <- new.env(hash=TRUE, parent=emptyenv())
globals$DEBUG <- FALSE

globals$DB_PATH <- system.file("extdata", "mogene_1_0_probe_chr19.db",
                               package="oligoMask")
if (nchar(globals$DB_PATH) == 0)
  stop("Unable to locate DB file")

initDbConnection <- function() {
    globals$dbCon <- dbConnect(SQLite(), dbname=globals$DB_PATH)
    globals$dbCon
}

getDb  <- function() {
    if (!is.null(globals$dbCon) && isIdCurrent(globals$dbCon))
      return(globals$dbCon)
    initDbConnection()
}

closeDb <- function() {
    ## FIXME: check for valid connection?
    if (isIdCurrent(globals$dbCon)){
        sapply(dbListResults(globals$dbCon), dbClearResult)
        dbDisconnect(globals$dbCon)
    }
    remove(list="dbCon", envir=globals)
}

.onAttach <- function(libname, pkgname) {
    globals$DB_PATH <- system.file("extdata", "mogene_1_0_probe_chr19.db",
                                   package="oligoMask",
                                   lib.loc=libname)
    if (nchar(globals$DB_PATH) == 0)
      stop("Unable to locate DB file")
    ## Establish a connection to the SQLite DB
    initDbConnection()
}

.onUnload <- function(libpath) {
    closeDb()
}

test.pd.1.0.st.v1 <- new("AffyGenePDInfo",
                    genomebuild="mm09",
                    getdb=getDb,
                    geometry=as.integer(strsplit("1050;1050", ";")[[1]]),
                    annotation="test.pd.1.0.st.v1")

setClass("GeneFeatureSetTest", contains=c("FeatureSet", "GeneFeatureSet"))

setMethod("show", "GeneFeatureSetTest", function(object){
    print("This is test data")
})

##end taken code

om.vcf.file <- function() system.file(file.path("extdata", "chr19_snvs_sort.vcf.gz"),package="oligoMask") 
om.tab.file <- function() system.file(file.path("extdata", "mogene_1_0_probe_chr19.tab.gz"), package="oligoMask")
om.db.file <- function() system.file(file.path("extdata", "NOD_B6_m37_chr19.db"), package="oligoMask")
om.lo.file <- function() system.file(file.path("extdata", "mogene_1_0_probe_chr19.tab.liftover.txt.gz"), package="oligoMask")

make.ref.dta <- function(vcf.list)
{
    annots <- vcf.list[["vcf_annot"]]
    annot.dta <- matrix(annots, nrow=1, dimnames=list(NULL, names(annots)))
    ref.dta <- do.call("rbind", lapply(vcf.list[["vcf_list"]], function(x)
                                       {
                                            vcf.dta <- as.data.frame(x$rowData)[,c("seqnames", "start", "end")]
                                            if (any(is.na(x$GENO$FI)) == FALSE && any(x$GENO$FI %in% c(0,1) == FALSE))
                                            {
                                                print(x$GENO$FI)
                                                stop("ERROR: Unexpected FI values found")
                                            }
                                            #also redefine the filter column to indicate whether all strain genotypes passed the filter or not
                                            x$GENO$FI[is.na(x$GENO$FI)] <- 0
                                            all.strains.pass <- apply(x$GENO$FI, 1, function(y) as.character(all(y==1)))
                                            vcf.dta$filter <- all.strains.pass
                                            return(vcf.dta)
                                       }))
    
    return(cbind(ref.dta, annot.dta[rep(1, nrow(ref.dta)),,drop=FALSE]))
}

make.allele.dta <- function(vcf.list)
{
    annots <- vcf.list[["vcf_annot"]]
    annot.dta <- matrix(annots, nrow=1, dimnames=list(NULL, names(annots)))
                               
    #Note that the ALT element was changed from comma-delimited to list-based in the current code in bioc 2.13, use the below hack for now which should be backwards compatible 
    alleles.list <- lapply(vcf.list[["vcf_list"]], function(x)
           {
                return(cbind(as.data.frame(x$rowData)[,c("seqnames", "start", "end")], REF=as.character(x$REF), ALT=sapply(x$ALT, paste, collapse=","), stringsAsFactors=FALSE))
           })
    
    allele.dta <- data.frame(do.call("rbind", alleles.list), stringsAsFactors=FALSE)
    
    use.allele.dta <- data.frame(do.call("rbind", lapply(1:nrow(allele.dta), function(x)
               {
                alt.alleles <- strsplit(as.character(allele.dta$ALT[x]), ",")[[1]]
                ref.alleles <- as.character(allele.dta$REF[x])
                
                alleles <- c(".", ref.alleles, alt.alleles)
                allele_num <- -1:length(alt.alleles)
                
                stopifnot(length(alleles) == length(allele_num))
                
                return(cbind(allele.dta[rep(x, length(alleles)),c("seqnames", "start", "end")], alleles=alleles, allele_num=allele_num, stringsAsFactors=FALSE))
               })), stringsAsFactors=FALSE)
    
    return(cbind(use.allele.dta, annot.dta[rep(1, nrow(use.allele.dta)),,drop=FALSE]))
}

make.genotype.dta <- function(vcf.list)
{
    annots <- vcf.list[["vcf_annot"]]
    annot.dta <- matrix(annots, nrow=1, dimnames=list(NULL, names(annots)))
    
    geno.list <- lapply(vcf.list[["vcf_list"]], function(x)
           {
                return(cbind(as.data.frame(x$rowData)[,c("seqnames", "start", "end")], GT=x$GENO$GT, stringsAsFactors=FALSE))
           })
    geno.dta <- data.frame(do.call("rbind", geno.list), stringsAsFactors=FALSE)
    
    melt.geno <- melt(geno.dta, id.vars=c("seqnames", "start", "end"))
    names(melt.geno) <- c("seqnames","start","end", "strain", "GT")
    melt.geno$strain <- sub("GT\\.", "", melt.geno$strain)
    melt.geno$GT[is.na(melt.geno$GT) | melt.geno$GT %in% c(".", "./.") ] <- "-1/-1"

    split.geno <- do.call("rbind", strsplit(as.character(melt.geno$GT), "/"))
    melt.geno <- cbind(melt.geno, chr_num=split.geno)
    
    use.geno <- melt(melt.geno, measure.vars=c("chr_num.1","chr_num.2"))
    use.geno <- use.geno[,-which(names(use.geno) == "GT")]
    names(use.geno)[5:6] <- c("geno_chr", "allele_num")    
    use.geno$geno_chr <- sub("chr_num\\.", "", use.geno$geno_chr)
    
    return(cbind(use.geno, annot.dta[rep(1, nrow(use.geno)),,drop=FALSE]))
}

granges.to.dta <- function(x)
    {
        #only keep the uniquely aligned reads for this purpose
        unique.reads <- x[["probe_info"]]$fasta_name[x[["probe_info"]]$align_status == "UniqueMapped"]
        x[["probe_align"]] <- x[["probe_align"]][values(x[["probe_align"]])$probe.name %in% unique.reads]
        
        temp.x <- as.data.frame(x[["probe_align"]])
        temp.x$fasta_name <- rownames(temp.x)
        sub.temp.x <- temp.x[,c("seqnames", "start", "end", "fasta_name")]
        names(sub.temp.x) <- c("probe_chr", "probe_start", "probe_end", "fasta_name")
        rownames(sub.temp.x) <- NULL
        return(sub.temp.x)
    }

make.probe.to.snp <- function(vcf.list)
{
    annots <- vcf.list[["vcf_annot"]]
    annot.dta <- matrix(annots, nrow=1, dimnames=list(NULL, names(annots)))
    
    vcf.list <- vcf.list[["vcf_list"]]
    
    p.s.list <- lapply(names(vcf.list), function(x)
           {
                probe.coords <- as.data.frame(matrix(strsplit(x, ":|-")[[1]], nrow=1, dimnames=list(NULL, c("probe_chr", "probe_start", "probe_end"))))
                ps.dta <- cbind(as.data.frame(vcf.list[[x]]$rowData), probe.coords[rep(1, length(vcf.list[[x]]$rowData)),])
                sub.ps.dta <- ps.dta[,c("seqnames", "start", "end", "probe_chr", "probe_start", "probe_end")]
                return(sub.ps.dta)
           })
    
    res.ps.dta <- data.frame(do.call("rbind", p.s.list), stringsAsFactors=FALSE)
    
    return(cbind(res.ps.dta, annot.dta[rep(1, nrow(res.ps.dta)),,drop=FALSE]))
}