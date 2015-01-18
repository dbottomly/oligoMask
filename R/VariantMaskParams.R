valid.VariantMaskParams <- function(object)
{
    #length(rm.unmap should be 1 and should be either TRUE or FALSE)
    #similar for rm.mult
    #geno.filt should be of length 1 and either TRUE or FALSE
    
    #tsl should contain search columns corresponding to 'mapping.status' and 'filter.status'
    
    if (object@rm.unmap %in% c(TRUE, FALSE) == FALSE)
    {
        return(FALSE)
    }
    
    if (object@rm.mult %in% c(TRUE, FALSE) == FALSE)
    {
        return(FALSE)
    }
    
    if (object@geno.filter %in% c(TRUE, FALSE) == FALSE)
    {
        return(FALSE)
    }
    
    if (object@mask.type != "static")
    {
        return(FALSE)
    }
    
    return(TRUE)
}

setClass(Class="VariantMaskParams", representation=list(mask.type="character", geno.filter="logical", rm.unmap="logical", rm.mult="logical", 
                                                        oligo.probe.id="character", var.db="VcfDB"),
                                    prototype=prototype(mask.type="static", geno.filter=FALSE, rm.unmap=TRUE, rm.mult=TRUE, oligo.probe.id ="fid", var.db=new("VcfDB")),
                                    validity=valid.VariantMaskParams)

setMethod("show", signature("VariantMaskParams"), function(object)
          {
                message("An object of class VariantMaskParams")
          })

setGeneric("maskDb", def=function(obj, ...) standardGeneric("maskDb"))
setMethod("maskDb", signature("VariantMaskParams"), function(obj)
          {
                return(dbFile(obj@var.db))
          })

VariantMaskParams <- function(var.db, geno.filter=FALSE, rm.unmap=TRUE, rm.mult=TRUE, mask.type="static")
{
    if (class(var.db) != "VcfDB")
    {
        stop("ERROR: var.db needs to be of class VcfDB")
    }
    
    return(new("VariantMaskParams", var.db=var.db, geno.filter=geno.filter, rm.unmap=rm.unmap, rm.mult=rm.mult, mask.type=mask.type))
}

setMethod("searchTables", signature("VariantMaskParams"), function(obj, name)
          {
            searchTables(obj@var.db, name=name)
          })

setMethod("searchCols", signature("VariantMaskParams"), function(obj, name, include.table=F)
          {
            searchCols(obj@var.db, name=name, include.table=include.table)
          })

setMethod("searchDict", signature("VariantMaskParams"), function(obj, name, value=NULL)
          {
            searchDict(obj@var.db, name=name, value=value)
          })

dict.to.where <- function(object, dict.name, values)
{
    search.vals <- searchDict(obj=object, name=dict.name, value=values)
    where.vals <- ifelse(sapply(search.vals, is.character), paste0("'", search.vals, "'"), search.vals)
    
    if (length(where.vals) > 1)
    {
        where.suff <- paste0("%in% ", "c(",paste(where.vals, collapse=",") , ")")
    }
    else if (length(where.vals) == 1)
    {
        where.suff <- paste0("== ", where.vals)
    }
    else
    {
        stop("ERROR: Unexpected length found for where.vals")
    }
    
    return(where.suff)
}

setGeneric("maskRMA", def=function(object, ...) standardGeneric("maskRMA"))
setMethod("maskRMA", signature("GeneFeatureSet"), function(object, background=TRUE, normalize=TRUE, subset=NULL, target="core", mask.type=c("before.rma", "before.summary"), apply.mask=FALSE, mask.params=NULL)
          {
            if (apply.mask == FALSE)
            {
                return(rma(object, background, normalize, subset, target))
            }
            else
            {
                if (missing(mask.params) || is.null(mask.params))
                {
                    mask.params <- new("VariantMaskParams")
                }
                else if (class(mask.params) != "VariantMaskParams")
                {
                    stop("ERROR: mask.params needs to be an object of class VariantMaskParams")
                }
                
                target <- match.arg(target, c("core", "probeset"))
                mask.type <- match.arg(mask.type)
             
                #need to get a data.frame with columns fid and fsetid corresponding to the probes and probesets/metaprobesets to be used
                featureInfo <- getProbeDf(mask.params, object, target)
                
                theClass <- class(exprs(object))
                pmi <- featureInfo[["fid"]]
                pnVec <- as.character(featureInfo[["fsetid"]])
                
                if ("matrix" %in% theClass) {
                    
                    if (mask.type == "before.rma")
                    {
                        pms <- exprs(object)[pmi, , drop = FALSE]
                        dimnames(pms) <- NULL
                        colnames(pms) <- sampleNames(object)
                        theExprs <- basicRMA(pms, pnVec, normalize, background)
                        rm(pms)
                    }
                    else if (mask.type == "before.summary")
                    {
                        if (background == TRUE)
                        {
                            bg.mat <- backgroundCorrect(exprs(object), method="rma", target=target)
                        }
                        else
                        {
                            bg.mat <- exprs(object)
                        }
                        
                        if (normalize == TRUE)
                        {
                            n.mat <- normalize(bg.mat, method="quantile", target=target)
                        }
                        else
                        {
                            n.mat <- bg.mat
                        }
                        
                        dimnames(n.mat) <- dimnames(exprs(object))
                        pms <- n.mat[pmi,,drop=FALSE]
                        dimnames(pms) <- NULL
                        colnames(pms) <- sampleNames(object)
                        theExprs <- summarize(pms, probes=pnVec, method="medianpolish")
                    }
                    else
                    {
                        stop("ERROR: Unexpected value for mask.type")
                    }
                    
                }
                else if ("ff_matrix" %in% theClass) {
                   stop("ERROR: oligoMask support for type 'ff_matrix' not currently supported")
                }
                else {
                    stop("basicRMA not implemented for '", theClass,"' objects.")
                }
                out <- new("ExpressionSet", assayData=assayDataNew(exprs = theExprs), phenoData=phenoData(object), featureData=oligo:::basicAnnotatedDataFrame(theExprs, byrow = TRUE),
                           protocolData=protocolData(object), annotation=object@annotation)

                if (validObject(out)) {
                    return(out)
                }
                else {
                    stop("Resulting object is invalid.")
                }
            }
          })

setGeneric("getProbeDf", def=function(object, ...) standardGeneric("getProbeDf"))
setMethod("getProbeDf", signature("VariantMaskParams"), function(object, gene.fs, target, sortBy="fsetid"){
    
    var.presence <- "contains_var"
    filter.presence <- "passes_filter"
    
    #get the probesets from oligoMask
    
    #get data corresponding to the requested columns
    
    #something like below
    #all.db <- select(test.db, probe_id, align_status, reference.ref_id, filter)
    
    all.db <- select_(object@var.db, object@var.db@var.mask.probe.id,
                      object@var.db@var.mask.var.id,
                      setNames(searchCols(object,  name="genotype.filter", include.table=T), NULL),
                      setNames(searchCols(object,  name="mapping.status", include.table=T), NULL))
    
    var.only <- sapply(strsplit(object@var.db@var.mask.var.id, "\\."), "[", 2)
    probe.only <- sapply(strsplit(object@var.db@var.mask.probe.id, "\\."), "[", 2)
    
    geno.counts <- paste("n(", var.only ,") > 0")

    geno.filt <- paste("sum(IFNULL(",searchCols(object, "genotype.filter"),dict.to.where(object, dict.name="genotype.filter", values="TRUE"), ",0)) > 0")
    
    probe.group <- group_by_(all.db, probe.only, searchCols(obj=object, name="mapping.status"))
    
    sum.by.probe <- summarize_(probe.group, .dots=setNames(list(geno.counts, geno.filt), c(var.presence, filter.presence)))
    
    #then add in the actual filters as below:
    
    if (object@geno.filter == TRUE)
    {
                    
        #where.base <- paste("WHERE (", var.presence, "= 1 AND", filter.presence, "= 1 AND", searchCols(obj=object@var.db, name="mapping.status"), dict.to.where(object, "mapping.status", "unique"), ")")
        probe.sum.filt <- filter_(sum.by.probe, paste(paste(var.presence, "== 1"), "&",
                                             paste(filter.presence, "== 1"), "&",
                                             paste(searchCols(obj=object, name="mapping.status"), dict.to.where(object, "mapping.status", "unique"))
                                             ))
    }
    else
    {
        #where.base <- paste("WHERE (", var.presence, "= 1 AND", searchCols(obj=object@var.db, name="mapping.status"), dict.to.where(object, "mapping.status", "unique"), ")")
        probe.sum.filt <- filter_(sum.by.probe, paste(paste(var.presence, "== 1"), "&",
                                             paste(searchCols(obj=object, name="mapping.status"), dict.to.where(object, "mapping.status", "unique"))
                                             ))
    }
    
    if ((object@rm.unmap == TRUE || object@rm.mult == TRUE))
    {
        if ((object@rm.unmap == TRUE && object@rm.mult == TRUE))
        {
            #where.suff <- dict.to.where(object, "mapping.status", c("multi", "non"))
            mapping.query <- paste( searchCols(obj=object, name="mapping.status"), dict.to.where(object, "mapping.status", c("multi", "non")))
        }
        else if (object@rm.unmap == FALSE && object@rm.mult == TRUE)
        {
            #where.suff <- dict.to.where(object, "mapping.status", "multi")
            mapping.query <- paste( searchCols(obj=object, name="mapping.status"), dict.to.where(object, "mapping.status", "multi"))
        }
        else if (object@rm.unmap == TRUE && object@rm.mult == FALSE)
        {
            #where.suff <- dict.to.where(object, "mapping.status", "non")
            mapping.query<- paste( searchCols(obj=object, name="mapping.status"), dict.to.where(object, "mapping.status", "non"))
        }
        else
        {
            stop("ERROR: rm.unmap and rm.mult should not be false here")
        }
        
        all.rm.probes <- dplyr::union(select_(probe.sum.filt, probe.only), select_(filter_(select(object@var.db, .tables="probe_info"), mapping.query), probe.only), copy=T)
        
    }else{
        
        all.rm.probes <- select_(probe.sum.filt, probe.only)
    }

    #get the probesets from oligo and remove the ones from above
    
    oligo.probe.dta <- oligo:::stArrayPmInfo(object = gene.fs, target = target, sortBy = "fsetid")
    
    return(oligo.probe.dta[as.character(oligo.probe.dta[,object@oligo.probe.id]) %in% as.character(as.data.frame(all.rm.probes)[,probe.only]) == F,])
})