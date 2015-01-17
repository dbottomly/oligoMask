#guess.core.query <- function()
#{
#    guess.query.func(oligo:::getFidMetaProbesetCore)
#}
#
#guess.probeset.query <- function()
#{
#    guess.query.func(oligo:::getFidProbeset)
#}

#guess.query.func <- function(object, target)
#{
#    if (target == "core")
#    {
#        use.func <- oligo:::getFidMetaProbesetCore
#    }else
#    {
#        use.func <- oligo:::getFidProbeset
#    }
#    
#    func.list <- as.list(body(use.func))
#    
#    which.sql <- sapply(func.list, function(x) any(grep("SELECT", x)))
#    
#    if (sum(which.sql) != 1)
#    {
#        stop("ERROR: Guess at query failed--multiple possible queries found")
#    }
#    
#    sql.val <- eval(func.list[which.sql][[1]])
#    
#    if (is.character(sql.val) == F || length(sql.val) != 1)
#    {
#        stop("ERROR: Guess at query failed--only one SQL query expected")
#    }
#    
#    return(sql.val)
#}

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

#should be replaced with poplite functionality...
#get.shortest.query.path <- function(var.mask.par, start=NULL, finish=NULL, reverse=TRUE, undirected=TRUE)
#{
#    if (class(var.mask.par) != "VariantMaskParams")
#    {
#        stop("ERROR: var.mask.par needs to be of class VariantMaskParams")
#    }
#    
#    if (missing(start) || is.null(start) || is.na(start))
#    {
#        start <- var.mask.par@var.db@start.var.table
#    }
#    
#    if (missing(finish) || is.null(finish) || is.na(finish))
#    {
#        finish <- var.mask.par@var.db@end.var.table
#    }
#    
#    tsl.graph <- tsl.to.graphNEL(var.mask.par@var.db@tbsl)
#    
#    if (undirected)
#    {
#        tsl.graph <- ugraph(tsl.graph)
#    }
#    
#    table.path <- sp.between(g=tsl.graph,start=start,finish=finish, detail=TRUE)
#    
#    #do it backwards as that is how we want to merge it with the pd tables
#    if (reverse==TRUE)
#    {
#        mask.query <- table.path[[1]]$path_detail[length(table.path[[1]]$path_detail):1]
#    }
#    else
#    {
#        mask.query <- table.path[[1]]$path_detail
#    }
#    
#    return(mask.query)
#}



#choose the set of probe IDs to remove from consideration
#setGeneric("validProbeQuery", def=function(object,...) standardGeneric("validProbeQuery"))
#setMethod("validProbeQuery", signature("VariantMaskParams"), function(object, target, should.add=TRUE, should.count=FALSE)
#          {
#                oligo.query <- guess.query.func(object, target)
#                
#                var.presence <- "contains_var"
#                filter.presence <- "passes_filter"
#                
#                #build a query in the snp mask by chasing down the foreign keys in a naive way starting at the specified starting table
#                #this first creates a directed graph representation and chases down the dependencies using the adjacent nodes
#                
#                #query.tables <- get.shortest.query.path(object)
#                
#                #query.vec <- check.and.add.tables(table.name=searchTables(object@var.db), query.tables=query.tables, object=object, default.join.type="NATURAL LEFT OUTER JOIN")
#                
#                #poplite:::join()
#                
#                #as there is always the possibility of multi variants per probe, the query will always have to aggregate to the probe level
#                
#                inner.query <- paste("(SELECT",object@var.db@var.mask.probe.id,", COUNT(",object@var.db@var.mask.var.id ,") > 0 AS",var.presence ,
#                                  ", SUM(IFNULL(", searchCols(object@var.db, "genotype.filter"),dict.to.where(object, dict.name="genotype.filter", values="TRUE"), ",0)) > 0 AS", filter.presence, ",",
#                                  searchCols(obj=object@var.db, name="mapping.status"), "FROM", names(query.vec)[1], paste(paste(query.vec[2:length(query.vec)], names(query.vec)[2:length(query.vec)]), collapse=" "), "GROUP BY", object@var.db@var.mask.probe.id, ")")
#                
#                #either way remove those that are uniquely mapping and that have a variant overlapping them
#                outer.query <- paste("SELECT", object@var.db@var.mask.probe.id, "FROM", inner.query)
#                
#                if (object@geno.filter == TRUE)
#                {
#                    #filter the result to only those with a filter value of TRUE, for this a NATURAL JOIN would suffice as we are filtering down
#                    
#                   where.base <- paste("WHERE (", var.presence, "= 1 AND", filter.presence, "= 1 AND", searchCols(obj=object@var.db, name="mapping.status"), dict.to.where(object, "mapping.status", "unique"), ")")
#                }
#                else
#                {
#                    where.base <- paste("WHERE (", var.presence, "= 1 AND", searchCols(obj=object@var.db, name="mapping.status"), dict.to.where(object, "mapping.status", "unique"), ")")
#                }
#                
#                outer.query <- paste(outer.query, where.base)
#                
#                if ((object@rm.unmap == TRUE || object@rm.mult == TRUE))
#                {
#                    if ((object@rm.unmap == TRUE && object@rm.mult == TRUE))
#                    {
#                        where.suff <- dict.to.where(object, "mapping.status", c("multi", "non"))
#                    }
#                    else if (object@rm.unmap == FALSE && object@rm.mult == TRUE)
#                    {
#                        where.suff <- dict.to.where(object, "mapping.status", "multi")
#                    }
#                    else if (object@rm.unmap == TRUE && object@rm.mult == FALSE)
#                    {
#                        where.suff <- dict.to.where(object, "mapping.status", "non")
#                    }
#                    else
#                    {
#                        stop("ERROR: rm.unmap and rm.mult should not be false here")
#                    }
#                    
#                    outer.query <- paste(outer.query, "OR (", searchCols(obj=object@var.db, name="mapping.status"), where.suff, ")")
#                }
#                
#                #if none are to be removed then don't add to the query 
#                
#                if (should.add == TRUE)
#                {
#                    ex.query <- paste(oligo.query, "JOIN (", outer.query, ") ON", object@oligo.probe.id, "=", object@var.db@var.mask.probe.id)
#                    
#                    return(paste(oligo.query, "EXCEPT", ex.query))
#                }
#                else
#                {
#                    return(outer.query)
#                }
#          })

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

#tsl.to.graphNEL <- function(tsl)
#{
#    edge.l <- lapply(tsl@tab.list, function(x)
#           {
#                cur.edges <- unlist(lapply(x$foreign.keys, "[[", "local.keys"))
#                if (is.null(cur.edges))
#                {
#                    return(list(edges=NULL))
#                }
#                else
#                {
#                    common.col <- intersect(x$db.cols, cur.edges)
#                    use.fk <- sapply(x$foreign.keys, function(x) x$local.keys %in% common.col)
#                    return(list(edges=names(use.fk)[use.fk == TRUE]))
#                }
#                    
#                #the old way...
#                #return(list(edges=names(x$foreign.keys)))
#           })
#    return(graphNEL(nodes=names(edge.l), edgeL=edge.l, edgemode='directed'))
#}

#not sure about this function, may no longer be necessary depending on poplite....
#check.and.add.tables <- function(table.name, query.tables, object, default.join.type="NATURAL JOIN")
#{
#    diff.tabs <- setdiff(table.name, query.tables)
#    
#    if (length(diff.tabs) > 0)
#    {
#        for (i in diff.tabs)
#        {
#            #attempt to add to query.tables by finding the shortest path between the end table and the specified table.name
#            new.query.path <- get.shortest.query.path(var.mask.par=object, start=i, finish=query.tables[length(query.tables)])
#            non.exist.tabs <- setdiff(new.query.path, query.tables)
#            query.tables <- append(query.tables, non.exist.tabs)
#        }
#    }
#    
#    query.joins <- rep(default.join.type, length(query.tables))
#    names(query.joins) <- query.tables
#    return(query.joins)
#}

#what the method is now...
#object="GeneFeatureSet"
#function (object, ...) 
#{
#    .local <- function (object, background = TRUE, normalize = TRUE, 
#        subset = NULL, target = "core") 
#    {
#        target <- match.arg(target, c("core", "probeset"))
#        featureInfo <- stArrayPmInfo(object, target = target)
#        theClass <- class(exprs(object))
#        pmi <- featureInfo[["fid"]]
#        pnVec <- as.character(featureInfo[["fsetid"]])
#        if ("matrix" %in% theClass) {
#            pms <- exprs(object)[pmi, , drop = FALSE]
#            dimnames(pms) <- NULL
#            colnames(pms) <- sampleNames(object)
#            theExprs <- basicRMA(pms, pnVec, normalize, background)
#            rm(pms)
#        }

#with:

#oligo:::stArrayPmInfo
#function (object, target = "core", sortBy = "fsetid") 
#{
#    target <- match.arg(target, c("probeset", "core", "full", 
#        "extended"))
#    theFun <- switch(target, probeset = getFidProbeset, core = getFidMetaProbesetCore, 
#        full = getFidMetaProbesetFull, extended = getFidMetaProbesetExtended)
#    theFun(object, sortBy = sortBy)
#}


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

#setMethod("getProbeDf", signature("VariantMaskParams"), function(object, gene.fs, target, sortBy="fsetid")
#          {
#                #check to see if the database already looks attached...
#                
#                all.tabs.pres <- sapply(names(object@var.db@tbsl@tab.list), function(x)
#                       {
#                            return(dbExistsTable(db(gene.fs), x))
#                       })    
#            
#                if (all(all.tabs.pres))
#                {
#                    dbGetQuery(db(gene.fs), paste0("DETACH DATABASE 'var_mask'"))
#                }
#                
#                #attach the snp mask database to the pd database and retrieve a data.frame of the form found in the rma method below
#                dbGetQuery(db(gene.fs), paste0("ATTACH DATABASE '",dbFile(object@var.db),"' AS var_mask"))
#                
#                fs.dta <- dbGetQuery(db(gene.fs), validProbeQuery(object, target))
#                
#                if (!is.null(sortBy)) {
#                fs.dta <- fs.dta[order(fs.dta[,sortBy]),]
#                rownames(fs.dta) <- NULL
#                }
#                
#                dbGetQuery(db(gene.fs), paste0("DETACH DATABASE 'var_mask'"))
#                
#                return(fs.dta)
#          })