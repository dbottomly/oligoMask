#Note that probe_to_snp has should.ignore set to TRUE as there are some completely overlapping probes
#ie dbGetQuery(db.con, "select * from probe_align where probe_start = 44230225") from test.make.vcf.table

default.tab.list <- function()
{
    return(list(probe_info=list(db.cols=c("probe_ind", "fasta_name", "probe_id", "align_status"),
                                db.schema=c("INTEGER PRIMARY KEY AUTOINCREMENT", "TEXT", "INTEGER", "TEXT"),
                                db.constr="CONSTRAINT probe_idx UNIQUE (fasta_name)",
                                dta.func=function(x)return(x[["probe_info"]]), should.ignore=FALSE, foreign.keys=NULL),
                probe_align=list(db.cols=c("probe_align_id", "probe_chr", "probe_start", "probe_end", "probe_ind"),
                                db.schema=c("INTEGER PRIMARY KEY AUTOINCREMENT", "TEXT", "INTEGER", "INTEGER", "INTEGER"),
                                db.constr="CONSTRAINT geno_idx UNIQUE (probe_chr, probe_start, probe_end, probe_ind)",
                                dta.func=granges.to.dta, should.ignore=FALSE, foreign.keys=list(probe_info=list(local.keys="probe_ind", ext.keys="fasta_name"))),
                vcf_annot=list(db.cols=c("vcf_annot_id", "vcf_name", "type"),
                               db.schema=c("INTEGER PRIMARY KEY AUTOINCREMENT", "TEXT", "TEXT"),
                               db.constr="CONSTRAINT probe_idx UNIQUE (vcf_name, type)",
                               dta.func=function(x) data.frame(vcf_name=x[["vcf_annot"]]["vcf_name"], type=x[["vcf_annot"]]["type"], stringsAsFactors=FALSE),
                               should.ignore=TRUE, foreign.keys=NULL),
                reference=list(db.cols=c("ref_id", "seqnames", "start", "end", "filter", "vcf_annot_id"),
                               db.schema=c("INTEGER PRIMARY KEY AUTOINCREMENT", "TEXT", "INTEGER", "INTEGER", "TEXT", "INTEGER"),
                               db.constr="CONSTRAINT ref_idx UNIQUE (seqnames, start, end)",
                               dta.func=make.ref.dta, should.ignore=TRUE, foreign.keys=list(vcf_annot=list(local.keys="vcf_annot_id", ext.keys=c("vcf_name", "type")))),
                allele=list(db.cols=c("allele_id", "alleles", "allele_num", "ref_id"),
                            db.schema=c("INTEGER PRIMARY KEY AUTOINCREMENT", "TEXT", "INTEGER", "INTEGER"),
                            db.constr="CONSTRAINT alelle_idx UNIQUE (alleles, allele_num, ref_id)",
                            dta.func=make.allele.dta, should.ignore=TRUE, foreign.keys=list(reference=list(local.keys="ref_id", ext.keys=c("seqnames", "start", "end")))),
                genotype=list(db.cols=c("geno_id", "geno_chr", "allele_num","strain", "ref_id"),
                              db.schema=c("INTEGER PRIMARY KEY AUTOINCREMENT", "INTEGER", "INTEGER", "TEXT", "INTEGER"),
                              db.constr="CONSTRAINT geno_idx UNIQUE (ref_id, strain, geno_chr, allele_num)",
                              dta.func=make.genotype.dta, should.ignore=TRUE, foreign.keys=list(reference=list(local.keys="ref_id", ext.keys=c("seqnames", "start", "end")))),
                probe_to_snp=list(db.cols=c("probe_snp_id", "ref_id", "probe_align_id"),
                                  db.schema=c("INTEGER PRIMARY KEY AUTOINCREMENT", "INTEGER", "INTEGER"),
                                  db.constr="CONSTRAINT p_s_idx UNIQUE (ref_id, probe_align_id)",
                                  dta.func=make.probe.to.snp, should.ignore=TRUE, foreign.keys=list(reference=list(local.keys="ref_id", ext.keys=c("seqnames", "start", "end")),
                                                                                    probe_align=list(local.keys="probe_align_id", ext.keys=c("probe_chr", "probe_start", "probe_end"))))))
}

default.search.cols <- function()
{
    return(list(mapping.status=list(table="probe_info", column="align_status", dict=c(unique="UniqueMapped", multi="MultiMapped", non="UnMapped")),
                genotype.filter=list(table="reference", column="filter", dict=c(`TRUE`="TRUE", `FALSE`="FALSE"))))
}

#need to add validitity checks to default.search.cols to below function
valid.TableSchemaList <- function(object)
{
    if (is.null(names(object@tab.list)))
    {
        return("The supplied tab.list needs to have names")
    }
    
    valid.list <- sapply(object@tab.list, function(x)
           {
                if (is.null(names(x)) == TRUE || all(names(x) %in% c( "db.cols","db.schema", "db.constr", "dta.func", "should.ignore", "foreign.keys")) == FALSE)
                {
                    return(FALSE)
                }
                else
                {
                    if (length(x$db.schema) == length(x$db.cols))
                    {
                        if (is.null(x$foreign.keys))
                        {
                            return(TRUE)
                        }
                        else if (class(x$foreign.keys) == "list" && is.null(names(x$foreign.keys)) == FALSE)
                        {
                            return(all(sapply(x$foreign.keys, function(x)
                                   {
                                        return(all(names(x) %in% c("local.keys", "ext.keys")))
                                   })))
                        }
                        else
                        {
                            return(FALSE)
                        }
                    }
                    else
                    {
                        return(FALSE)
                    }
                }
           })
    
    if (all(valid.list) == TRUE)
    {
        return(TRUE)
    }
    else
    {
        return(paste("Invalid input for: ", names(valid.list)[valid.list == FALSE]))
    }
}

setClass(Class="TableSchemaList", representation=list(tab.list="list", search.cols="list"), prototype=prototype(tab.list=default.tab.list(), search.cols=default.search.cols()), validity=valid.TableSchemaList)

setMethod("show", signature("TableSchemaList"), function(object)
          {
                message("An object of class TableSchemaList")
          })

setMethod("subset", signature("TableSchemaList"), function(x, table.name)
          {
            if (all(table.name %in% names(x@tab.list)) == FALSE)
            {
                stop("ERROR: Only valid names can be used for subsetting")
            }
            
            return(new("TableSchemaList", tab.list=x@tab.list[table.name]))
          })

SangerTableSchemaList <- function()
{
    new("TableSchemaList")
}

return.element <- function(use.obj, name)
{
    return(sapply(use.obj@tab.list, "[[", name))
}

setGeneric("searchTables", def=function(obj, ...) standardGeneric("searchTables"))
setMethod("searchTables", signature("TableSchemaList"), function(obj, name)
          {
            return(sapply(obj@search.cols[name], "[[", "table"))
          })

setGeneric("searchCols", def=function(obj, ...) standardGeneric("searchCols"))
setMethod("searchCols", signature("TableSchemaList"), function(obj, name)
          {
            return(sapply(obj@search.cols[name], "[[", "column"))
          })

setGeneric("searchDict", def=function(obj, ...) standardGeneric("searchDict"))
setMethod("searchDict", signature("TableSchemaList"), function(obj, name, value=NULL)
          {
            if (missing(value) || is.null(value) || is.na(value))
            {
                return(lapply(obj@search.cols[name], "[[", "dict"))
            }
            else
            {
                return(sapply(lapply(obj@search.cols[name], "[[", "dict"), "[", value))
            }
          })

#need to go through, either here or in the object validation function to ensure that the expected foreign keys match up to the observed
setGeneric("foreignExtKeyCols", def=function(obj, ...) standardGeneric("foreignExtKeyCols"))
setMethod("foreignExtKeyCols", signature("TableSchemaList"), function(obj, table.name)
          {
            as.character(sapply(return.element(obj, "foreign.keys")[table.name], function(x)
                   {
                        as.character(sapply(names(x), function(y)
                               {
                                    return(x[[y]]$ext.keys)
                               }))
                   }))
          })

setGeneric("foreignLocalKeyCols", def=function(obj, ...) standardGeneric("foreignLocalKeyCols"))
setMethod("foreignLocalKeyCols", signature("TableSchemaList"), function(obj, table.name)
          {
            as.character(sapply(return.element(obj, "foreign.keys")[table.name], function(x)
                   {
                        as.character(sapply(names(x), function(y)
                               {
                                    return(x[[y]]$local.keys)
                               }))
                   }))
          })

setGeneric("foreignExtKeySchema", def=function(obj, ...) standardGeneric("foreignExtKeySchema"))
setMethod("foreignExtKeySchema", signature("TableSchemaList"), function(obj, table.name)
          {
            as.character(sapply(return.element(obj, "foreign.keys")[table.name], function(x)
                   {
                        as.character(sapply(names(x), function(y)
                               {
                                    colSchema(obj, y, mode="normal")[match(x[[y]]$ext.keys, colNames(obj, y, mode="normal"))]
                               }))
                   }))
          })

setGeneric("bindDataFunction", def=function(obj, ...) standardGeneric("bindDataFunction"))
setMethod("bindDataFunction", signature("TableSchemaList"), function(obj, table.name, bind.vals, mode=c("normal", "merge"))
        {
            table.mode <- match.arg(mode)
            
            vcf.dta <- return.element(subset(obj, table.name), "dta.func")[[1]](bind.vals)
            
            cur.cols <- colNames(obj, table.name, mode=table.mode)
            cur.schema <- colSchema(obj, table.name, mode=table.mode)
            
            #don't need to supply columns which are autoincremented, they will be automatically added to the data.frame
            auto.col <- cur.cols[cur.schema == "INTEGER PRIMARY KEY AUTOINCREMENT"]
            #stopifnot(length(auto.col) == 1)
            
            diff.cols <- setdiff(cur.cols, colnames(vcf.dta))
            
            if (length(diff.cols) == 0)
            {
                return(vcf.dta[,cur.cols])
            }
            else if (length(diff.cols) == 1 && diff.cols == auto.col)
            {
                temp.vcf.dta <- cbind(vcf.dta, NA_integer_)
                names(temp.vcf.dta) <- c(names(vcf.dta), auto.col)
                return(temp.vcf.dta[,cur.cols])
            }
            else
            {
                nf.cols <- setdiff(cur.cols, colnames(vcf.dta))
                stop(paste("ERROR: Cannot find column(s)", paste(nf.cols, collapse=",")))
            }
            
        })

setGeneric("shouldIgnore", def=function(obj, ...) standardGeneric("shouldIgnore"))
setMethod("shouldIgnore", signature("TableSchemaList"), function(obj, table.name)
        {
            return(return.element(subset(obj, table.name), "should.ignore"))
        })

setGeneric("shouldMerge", def=function(obj, ...) standardGeneric("shouldMerge"))
setMethod("shouldMerge", signature("TableSchemaList"), function(obj, table.name=NULL)
        {
            if (missing(table.name) || is.null(table.name))
            {
                sub.obj <- obj
            }
            else
            {
                sub.obj <- subset(obj, table.name)
            }
            
            return(any(sapply(return.element(sub.obj, "foreign.keys"), is.null) == FALSE))
        })

setGeneric("schemaNames", def=function(obj, ...) standardGeneric("schemaNames"))
setMethod("schemaNames", signature("TableSchemaList"), function(obj)
        {
            names(obj@tab.list)
        })

setGeneric("tableConstr", def=function(obj, ...) standardGeneric("tableConstr"))
setMethod("tableConstr", signature("TableSchemaList"), function(obj, table.name, mode=c("normal", "merge"))
        {
            table.mode <- match.arg(mode)
            
            ret.el <- return.element(subset(obj, table.name), "db.constr")
            
            if (is.null(ret.el) || is.na(ret.el) || table.mode == "merge")
            {
                return("")
            }
            else
            {
                return(ret.el)
            }
        })

setGeneric("tableName", def=function(obj, ...) standardGeneric("tableName"))
setMethod("tableName", signature("TableSchemaList"), function(obj, table.name, mode=c("normal", "merge"))
          {
                table.mode <- match.arg(mode)
            
                if (table.name %in% names(obj@tab.list) == FALSE)
                {
                    stop("ERROR: Invalid table.name supplied")
                }
                
                return(switch(table.mode, normal=table.name, merge=paste(table.name, "temp", sep="_")))
          })

setGeneric("colNames", def=function(obj, ...) standardGeneric("colNames"))
setMethod("colNames", signature("TableSchemaList"), function(obj, table.name, mode=c("normal", "merge"))
          {
                table.mode <- match.arg(mode)
                sub.obj <- subset(obj, table.name)
                base.cols <- as.character(return.element(sub.obj, "db.cols"))
                base.schema <- as.character(return.element(sub.obj, "db.schema"))
                foreign.cols <- foreignExtKeyCols(obj, table.name)
                
                #also remove the columns that will be present in the final table but not part of the initial table
                rm.cols <- foreignLocalKeyCols(obj, table.name)
                
                return(switch(table.mode, normal=base.cols, merge=c(foreign.cols, base.cols[base.schema != "INTEGER PRIMARY KEY AUTOINCREMENT" & base.cols %in% rm.cols == FALSE])))
          })

setGeneric("colSchema", def=function(obj, ...) standardGeneric("colSchema"))
setMethod("colSchema", signature("TableSchemaList"), function(obj, table.name, mode=c("normal", "merge"))
          {
                table.mode <- match.arg(mode)
                sub.obj <- subset(obj, table.name)
                base.schema <- as.character(return.element(sub.obj, "db.schema"))
                base.cols <- as.character(return.element(sub.obj, "db.cols"))
                foreign.schema <- foreignExtKeySchema(obj, table.name)
                
                #also remove the columns that will be present in the final table but not part of the initial table
                rm.cols <- foreignLocalKeyCols(obj, table.name)
                rm.schema <- base.cols %in% rm.cols
                
                return(switch(table.mode, normal=base.schema, merge=c(foreign.schema, base.schema[base.schema != "INTEGER PRIMARY KEY AUTOINCREMENT" & rm.schema == FALSE])))
          })

setGeneric("createTable", def=function(obj, ...) standardGeneric("createTable"))
setMethod("createTable", signature("TableSchemaList"), function(obj, table.name, mode=c("normal", "merge"))
          {
                table.mode <- match.arg(mode)
                if (shouldMerge(obj, table.name) == TRUE && table.mode == "merge")
                {
                    temp.str <- "TEMPORARY"
                }
                else if (table.mode == "normal")
                {
                    temp.str <- ""
                }
                else
                {
                    stop("ERROR: Cannot generate statement with mode set to 'merge' and a NULL foreign.key element")
                }
                
                use.cols <- colNames(obj, table.name, mode=table.mode)
                use.schema <- colSchema(obj, table.name, mode=table.mode)
                
                tab.constr <- tableConstr(obj, table.name, mode=table.mode)
                tab.constr <- ifelse(tab.constr == "", tab.constr, paste0(",", tab.constr))
                
                return(paste("CREATE",temp.str,"TABLE", tableName(obj, table.name, table.mode), "(", paste(paste(use.cols, use.schema), collapse=","), tab.constr, ")"))
          })

setGeneric("mergeStatement", def=function(obj, ...) standardGeneric("mergeStatement"))
setMethod("mergeStatement", signature("TableSchemaList"), function(obj, table.name)
          {
                #currently, probably the temporary table
                cur.db <- tableName(obj, table.name, mode="merge")
                #table trying to create
                target.db <- tableName(obj, table.name, mode="normal")
                
                target.cols <- colNames(obj, table.name, mode="normal")
                target.schema <- colSchema(obj, table.name, mode="normal")
                #remove the autoincrement column first
                target.cols <- target.cols[target.schema != "INTEGER PRIMARY KEY AUTOINCREMENT"]
                paste.targs <- paste(target.cols, collapse=",")
                
                #create the join statement using the foreign.keys slot
                
                fk <- return.element(obj, "foreign.keys")[[table.name]]
                
                if (is.null(fk))
                {
                    stop("ERROR: Cannot generate statement if the foreign key element is NULL")
                }
                
                keys <- sapply(names(fk), function(y)
                       {
                            return(paste(fk[[y]]$ext.keys, collapse=","))
                       })
                        
                join.statement <- paste(paste("JOIN", names(keys), "USING", paste0("(", keys,")")), collapse=" ")
                
                if (shouldIgnore(obj, table.name))
                {
                    ignore.str <- "OR IGNORE"
                }
                else
                {
                    ignore.str <- ""
                }
                
                return(paste("INSERT",ignore.str,"INTO", target.db, "(", paste.targs,") SELECT", paste.targs,"FROM", cur.db , join.statement))
          })

setGeneric("insertStatement", def=function(obj, ...) standardGeneric("insertStatement"))
setMethod("insertStatement", signature("TableSchemaList"), function(obj, table.name, mode=c("normal", "merge"))
          {
                table.mode <- match.arg(mode)
                
                use.cols <- colNames(obj, table.name, mode=table.mode)
                
                if (shouldIgnore(obj, table.name))
                {
                    ignore.str <- "OR IGNORE"
                }
                else
                {
                    ignore.str <- ""
                }
                
                if (any(use.cols == "character(0)"))
                {
                    stop("ERROR: Cannot run statement with mode 'merge' and a NULL foreign.key element")
                }
                
                return(paste("INSERT",ignore.str,"INTO", tableName(obj, table.name, table.mode), "VALUES (", paste(paste0(":", use.cols), collapse=","), ")"))
          })

