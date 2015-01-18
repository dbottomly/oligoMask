#TODO: if no strains/sample names are specified, parse based on the REF and ALT and maybe another field that would indicate populations

#old definition
#setClass(Class="VcfDB", representation=list(db.path="character", tbsl="TableSchemaList", start.var.table="character", end.var.table="character", var.mask.probe.id="character", var.mask.var.id="character"),
#         prototype=prototype(start.var.table="reference", end.var.table="probe_info", var.mask.probe.id="probe_id", var.mask.var.id="ref_id"))

sanger.search.cols <- function()
{
    return(list(mapping.status=list(table="probe_info", column="align_status", dict=c(unique="UniqueMapped", multi="MultiMapped", non="UnMapped")),
                genotype.filter=list(table="reference", column="filter", dict=c(`TRUE`="TRUE", `FALSE`="FALSE"))))
}

setClass(Class="VcfDB", contains="Database", representation=list(start.var.table="character", end.var.table="character", var.mask.probe.id="character",
                                                                 var.mask.var.id="character",search.cols="list"),
         prototype=prototype(var.mask.probe.id="probe_info.probe_id", var.mask.var.id="reference.ref_id", search.cols=sanger.search.cols(), tbsl=SangerTableSchemaList()))


setGeneric("searchTables", def=function(obj, ...) standardGeneric("searchTables"))
setMethod("searchTables", signature("VcfDB"), function(obj, name)
          {
            return(sapply(obj@search.cols[name], "[[", "table"))
          })

setGeneric("searchCols", def=function(obj, ...) standardGeneric("searchCols"))
setMethod("searchCols", signature("VcfDB"), function(obj, name, include.table=F)
          {
	    if (include.table == T)
	    {
		return(sapply(obj@search.cols[name], function(x)  paste(x[["table"]], x[["column"]], sep=".")))
	    }else{
		return(sapply(obj@search.cols[name], "[[", "column"))
	    }
          })

setGeneric("searchDict", def=function(obj, ...) standardGeneric("searchDict"))
setMethod("searchDict", signature("VcfDB"), function(obj, name, value=NULL)
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

###replace with dbFile
#setGeneric("vcfDb", def=function(obj, ...) standardGeneric("vcfDb"))
#setMethod("vcfDb", signature("VcfDB"), function(obj)
#          {
#                return(obj@db.file)
#          })

setGeneric("getProbeVars", def=function(obj, ...) standardGeneric("getProbeVars"))
setMethod("getProbeVars", signature("VcfDB"), function(obj, probe.ids)
          {
			if (is.null(probe.ids) || is.na(probe.ids) || is.character(probe.ids) == FALSE || length(probe.ids) == 0)
			{
				stop("ERROR: probe.ids needs to be a non-empty character vector")
			}
			
			temp.probes <- select(test.db, probe_id, fasta_name, align_status, probe_chr, probe_start, probe_end, seqnames, start,
						end, filter, geno_chr, genotype.allele_num, strain)
			
			filt.probes <- filter(temp.probes, probe_id %in% probe.id)
			
			return(filt.probes)
          })

###replace with schema...
#setGeneric("tbsl", def=function(obj, ...) standardGeneric("tbsl"))
#setMethod("tbsl", signature("VcfDB"), function(obj)
#          {
#                return(obj@tbsl)
#          })

setMethod("show", signature("VcfDB"), function(object)
        {
			#a hack for now to determine whether the object is a package
			split.path <- strsplit(dbFile(object), .Platform$file.sep)[[1]]
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
			
			#db.con <- dbConnect(SQLite(), dbFile(object))
			
			if (isTRUE(all.equal(schema(object), SangerTableSchemaList())))
			{
				print('in progress...')
				#probe.count <- summarize(group_by(select(object, probe_info.probe_id), probe_id), n_distinct(probe_id))
				#probe.count <- dbGetQuery(db.con, "SELECT COUNT(DISTINCT(probe_id)) FROM probe_info")[,1]
				#message(paste("Containing alignments for", probe.count, "probes"))
				#
				#var.count <- dbGetQuery(db.con, "SELECT COUNT(DISTINCT(ref_id)) FROM reference")[,1]
				#strain.count <- dbGetQuery(db.con, "SELECT COUNT(DISTINCT(strain)) FROM genotype")[,1]
				#message(paste("Containing", var.count , "variants from", strain.count, "inbred strains"))
			}
			else
			{
				stop("ERROR: Unknown type of TableSchemaList")
			}
			
			#invisible(dbDisconnect(db.con))
        })

#summary method
#probe.cat.counts <- dbGetQuery(db.con, "SELECT align_status AS Alignment_Status, COUNT(DISTINCT(probe_ind)) AS Count from probe_info GROUP BY align_status")
#				format(probe.cat.counts)
