
#' Query UCSC for Recombination data 
#' 
#' Adds recombination data to a 'locus' object by querying UCSC genome browser.
#' 
#' @param loc Object of class 'locus' generated by [locus()]
#' @param genome Either `"hg38"` or `"hg19"`
#' @param table Optional character value specifying which recombination table to
#'   use.
#' @details
#' Uses the `rtracklayer` package to query UCSC genome browser for recombination
#' rate data. The results are cached using `memoise` to reduce API requests.
#' 
#' Possible options for `table` for hg19 are `"hapMapRelease24YRIRecombMap"`,
#' `"hapMapRelease24CEURecombMap"`, `"hapMapRelease24CombinedRecombMap"` (the
#' default).
#' 
#' The only option for `table` for hg38 is `"recomb1000GAvg"` (the default).
#' 
#' @returns A list object of class 'locus'. Recombination data is added as list
#'   element `recomb`.
#' @importFrom GenomeInfoDb genome<-
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom rtracklayer browserSession ucscTableQuery getTable
#' @export
#'
link_recomb <- function(loc, genome = "hg38", table = NULL) {
  if (!inherits(loc, "locus")) stop("Not a locus object")
  loc$recomb <- mem_query_recomb(genome, loc$xrange, loc$seqname, table)
  loc
}


query_recomb <- function(gen, xrange, seqname, table = NULL) {
  if (is.null(table)) {
    table <- if (gen == "hg38") {"recomb1000GAvg"
    } else if (gen == "hg19") "hapMapRelease24CombinedRecombMap"
  }
  if (!grepl("chr", seqname)) seqname <- paste0("chr", seqname)
  gr <- GRanges(ranges = IRanges(start = xrange[1], end = xrange[2]),
                seqnames = seqname)
  message("Retrieving recombination data from UCSC")
  session <- browserSession("UCSC")
  genome(session) <- gen
  query <- ucscTableQuery(session, table = table, range = gr)
  getTable(query)
}

# use memoise to reduce calls to UCSC API
mem_query_recomb <- memoise(query_recomb)