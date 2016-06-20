#' Create a phyloseq object
#' 
#' \code{make.phyloseq} creates a \code{\link[phyloseq]{phyloseq}} object starting
#'   from \code{data.proc} output. Sample metadata and taxonomic assignments can
#'   be optionally included. If a minimum of two elements are passed to 
#'   \code{make.phyloseq}, then a phyloseq objext is returned, otherwise, 
#'   \code{data.proc} output will be returned as 
#'   \code{\link[phyloseq]{otu_table-class}}.
#'    
#'   
#' A \code{data.frame} can be passed with \code{sample.table} to include sample 
#'   metadata. In this case, the row names have to be identical to the sample 
#'   names included in the  \code{data.proc} output. The safest way to make sure
#'   this is the case is to use:
#'   \code{row.names(sampleTable) <- row.names(dproc$stable)} where 
#'   \code{sampleTable} is the  \code{data.frame} with the metadata. Note that 
#'   this assumes that the samples are in the same order as in the the sequence
#'   table.
#' @param dproc The \code{\link{data.proc}} output
#' @param sample.table The data.frame with sample metadata 
#' @param tax.table A matrix with taxonomic assignments (see
#'   \code{\link[phyloseq]{phyloseq}} for details)
#' @param phy.tree A phylogenetic tree (see
#'   \code{\link[phyloseq]{phyloseq}} for details)
#' @seealso \code{\link{data.proc}}, \code{\link[phyloseq]{phyloseq}}
#' @export
#' @importFrom phyloseq phyloseq
make.phyloseq <- function(dproc, sample.table=NULL, tax.table=NULL, 
                          phy.tree=NULL ) {
  suppressWarnings(library(phyloseq, quietly=TRUE, warn.conflicts=FALSE))
  stable <- dproc$stable
  colnames(stable) <- dproc$seq_list$sequence
  ps <- phyloseq(otu_table(stable, taxa_are_rows=FALSE), 
                 if(!is.null(sample.table)) sample_data(sample.table), 
                 if(!is.null(tax.table)) tax_table(tax.table),
                 if(!is.null(phy.tree)) phy_tree(phy.tree))
  
  return(ps)
}

