#' Write unique sequences to a fasta file for each sample where sequence names  
#' match the sequence table and include abundance data
#' 
#' \code{table2fasta} takes the output from 
#' \code{\link[dada2]{makeSequenceTable}} or the table element from 
#' \code{data.proc} as input and generates a fasta file for each sample with 
#' unique sequences. It uses internally \code{\link[dada2]{uniquesToFasta}}. 
#' However, sequences are named with a sequential number, which is consistent 
#' across all samples and abundance data is appended to the sequence names.
#' 
#' The input can be either the \code{\link[dada2]{makeSequenceTable}} or 
#' \code{data.proc} output. Note that \code{data.proc} output is actually a 
#' modified \code{\link[dada2]{makeSequenceTable}} output where the column 
#' names, which normally are the actual sequeces, have been replaced for 
#' convinience. \code{data.proc} reports the sequence IDs and matching sequences 
#' in a separated elelment of the returned list. If the input is 
#' from \code{data.proc}, then also the \code{seq_list} element needs to be 
#' passed to \code{table2fasta}.
#' 
#' In each file, sequences are named with a sequential number, which is 
#' consistent across all samples. For example, if sample1 has sequence seq1 and 
#' seq2 and sample2 has seq1 and seq3, the sequences in the fasta file for 
#' sample1 will be named seq1 and seq2, and those in the fasta file for sample2 
#' will be named seq1 and seq3, so that it will be possible to match back the 
#' same sequences in different samples in down stream analysis.
#' 
#' 
#' 
#' @param stable The sequence table (output of \code{dada2::makeSequenceTable})
#'   or \code{data.proc}
#' @param seq.list The seq_list element from \code{data.proc}. Required if the 
#'   input table is the output from \code{data.proc}. If NULL (default), it is 
#'   assumed that the input table is from \code{dada2::makeSequenceTable} and 
#'   sequences are named with a seq# pattern, where # is a sequential number 
#'   respect to the sequence in the table.
#' @param dir.out The directory where to save fasta files
#' @inheritParams data.proc
#' @seealso \code{\link{data.proc}},
#'   \code{\link[dada2]{uniquesToFasta}}, \code{\link[dada2]{makeSequenceTable}}
#' @return A fasta file for each sample
#' @export
#' 
table2fasta <- function(stable, seq.list=NULL, dir.out=NULL, verbose=TRUE) {
  #----------------------------------------------------------------------------#
  if (!requireNamespace("dada2", quietly = TRUE)) {
    stop("Package 'dada2' needed for this function to work. Please install it 
         either manually or using the function amplicR::setup().",
         call. = FALSE)
  }  
  #----------------------------------------------------------------------------#
  # Helper functions
  #----------------------------------------------------------------------------#
  w.fasta <- function(i, stable, seqs, seq_names, verbose) {
    uniq_vect <- as.integer(stable[i,])
    retain <- as.logical(unname(uniq_vect))
    uniq_vect <- uniq_vect[retain]
    names(uniq_vect) <- seqs[retain]
    seq_names_sub <- paste0(seq_names[retain], ";size=", unname(uniq_vect), ";")
    if(verbose) 
      message(paste("Writing fasta file", paste0(row.names(stable)[i], ".fasta")))
    dada2::uniquesToFasta(uniq_vect, 
                   paste(dir.out, paste0(row.names(stable)[i], ".fasta"), sep="/"),
                   ids=seq_names_sub)
  }
  
  #----------------------------------------------------------------------------#
  if(is.null(dir.out)) choose.dir(caption="Selec folder where to save fasta files")
  dir.create(dir.out, showWarnings=FALSE, recursive=TRUE)
  
  if(is.null(seq.list)) {
    seqs <- colnames(stable)
    seq_names <- paste0("seq", 1:dim(stable)[2])
  } else {
    seqs <- seq.list$sequence
    seq_names <- seq.list$seq_names
  }
  lapply(1:dim(stable)[1], w.fasta, stable, seqs, seq_names, verbose)
}

































