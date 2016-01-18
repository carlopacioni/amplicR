#' Detect a reference sequence
#' 
#' \code{detect} takes in either the output of \code{data.proc}, or a list of 
#' sequences from fasta files, and compare them with a reference sequence 
#' reporting the number of mismatch using \code{srdistance} from the package 
#' \code{shortRead}.
#' 
#' The output from \code{data.proc} can be passe with \code{data}. If no data is
#' passed to \code{detect}, then it will read fasta files in the directory 
#' passed with \code{dir.in}. If \code{dir.in=NULL}, then an interactive window 
#' will open to select the location of the files. As for \code{data.proc}, 
#' \code{detect} assumes that each file represents a sample.
#' 
#' A summary of the number of sequences found and the number of mismatch is 
#' returned as \code{data.frame} as well as being written to disk (Summary.csv) 
#' together with the alignments of the sequences provided with the reference 
#' sequence (in the folder "Final_alns"). The alignment is built using 
#' \code{PairwiseAlignments} from the package \code{Biostrings}.
#' 
#' @param data The output from \code{data.proc}
#' @param ext A character vector (of length=1) with the extension of the fasta
#'   files. (Default "fasta")
#' @param ref_seq A character vector with the reference sequence  
#' @inheritParams data.proc
#' @return A list where each element is an alignment of the sequences within 
#'   each sample with the reference sequence. The alignments and the number of 
#'   sequences found and mismatch is also written to file (see Details)
#' @import data.table
#' @export
#' @examples 
#' # Select the directory where the example data are stored
#' example.data <- system.file("extdata", "HTJ", package="amplicR")
#' # Select a temporary directory where to store the outputs
#' out <- tempdir()
#' # Process raw data
#' HTJ.test <- data.proc(example.data, out, bp=140)
#' # Referece Mycobacteriumavium subspecies paratuberculosis sequence 
#' HTJ <- "CTGCGCGCCGGCGATGACATCGCAGTCGAGCTGCGCATCCTGACCAGCCGACGTTCCGATCTGGTGGCTGATCGGACCCGGGCGATCGAACCGAATGCGCGCCCAGCTGCTGGAATACTTTCGGCGCTGGAACGCGCCTT"
#' 
#' 
#' # Use 'Detect' to verify the presence of Mycobacteriumavium subspecies 
#' # paratuberculosis
#' 
#' det <- detect(HTJ.test, dir.out=out, ref_seq=HTJ)
#' # To clean up the temp directory
#' unlink(out, recursive=TRUE)

detect <- function(data=NULL, dir.in=NULL, dir.out=NULL, ext="fasta", ref_seq) {
  #----------------------------------------------------------------------------#
  library(dada2)
  library(ShortRead)
  #----------------------------------------------------------------------------#
  # Helper functions
  #----------------------------------------------------------------------------#
  
  diseased <- function(i, lDNAstr, s_names, ref_seq) {
    nDiff <- unlist(srdistance(lDNAstr[[i]], ref_seq))
    Seq <- 1:length(nDiff)
    Sample <- rep(s_names[[i]], length(nDiff))
    df <- data.frame(Sample, Seq, nDiff)
    return(df)
  }
  
  make.aln <- function(DNAstr, ref_seq) {
    aln <- pairwiseAlignment(DNAstr, ref_seq)
    return(aln)
  }
  
  w.aln <- function(i, lalns, sample_names, dir.out, aln_fold) {
    writePairwiseAlignments(lalns[[i]], paste(dir.out, aln_fold, 
                                            paste0(sample_names[[i]], ".fasta"), 
                                                  sep="/"), 
                                            block.width=250)
  }
    
    #----------------------------------------------------------------------------#
  ref_seq <- DNAString(ref_seq)
  
  if(is.null(data)) {
    if(is.null(dir.in)) {
      dir.in <- choose.dir(caption="Please, select the directory where the fasta
                         files are located")
    } else {
      fns <- list.files(path=dir.in, full.names=TRUE)
      fasta_files <- fns[grepl(paste0(".", ext, "$"), fns)]
      lDNAstr <- lapply(fasta_files, readFasta)
      sample_names <- sub(paste0(".", ext, "$"), "", basename(fasta_files))
      names(lDNAstr) <- sample_names
      nSeq <- sapply(lDNAstr, length)
      lsummary <- list(data.frame(Sample=names(lDNAstr), nSeq))
      }
  } else {
    lDNAstr <- data[[1]]
    lsummary <- data[[2]]
    }
  
  if(is.null(dir.out)) {
    if(is.null(dir.in)) {
      dir.out <- choose.dir(caption="Please, select the directory where to save 
                            results") 
    } else {
      dir.out <- dir.in
    }
  }
      
  dir.create(dir.out, showWarnings=FALSE, recursive=TRUE)
  
    #### Detect ####
    el <- length(lsummary)
    lm <- lapply(seq_along(lDNAstr), diseased, lDNAstr, names(lDNAstr), ref_seq)
    res <- rbindlist(lm)
    write.csv(res, file=paste(dir.out, "Results.csv", sep="/"), row.names=FALSE)
    el <-el + 1
    lsummary[[el]] <- res[, .(nDiff=min(nDiff)), by="Sample"]
    
    summary <- plyr::join_all(lsummary, by='Sample', type="left")
    write.csv(summary, file=paste(dir.out, "Summary.csv", sep="/"), 
              row.names=FALSE)
    
    alns <- lapply(lDNAstr, make.aln, ref_seq)
    aln_fold <- "Final_alns"
    dir.create(paste(dir.out, aln_fold, sep="/"), showWarnings=FALSE, 
               recursive=TRUE)
    
    lapply(seq_along(alns), w.aln, alns, names(lDNAstr), dir.out, aln_fold)
    return(summary)
  }
