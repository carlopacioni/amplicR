#' Detect a reference sequence
#' 
#' \code{detect} takes in either the output of \code{data.proc}, or a list of 
#' sequences from fasta files, and compare them with a reference sequence 
#' reporting the number of mismatch using \code{srdistance} from the package 
#' \code{shortRead}.
#' 
#' The output from \code{data.proc} can be passed with \code{data}. If no data 
#' is passed to \code{detect}, then it will read fasta files in the directory 
#' passed with \code{dir.in}. If \code{dir.in=NULL}, then an interactive window 
#' will open to select the location of the files. As for \code{data.proc}, 
#' \code{detect} assumes that each file represents a sample.
#' 
#' If both \code{dir.out=NULL} and \code{dir.in=NULL}, then the path where to 
#' save the results will be asked with an interactive window, otherwise  
#' \code{dir.out <- dir.in}
#' 
#' A summary of the number of sequences found and the minimum number of mismatch
#' within each sample is returned as \code{data.frame}. Each comlumn reporting 
#' the number of mismatch is named with the named of the character vector passed
#' with \code{ref_seqs}. These results are also written to disk (Summary.csv) 
#' together with the alignments of the sequences provided with the reference 
#' sequence (in the folder "Final_alns"). The alignment is built using 
#' \code{PairwiseAlignments} from the package \code{Biostrings}.
#' 
#' @param data The output from \code{data.proc}
#' @param ext A character vector (of length=1) with the extension of the fasta 
#'   files. (Default "fasta")
#' @param ref_seqs A named character vector with the reference sequence(s)
#' @inheritParams data.proc
#' @return A \code{data.frame} with the number of sequences found in each samle
#'   and the minimum mismatch count for each reference sequence. These results
#'   are also also written to file (see Details)
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
#' # Naming the reference sequence
#' names(HTJ) <- "HTJ" 
#' 
#' # Use 'Detect' to verify the presence of Mycobacteriumavium subspecies 
#' # paratuberculosis
#' 
#' det <- detect(HTJ.test, dir.out=out, ref_seqs=HTJ)
#' # To clean up the temp directory
#' unlink(out, recursive=TRUE)

detect <- function(data=NULL, dir.in=NULL, dir.out=NULL, ext="fasta", ref_seqs) {
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
    dt <- data.table(Sample, Seq, nDiff)
    setnames(dt, "nDiff", names(ref_seq))
    return(dt)
  }
  
  make.aln <- function(DNAstr, ref_seq) {
    aln <- pairwiseAlignment(DNAstr, ref_seq)
    return(aln)
  }
  
  w.aln <- function(i, lalns, sample_names, dir.out, aln_fold) {
    dir.create(paste(dir.out, aln_fold, sep="/"), 
               showWarnings=FALSE, 
               recursive=TRUE)
    writePairwiseAlignments(lalns[[i]], paste(dir.out, aln_fold, 
                                            paste0(sample_names[[i]], ".fasta"), 
                                                  sep="/"), 
                                            block.width=250)
  }
    
    #----------------------------------------------------------------------------#
  if(is.null(data)) {
    if(is.null(dir.in)) {
      dir.in <- choose.dir(caption="Please, select the directory where the fasta
                         files are located")
    }
  }
  
  if(is.null(data)) {
    fns <- list.files(path=dir.in, full.names=TRUE)
    fasta_files <- fns[grepl(paste0(".", ext, "$"), fns)]
    lDNAstr <- lapply(fasta_files, readFasta)
    sample_names <- sub(paste0(".", ext, "$"), "", basename(fasta_files))
    names(lDNAstr) <- sample_names
    nSeq <- sapply(lDNAstr, length)
    lsummary <- list(data.frame(Sample=names(lDNAstr), nSeq))
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
      
    #### Detect ####
  el <- length(lsummary)
  lres <- list()
  nref_seqs <- seq_along(ref_seqs)
  aln.out <- paste(dir.out, "Final_alns", sep="/")
  
  for (r in nref_seqs) {
    lm <- lapply(seq_along(lDNAstr), diseased, lDNAstr, names(lDNAstr), 
                 ref_seq=ref_seqs[r])
    lres[[r]] <- rbindlist(lm)
    
    alns <- lapply(lDNAstr, make.aln, ref_seqs[r])
    lapply(seq_along(alns), w.aln, alns, names(lDNAstr), aln.out, names(ref_seqs[r]))
  }
  
  results <- plyr::join_all(lres, by=c('Sample', "Seq"), type="left")
  write.csv(results, 
            file=paste(dir.out, "Results.csv", sep="/"), 
            row.names=FALSE)
  
  el <-el + 1
  results <- data.table(results)
  lsummary[[el]] <- results[, lapply(.SD, min), by="Sample", .SDcols=names(ref_seqs)]
  
  summary <- plyr::join_all(lsummary, by='Sample', type="left")
  write.csv(summary, file=paste(dir.out, "Summary.csv", sep="/"), 
            row.names=FALSE)
  
  
  return(summary)
  }
