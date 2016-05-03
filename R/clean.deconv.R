#' Remove end adaptor
#' 
#' This function is used to remove end adaptors, starting from a fastq file.
#' 
#' As mentioned in the general description of this package, most functions are 
#' tailored to Illumina architecture. \code{rmEndAdaptor} was developed to 
#' remove the P7 adaptor at the end of single-reads. However, it can be actually
#' used to remove any 'tail'. Reads are trimmed at the first position of the 
#' match with the passed. \code{EndAdaptor} pattern.  The sequence of the 
#' \code{EndAdaptor} is passed (as character vector) in a 5' to 3' direction and
#' it is internally reversed and complemented. Other than the sequence, it is 
#' possible to pass the character vector "P7" or "P7_last10". With the first, 
#' the sequence of the P7 adaptor is selected (CAAGCAGAAGACGGCATACGAGAT). With
#' the latter a partial match is searched for (the last 10 bp: CATACGAGAT).
#' 
#' Matches are searched using \code{\link[Biostrings]{vmatchPattern}}, with 
#' \code{adaptor.mismatch} used for max.mismatch (min.mismatch is fixed to 
#' zero).
#' 
#' The search is conducted with \code{fixed=TRUE}, which means (from Biostring):
#' "an IUPAC ambiguity code in the pattern can only match the same code in the 
#' subject, and vice versa".
#' 
#' @param fn Fully qualified name (i.e. the complete path) of the fastq file
#' @param nRead The number of bytes or characters to be read at one time. See 
#'   \code{\link[shortRead]{FastqStreamer}} for details
#' @param EndAdaptor A character vector with the sequence of the end adaptor. 
#'   See details.
#' @param adaptor.mismatch The maximum number of allowed mismatch. See details.
#' @export
#' @return A fasta file with the reads where the end adaptor was found (and 
#'   removed) saved in the same location where the input data was located. The
#'   file is named with the suffix "_EndAdRm". A list with the
#'   total number of reads that were processed and retained is also returned.
#'   
rmEndAdaptor <- function(fn, nRead=1e8, EndAdaptor="P7_last10", adaptor.mismatch=0) {
  library(ShortRead)
  #----------------------------------------------------------------------------#
  # Helper functions
  #----------------------------------------------------------------------------#
  
  putzero <- function(x){
    if(is.null(x)) x <- 0 else x <- x - 1
    return(x)
  }
  #----------------------------------------------------------------------------#
  
  if(EndAdaptor == "P7") EndAdaptor <- "CAAGCAGAAGACGGCATACGAGAT"
  if(EndAdaptor == "P7_last10") EndAdaptor <- "CATACGAGAT"  
  
  stream <- FastqStreamer(fn, nRead)
  on.exit(close(stream))
  nRetained <- 0
  
  while (length(fq <- yield(stream))) {
    seqs <- sread(fq)
    qual <- quality(fq)
    qual <- quality(qual)
    P7_hits <- vmatchPattern(pattern=reverseComplement(DNAString(EndAdaptor)), 
                             subject=seqs,
                             max.mismatch=adaptor.mismatch, 
                             min.mismatch=0,
                             with.indels=FALSE, fixed=TRUE,
                             algorithm="auto")
    
    if(sum(elementLengths(P7_hits) > 1) > 0) {
      message(cat(paste("Warning: some reads in", fn, "have more than one match with:", 
                        EndAdaptor, "The first match from the left is used",
                        sep="\n")))
    }
    ends <- startIndex(P7_hits)
    ends <- lapply(ends, "[", 1)
    ends <- lapply(ends, putzero)
    
    seqs <- DNAStringSet(seqs, start=1, end=unlist(ends))
    qual <- BStringSet(qual, start=1, end=unlist(ends))
    qual <- SFastqQuality(qual) 
    trimmed <- ShortReadQ(sread=seqs, quality=qual, id=id(fq))
    trimmed <- trimmed[width(trimmed) > 0]
    
    fname <- paste0(substr(fn, start=1, stop=regexpr(".fastq", fn)[1] - 1), 
                    "_EndAdRm.fastq.gz")
    if(file.exists(fname)) message(cat(paste("Warning: the file", fname, 
                                             "already exists\nSequences were appended")))
    writeFastq(trimmed, fname, mode="a")
    nRetained <- nRetained + length(trimmed)
  }
  nInitial <- stream$status()["total"]
  message(paste("Processed", nInitial, "- retained", nRetained))
  return(rmEndAd <- list(nInitialReads=nInitial, nRetained=nRetained))
}
