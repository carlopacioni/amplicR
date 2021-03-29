#' Extracts rows from a matrix into a vector
#' 
#' The row rn becomes an element of a vector. This function is used internally by amplicR
#' @param m A matrix
#' @param rn The row number to extract
#' @return A vector of same type as the matrix
#' @export
#' 
getRow <- function(m, rn) {
  r <- as.integer(round(m[rn,], 0))
  qual <- Biostrings::quality(NumericQuality(r))
  qual <- IlluminaQuality(qual)
  return(qual)
}

#' Write a nexus file from a DNAStringSet object
#'
#' Takes a \code{character} vector or a \code{DNAStringSet} and write a nexus
#' file to disk.
#'
#' If \code{x} is a \code{character} vector, it will convert it to a
#' \code{DNAStringSet}. if \code{aln=FALSE}, \code{write.nexus} will align the
#' sequences first and then write the nexus file. If \code{aln=TRUE}, it will
#' assume that all sequences have the same length and are aligned.
#'
#' If \code{charset=TRUE}, \code{write.nexus} will append a character block at
#' the end of the file. This is sometimes used to partition the alignment when
#' it is imported in software like BEAST2
#'
#' @param x Either a \code{character} vector or a \code{DNAStringSet}
#' @param aln Logocal. Whether x is an alignment
#' @param gapOening Integer (negative). Penalty for opening a gap in the
#'   alignment
#' @param gapExtension Integer (negative). Penalty for extending a gap in the
#'   alignment
#' @param dir.out The path where to save the output
#' @param fn The output file name
#' @param charset Whether a charset block should be added at the end of the
#'   nexus file
#' @param locusIDs \code{character} vector with the names of the loci (if
#'   \code{charset=TRUE}) in the same order as they appear in the sequences
#' @param locusLength Integervector with the length of each locus (if
#'   \code{charset=TRUE}) in the same order as they appear in the sequences
#' @return Writes a nexus file to disk
#' @export
#' @examples seqs <- c("AAATTTT", "GAATTCT")
#'         names(seqs) <- c("seq1", "seq2")
#'         x <- Biostrings::DNAStringSet(seqs)
#'         locusNames <- c("locus1", "locus2")
#'         locusLen <- c(3, 4)
#'         tmpDir <- tempdir()
#'         tmpFile <- tempfile(tmpdir = tmpDir)
#'         write.nexus(x, aln=TRUE, dir.out=tmpDir, fn=tmpFile, charset=TRUE, 
#'                 locusIDs=locusNames, locusLength=locusLen)

write.nexus <- function(x, aln=FALSE, gapOpening=c(-18, -16), gapExtension=c(-2, -1), 
                        dir.out, fn, charset=FALSE, locusIDs, locusLength, verbose=FALSE) {
  whatsx <- class(x)
  if(whatsx == "character") x <- Biostrings::DNAStringSet(x)
  if(whatsx == "DNAStringSet" & aln == FALSE) {
    x <- DECIPHER::AlignSeqs(x, gapOpening=gapOpening, gapExtension=gapExtension) 
  } else {
    if(whatsx == "DNAStringSet" & aln == TRUE & verbose == TRUE) 
      message("'aln' is set to TRUE, assuming that 'x' is an alignment")
  }
  
  nex <- c("#NEXUS", 
           "begin taxa;",
           "dimensions ntax=", length(x), ";",
           "taxlabels", names(x), ";",
           "End;",
           "begin characters;",
           "dimensions nchar=", IRanges::width(x)[1], ";",
           "Format datatype=dna missing=? gap=- matchchar=.;",
           "Matrix",
           paste(names(x), x),
           ";",
           "End;")
 
  charsets <- NULL
  if(charset == TRUE) {
    pos <- 0
    charsets <- character(length(locusIDs))
    for(i in seq_along(charsets)) {
      charsets[i] <- paste("charset", locusIDs[i], "=", 1 + pos, "-", 
                           pos + locusLength[i], ";")
      pos <- pos + locusLength[i]
      charblock <- c("begin assumptions;", charsets, "end;")
    }
    
    
    
  }
  writeLines(c(nex, if(charset == TRUE) charblock), con=file.path(dir.out, fn))
}

#' Extract sequences from a derep object
#' 
#' This function is used internally to batch-extract sequences from a derep object
#' 
#' @param derepObj The derep-element (that is one element from the list that constitutes the dereb list)
#' @return A character vector with the sequences
#' @export
#' 
getSeqFromDerep <- function(derepObj) {
  reads <- list(sequence=names(derepObj$uniques))
  return(reads)
}
