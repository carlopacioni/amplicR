#' Remove end adaptor
#' 
#' This function is used to remove end adaptors, starting from a fastq file.
#' 
#' As mentioned in the general description of this package, most functions are 
#' tailored to Illumina architecture. \code{rmEndAdaptor} was developed to 
#' remove the P7 adaptor at the end of single-reads. However, it can be actually
#' used to remove any 'tail'. Reads are trimmed at the first position of the 
#' match with the passed \code{EndAdaptor} pattern.  The sequence of the 
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
#' @param EndAdaptor A character vector with the sequence of the end adaptor,
#'   "P7" or "P7_last10" (See details)
#' @param adaptor.mismatch The maximum number of allowed mismatch (See details)
#' @export
#' @return A fastq file with the reads where the end adaptor was found (and 
#'   removed) saved in the same location where the input data was located. The 
#'   file is named with the suffix "_EndAdRm". A list with the total number of
#'   reads that were processed and retained is also returned.
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
      message(cat(paste("Note: some reads in", fn, "have more than one match with:", 
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
    if(file.exists(fname)) message(cat(paste("Note: the file", fname, 
                                             "already exists\nSequences were appended")))
    writeFastq(trimmed, fname, mode="a")
    nRetained <- nRetained + length(trimmed)
  }
  nInitial <- stream$status()["total"]
  message(paste("Processed", nInitial, "- retained", nRetained))
  return(rmEndAd <- list(nInitialReads=nInitial, nRetained=nRetained))
}

#' Separate reads by genes and deconvolute them based on barcodes
#' 
#' \code{deconv} takes a fastq file and will search for the forward primer and 
#' use this to separate the reads. That is, amplicons from multiple genes will 
#' be separated. Within each gene, reads are then separated based on forward 
#' and/or reverse index.  The end products are several fastq files - one for 
#' each samples, in as many folders as genes - where primers and indexes were 
#' removed.
#' 
#' This function applies only to reads with in-line indexes. That is, where the 
#' architecture of the reads is as follows:
#' 
#' F_index---F_primer---Target_sequence---R_primer---R_index
#' 
#' Note that the P7 adapter can be removed with \code{\link{rmEndAdaptor}}.
#' 
#' It is possible to control the number of mismatch, and IUPAC ambiguities codes
#' can be used in primers (i.e.the search for the primers is conducted with 
#' \code{fixed=FALSE}, which means (from Biostring): "an IUPAC ambiguity code in
#' the pattern can match any letter in the subject that is associated with the 
#' code, and vice versa". Note that indexes are searched with \code{fixed=TRUE}.
#' 
#' Information about the reads are passed with a comma separated file (CSV), 
#' whose path and name is passed with \code{info.file}. This must contain a 
#' column for each: the foward index, the reverse index, the forward primer, the
#' reverse primer, the sample IDs, and the gene (or an identifier of gene) being
#' amplified. The column headings are passed with the function arguments. While 
#' it is possible to include other columns where the users can record additional
#' information, these are effectively ignored. It is not mandatory to have both,
#' the forward and reverse indexes, but if one is not used, there is still the 
#' need to include a blank column in  \code{info.file} and indicate the column 
#' heading.
#' 
#' All sequences for indexes and primers are passed (as character vector) in 5' 
#' to 3' direction and are internally reversed and complemented when necessary.
#' 
#' \code{deconv} initially searches for the forward primer and separates the 
#' reads creating (if not existing already) a folder named as for the relevant 
#' information provided in the column \code{gene}, which is typically an 
#' identifier of the targeted genes. It is possible to use the column 
#' \code{gene} to group genes in other logical way than genes, but all identical
#' forward primer should have the same \code{gene} information. This is because 
#' (for efficiency) \code{deconv} uses only the first line for each forward 
#' primer to identify where the processed data should be saved and if multiple 
#' codes are used for \code{gene} for the same forward primer, these are 
#' actually ignored.
#' 
#' When relevant, a warning is reported and a text file with the sequence IDs 
#' that had multiple hits in the preliminary search for the forward primer is
#' saved.
#' 
#' After reads are separated based on the forward primer, indexes and primers 
#' are removed and processed samples are written to fastq files. Only reads 
#' where both primers (forward and reverse) and indexes (if there is information
#' for both in \code{info.file}) were found are retained.
#' 
#' @param info.file Fully qualified name (i.e. the complete path) of the CSV 
#'   file with the information needed on primers, indexes etc. (See details)
#' @param sample.IDs A character vector with the name of the column in info.file
#'   containing the sample IDs
#' @param Fprimer,Rprimer A character vector with the name of the column in 
#'   info.file containing the forward and reverse primer sequence, respectively
#' @param primer.mismatch The maximum number of primer mismatch
#' @param Find, Rind A character vector with the name of the column in info.file
#'   containing the forward and reverse index sequence respectively
#' @param index.mismatch The maximum number of index mismatch
#' @param gene A character vector with the name of the column in info.file 
#'   containing the name of the gene or other group idenifiers (see details)
#' @param dir.out The directory where to save the results. If NULL (default) 
#'   then it will be set the same location where the input data was located
#' @return  A fastq file with the reads that were retained after removing the 
#'   indexes (with the suffix "_IndRm") and after removing the primers (with the
#'   suffix "_Ind_primerRm", in a folder named "Final") were end adaptor was 
#'   found (and removed) saved. A list with the total number of reads that were 
#'   processed and retained is also returned.
#'   
#'   When relevant, a text file with the sequence IDs that had multiple hits in 
#'   the preliminary search for the forward primer.
#' @export
deconv <- function(fn, nRead=1e8, info.file, 
                   sample.IDs, Fprimer, Rprimer, primer.mismatch=0,
                   Find, Rind, index.mismatch=0,
                   gene, dir.out=NULL) {
  library(ShortRead)
  
  info_table <- read.csv(info.file)
  primers <- unique(info_table[, Fprimer])
  if (is.null(dir.out)) dir.out  <- dirname(fn)
  nRetained <- 0
  stream <- FastqStreamer(fn, nRead)
  on.exit(close(stream))
  nIndRet <- 0
  nPrimRet <- 0
  while (length(fq <- yield(stream))) {
    for(primer in primers) {
      seqs <- sread(fq)
      qual <- quality(fq)
      qual <- quality(qual)
      sel <- info_table[, Fprimer] == primer
      sub_info_table <- info_table[sel, ]
      gene.out <- paste(dir.out, sub_info_table[, gene][1], sep="/") 
      dir.create(path=gene.out, 
                 showWarnings=FALSE, recursive=TRUE)
      
      
      # Search for earch primers and identify gene
      primer_hits <- vmatchPattern(pattern=DNAString(primer), 
                                   subject=seqs,
                                   max.mismatch=primer.mismatch, 
                                   min.mismatch=primer.mismatch,
                                   with.indels=FALSE, fixed=FALSE,
                                   algorithm="auto")
      
      if(sum(elementLengths(primer_hits) > 1) > 0) {
        message(cat(paste("Note: some reads in", fn, 
                          "have more than one match with the foward primer:", 
                          primer, "Gene:",
                          sub_info_table[, gene][1], 
                          sep="\n")))
        mhits <- paste(gene.out, 
                       paste0(sub_info_table[, gene][1], "_primer_mult_hits.txt"), 
                       sep="/")
        message(cat(paste("Details of reads with multiple primer hits are saved in the txt file:",
                          mhits, sep="\n")))
        capture.output(id(fq)[elementLengths(primer_hits) > 1], file=mhits)
      }
      retain <- as.logical(elementLengths(primer_hits))
      seqs <- seqs[retain]
      qual <- qual[retain]
      ids <- id(fq)[retain]
      
      for(row in 1:dim(sub_info_table)[1]) {
        
        # Removing indexes
        trimCoords <- trimLRPatterns(Lpattern=DNAString(sub_info_table[row, Find]), 
                 Rpattern=reverseComplement(DNAString(sub_info_table[row, Rind])), 
                 subject=seqs,
                 max.Lmismatch=index.mismatch, max.Rmismatch=index.mismatch,
                 with.Lindels=FALSE, with.Rindels=FALSE,
                 Lfixed=TRUE, Rfixed=TRUE, ranges=TRUE)
        retain <- width(trimCoords) == width(seqs) - 
          (nchar(as.vector(sub_info_table[row, Find])) + 
             nchar(as.vector(sub_info_table[row, Rind])))
        seqs_NoInd <- DNAStringSet(seqs, start=start(trimCoords), 
                                   end=end(trimCoords))
        qual_NoInd <- BStringSet(qual, start=start(trimCoords), 
                                 end=end(trimCoords))
        qual_NoInd <- SFastqQuality(qual_NoInd) 
        trimmed <- ShortReadQ(sread=seqs_NoInd, quality=qual_NoInd, id=ids)
        trimmed <- trimmed[retain]
        nIndRet <- nIndRet + length(trimmed)
        if(length(trimmed) > 0) {
          fname <- paste(gene.out, 
                         paste0(sub_info_table[row, sample.IDs], "_IndRm.fastq.gz"), 
                         sep="/")
          writeFastq(trimmed, fname, mode="a")
        }
        
        # Removing primers
        seqs_NoInd <- sread(trimmed)
        qual_NoInd <- quality(trimmed)
        qual_NoInd <- quality(qual_NoInd)
        trimCoords <- trimLRPatterns(Lpattern=DNAString(sub_info_table[row, Fprimer]), 
                 Rpattern=reverseComplement(DNAString(sub_info_table[row, Rprimer])), 
                 subject=seqs_NoInd,
                 max.Lmismatch=primer.mismatch, max.Rmismatch=primer.mismatch,
                 with.Lindels=FALSE, with.Rindels=FALSE,
                 Lfixed=FALSE, Rfixed=FALSE, ranges=TRUE)
        retain <- width(trimCoords) == width(seqs_NoInd) - 
          (nchar(as.vector(sub_info_table[row, Fprimer])) + 
             nchar(as.vector(sub_info_table[row, Rprimer])))
        seqs_NoPrim <- 
          DNAStringSet(seqs_NoInd, start=start(trimCoords), end=end(trimCoords))
        qual_NoPrim <- 
          BStringSet(qual_NoInd, start=start(trimCoords), end=end(trimCoords))
        qual_NoPrim <- SFastqQuality(qual_NoPrim) 
        trimmed <- 
          ShortReadQ(sread=seqs_NoPrim, quality=qual_NoPrim, id=id(trimmed))
        trimmed <- trimmed[retain]
        nPrimRet <- nPrimRet + length(trimmed)
        if(length(trimmed) > 0) {
          dir.create(path=paste(gene.out, "Final", sep="/"),
                     showWarnings=FALSE, recursive=TRUE)
          fname <- paste(gene.out, "Final",
                         paste0(sub_info_table[row, sample.IDs], "_Ind_primerRm.fastq.gz"), 
                         sep="/")
          writeFastq(trimmed, fname, mode="a")
        }
        
      }
    }
    
  }
  nInitial <- stream$status()["total"]
  message(paste("Processed", nInitial, "reads - retained", nIndRet, 
                "after having removed the index."))
  message(paste(nPrimRet, "were retained after removing the primers."))
  return(list(Read=unname(nInitial), IndexRm=nIndRet, PrimerRm=nPrimRet))
}
