#' Remove end adapter
#' 
#' This function is used to remove end adapters, starting from a fastq file.
#' 
#' As mentioned in the general description of this package, most functions are 
#' tailored to Illumina architecture. \code{rmEndAdapter} was developed to 
#' remove the P7 adapter at the end of single-reads. However, it can be actually
#' used to remove any 'tail'. Reads are trimmed at the first position of the 
#' match with the passed \code{EndAdapter} pattern.  The sequence of the 
#' \code{EndAdapter} is passed (as character vector) in a 5' to 3' direction and
#' it is internally reversed and complemented. Other than the sequence, it is 
#' possible to pass the character vector "P7" or "P7_last10". With the first, 
#' the sequence of the P7 adapter is selected (CAAGCAGAAGACGGCATACGAGAT). With 
#' the latter a partial match is searched for (the last 10 bp: CATACGAGAT).
#' 
#' Matches are searched using \code{\link[Biostrings]{vmatchPattern}}, with 
#' \code{adapter.mismatch} used for max.mismatch (min.mismatch is fixed to 
#' zero).
#' 
#' The search is conducted with \code{fixed=TRUE}, which means (from Biostring):
#' "an IUPAC ambiguity code in the pattern can only match the same code in the 
#' subject, and vice versa".
#' 
#' @param fn Fully qualified name (i.e. the complete path) of the fastq file
#' @param nRead The number of bytes or characters to be read at one time. See 
#'   \code{\link[ShortRead]{FastqStreamer}} for details
#' @param EndAdapter A character vector with the sequence of the end adapter,
#'   "P7" or "P7_last10" (See details)
#' @param adapter.mismatch The maximum number of allowed mismatch (See details)
#' @export
#' @return A fastq file with the reads where the end adapter was found (and 
#'   removed) saved in the same location where the input data was located. The 
#'   file is named with the suffix "_EndAdRm". A list with the total number of
#'   reads that were processed and retained is also returned.
#'   
rmEndAdapter <- function(fn, nRead=1e8, EndAdapter="P7_last10", adapter.mismatch=0) {
  if (!requireNamespace("ShortRead", quietly = TRUE)) {
    stop("Package 'ShortRead' needed for this function to work. Please install it 
         either manually or using the function amplicR::setup().",
         call. = FALSE)
  }
  #----------------------------------------------------------------------------#
  # Helper functions
  #----------------------------------------------------------------------------#
  
  putzero <- function(x){
    if(is.null(x)) x <- 0 else x <- x - 1
    return(x)
  }
  #----------------------------------------------------------------------------#
  
  if(EndAdapter == "P7") EndAdapter <- "CAAGCAGAAGACGGCATACGAGAT"
  if(EndAdapter == "P7_last10") EndAdapter <- "CATACGAGAT"  
  
  stream <- ShortRead::FastqStreamer(fn, nRead)
  on.exit(close(stream))
  nRetained <- 0
  
  while (length(fq <- ShortRead::yield(stream))) {
    seqs <- ShortRead::sread(fq)
    qual <- Biostrings::quality(fq)
    qual <- Biostrings::quality(qual)
    P7_hits <- Biostrings::vmatchPattern(
                             pattern=Biostrings::reverseComplement(
                                             Biostrings::DNAString(EndAdapter)), 
                             subject=seqs,
                             max.mismatch=adapter.mismatch, 
                             min.mismatch=0,
                             with.indels=FALSE, fixed=TRUE,
                             algorithm="auto")
    
    if(sum(S4Vectors::elementNROWS(P7_hits) > 1) > 0) {
      message(cat(paste("Note: some reads in", fn, "have more than one match with:", 
                        EndAdapter, "The first match from the left is used",
                        sep="\n")))
    }
    ends <- Biostrings::startIndex(P7_hits)
    ends <- lapply(ends, "[", 1)
    ends <- lapply(ends, putzero)
    
    seqs <- Biostrings::DNAStringSet(seqs, start=1, end=unlist(ends))
    qual <- Biostrings::BStringSet(qual, start=1, end=unlist(ends))
    qual <- ShortRead::SFastqQuality(qual) 
    trimmed <- ShortRead::ShortReadQ(sread=seqs, quality=qual, 
                                                           id=ShortRead::id(fq))
    trimmed <- trimmed[IRanges::width(trimmed) > 0]
    
    fname <- paste0(substr(fn, start=1, stop=regexpr(".fastq", fn)[1] - 1), 
                    "_EndAdRm.fastq.gz")
    if(file.exists(fname)) message(cat(paste("Note: the file", fname, 
                                    "already exists\nSequences were appended")))
    ShortRead::writeFastq(trimmed, fname, mode="a")
    nRetained <- nRetained + length(trimmed)
  }
  nInitial <- stream$status()["total"]
  message(paste("Processed", nInitial, "- retained", nRetained))
  return(rmEndAd <- list(Read=unname(nInitial), nRetained=nRetained))
}

#' Separate reads by genes and deconvolute them based on barcodes
#' 
#' \code{deconv} takes a fastq file and will search for the forward primer and 
#' use this to separate the reads. That is, different PCR products ('genes') 
#' will be separated based on the forward primer. Within each gene, reads are 
#' then separated based on forward and/or reverse index.  The end products are 
#' several fastq files - one for each samples, in as many folders as how many 
#' gene identifiers were provided with the \code{info.file} - where primers and 
#' indexes were removed.
#' 
#' This function applies only to reads with in-line indexes. That is, where the 
#' architecture of the reads is as follows:
#' 
#' F_index---F_primer---Target_sequence---R_primer---R_index
#' 
#' Note that the P7 adapter can be removed with \code{\link{rmEndAdapter}}.
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
#' reverse primer, the sample IDs, and an identifier of the PCR product being 
#' amplified, typically the gene's name. The column headings where these 
#' information are stored in \code{info.file} are passed with the function 
#' arguments. While it is possible to include other columns where the users can 
#' record additional information, these are effectively ignored. It is not 
#' mandatory to have both, the forward and reverse indexes, but if one is not 
#' used, there is still the need to include a blank column in \code{info.file} 
#' and indicate the column heading. Note that, when importing \code{info.file}, 
#' R will automatically convert illegal characters (e.g. sapces, paranthesis) in
#' dots ('.'), so it is probably safer to only use alpha-numeric characters 
#' and/or dots or underscores ('_') in the function's arguments.
#' 
#' All sequences for indexes and primers are passed (as character vector) in 5' 
#' to 3' direction and are internally reversed and complemented when necessary.
#' 
#' \code{deconv} initially searches for the forward primer and separates the 
#' reads creating (if not existing already) a folder named as for the relevant 
#' information provided in the column \code{gene}, which is typically an 
#' identifier of the targeted genes. It is possible to use the column 
#' \code{gene} to group PCR products in other logical way than genes, but all
#' identical forward primer should have the same \code{gene} information. This
#' is because (for efficiency) \code{deconv} uses only the first line for each
#' forward primer to identify where the processed data should be saved and if
#' multiple codes are used for \code{gene} for the same forward primer, these
#' are actually ignored.
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
#' @param Find,Rind A character vector with the name of the column in info.file 
#'   containing the forward and reverse index sequence respectively
#' @param index.mismatch The maximum number of index mismatch
#' @param gene A character vector with the name of the column in info.file 
#'   containing the name of the gene or other group idenifiers (see details)
#' @param dir.out The directory where to save the results. If NULL (default) 
#'   then it will be set the same location where the input data was located
#' @inheritParams rmEndAdapter
#' @return  A fastq file with the reads that were retained after removing the 
#'   indexes (with the suffix "_IndRm") and after removing the primers (with the
#'   suffix "_Ind_primerRm", in a folder named "Final") were end adapter was 
#'   found (and removed) saved. A list with the total number of reads that were 
#'   processed and retained is also returned.
#'   
#'   When relevant, a text file with the sequence IDs that had multiple hits in 
#'   the preliminary search for the forward primer.
#' @export
deconv <- function(fn, nRead=1e8, info.file, sample.IDs="Sample_IDs", 
                   Fprimer="F_Primer", Rprimer="R_Primer", primer.mismatch=0,
                   Find="F_ind", Rind="R_ind", index.mismatch=0,
                   gene="Gene", dir.out=NULL) {
  if (!requireNamespace("ShortRead", quietly = TRUE)) {
    stop("Package 'ShortRead' needed for this function to work. Please install it 
         either manually or using the function amplicR::setup().",
         call. = FALSE)
  }
  
  
  info_table <- read.csv(info.file)
  primers <- unique(info_table[, Fprimer])
  if (is.null(dir.out)) dir.out  <- dirname(fn)
  nRetained <- 0
  stream <- ShortRead::FastqStreamer(fn, nRead)
  on.exit(close(stream))
  nIndRet <- 0
  nPrimRet <- 0
  while (length(fq <- ShortRead::yield(stream))) {
    for(primer in primers) {
      seqs <- ShortRead::sread(fq)
      qual <- Biostrings::quality(fq)
      qual <- Biostrings::quality(qual)
      sel <- info_table[, Fprimer] == primer
      sub_info_table <- info_table[sel, ]
      gene.out <- paste(dir.out, sub_info_table[, gene][1], sep="/") 
      dir.create(path=gene.out, 
                 showWarnings=FALSE, recursive=TRUE)
      
      
      # Search for earch primers and identify gene
      primer_hits <- Biostrings::vmatchPattern(
                                   pattern=Biostrings::DNAString(primer), 
                                   subject=seqs,
                                   max.mismatch=primer.mismatch, 
                                   min.mismatch=0,
                                   with.indels=FALSE, fixed=FALSE,
                                   algorithm="auto")
      
      if(sum(S4Vectors::elementNROWS(primer_hits) > 1) > 0) {
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
        capture.output(
          ShortRead::id(fq)[S4Vectors::elementNROWS(primer_hits) > 1], file=mhits)
      }
      retain <- as.logical(S4Vectors::elementNROWS(primer_hits))
      seqs <- seqs[retain]
      qual <- qual[retain]
      ids <- ShortRead::id(fq)[retain]
      
      for(row in 1:dim(sub_info_table)[1]) {
        
        # Removing indexes
        trimCoords <- Biostrings::trimLRPatterns(
                      Lpattern=Biostrings::DNAString(sub_info_table[row, Find]), 
                      Rpattern=Biostrings::reverseComplement(
                              Biostrings::DNAString(sub_info_table[row, Rind])), 
                 subject=seqs,
                 max.Lmismatch=index.mismatch, max.Rmismatch=index.mismatch,
                 with.Lindels=FALSE, with.Rindels=FALSE,
                 Lfixed=TRUE, Rfixed=TRUE, ranges=TRUE)
        retain <- IRanges::width(trimCoords) == IRanges::width(seqs) - 
          (nchar(as.vector(sub_info_table[row, Find])) + 
             nchar(as.vector(sub_info_table[row, Rind])))
        seqs_NoInd <- Biostrings::DNAStringSet(seqs, 
                                          start=BiocGenerics::start(trimCoords), 
                                          end=BiocGenerics::end(trimCoords))
        qual_NoInd <- Biostrings::BStringSet(qual, 
                                          start=BiocGenerics::start(trimCoords), 
                                          end=BiocGenerics::end(trimCoords))
        qual_NoInd <- ShortRead::SFastqQuality(qual_NoInd) 
        trimmed <- ShortRead::ShortReadQ(sread=seqs_NoInd, quality=qual_NoInd, id=ids)
        trimmed <- trimmed[retain]
        nIndRet <- nIndRet + length(trimmed)
        if(length(trimmed) > 0) {
          fname <- paste(gene.out, 
                         paste0(sub_info_table[row, sample.IDs], "_IndRm.fastq.gz"), 
                         sep="/")
          ShortRead::writeFastq(trimmed, fname, mode="a")
        }
        
        # Removing primers
        seqs_NoInd <- ShortRead::sread(trimmed)
        qual_NoInd <- Biostrings::quality(trimmed)
        qual_NoInd <- Biostrings::quality(qual_NoInd)
        trimCoords <- Biostrings::trimLRPatterns(
                   Lpattern=Biostrings::DNAString(sub_info_table[row, Fprimer]), 
                   Rpattern=Biostrings::reverseComplement(
                           Biostrings::DNAString(sub_info_table[row, Rprimer])), 
                 subject=seqs_NoInd,
                 max.Lmismatch=primer.mismatch, max.Rmismatch=primer.mismatch,
                 with.Lindels=FALSE, with.Rindels=FALSE,
                 Lfixed=FALSE, Rfixed=FALSE, ranges=TRUE)
        retain <- IRanges::width(trimCoords) == IRanges::width(seqs_NoInd) - 
          (nchar(as.vector(sub_info_table[row, Fprimer])) + 
             nchar(as.vector(sub_info_table[row, Rprimer])))
        seqs_NoPrim <- 
          Biostrings::DNAStringSet(seqs_NoInd, 
                                          start=BiocGenerics::start(trimCoords), 
                                          end=BiocGenerics::end(trimCoords))
        qual_NoPrim <- 
          Biostrings::BStringSet(qual_NoInd, 
                                         start=BiocGenerics::start(trimCoords), 
                                         end=BiocGenerics::end(trimCoords))
        qual_NoPrim <- ShortRead::SFastqQuality(qual_NoPrim) 
        trimmed <- 
          ShortRead::ShortReadQ(sread=seqs_NoPrim, quality=qual_NoPrim, 
                                                      id=ShortRead::id(trimmed))
        trimmed <- trimmed[retain]
        nPrimRet <- nPrimRet + length(trimmed)
        if(length(trimmed) > 0) {
          dir.create(path=paste(gene.out, "Final", sep="/"),
                     showWarnings=FALSE, recursive=TRUE)
          fname <- paste(gene.out, "Final",
                         paste0(sub_info_table[row, sample.IDs], 
                                "_Ind_primerRm.fastq.gz"), 
                         sep="/")
          ShortRead::writeFastq(trimmed, fname, mode="a")
        }
        
      }
    }
    
  }
  nInitial <- stream$status()["total"]
  message(paste("Processed", nInitial, "reads - retained", nIndRet, 
                "after having removed the index(es)."))
  message(paste(nPrimRet, "reads were retained after removing the primers."))
  return(list(Read=unname(nInitial), IndexRm=nIndRet, PrimerRm=nPrimRet))
}

#' From raw data to data.proc()
#' 
#' This function is a wrapper for \code{\link{rmEndAdapter}}, 
#' \code{\link{deconv}} and \code{\link{data.proc}}. It takes in a raw fastq 
#' file, removes the end adapter, separates the reads based on their forward 
#' primers. Within each of the iodentified group, separates the reads based on 
#' barcodes (indexes) and eventually calls \code{\link{data.proc}} to process 
#' (quality checking, denoising and chimeras filtering) the retained data from 
#' the NGS run.
#' 
#' Note that the amplicon size for \code{\link{data.proc}} is obtained from the 
#' comma delimited file \code{info.file}, searching in the column with the 
#' heading indicated in \code{amplic.size}. Zeros can be used in this column if 
#' no truncation is wanted. For each entry in the column 
#' indicated with the argument \code{gene}, the function will use the first 
#' entry found in \code{amplic.size} for the relevant \code{gene}. If the same 
#' gene identifier is used for multiple forward primers (see documentation for 
#' the \code{\link{deconv}} to see how multiple PCR product can be grouped 
#' together using the \code{gene} column), then these have to have all the same 
#' amplicon length to use \code{raw2data.proc}, otherwise the three functions 
#' (\code{\link{rmEndAdapter}}, \code{\link{deconv}} and 
#' \code{\link{data.proc}}) need to be called manually, rather than with 
#' \code{raw2data.proc}. This is because \code{\link{data.proc}} needs a 
#' gene-specific amplicon length to correctly process the data.
#' 
#' By default, \code{dir.out} is set to the location where the input file is and 
#'   \code{verbose=FALSE} for \code{\link{data.proc}}.
#' 
#' 
#' Please, see documentations for each functions for more information.
#' 
#' @inheritParams rmEndAdapter
#' @inheritParams deconv
#' @inheritParams data.proc
#' @param amplic.size A character vector with the name of the column in 
#'   info.file containing the amplicon size of the PCR product
#' @return A list that has for elements the output of \code{\link{data.proc}} 
#'   for each PCR product
#'   
#'   Also, in addition to the output files described in the documentations for
#'   \code{\link{rmEndAdapter}}, \code{\link{deconv}} and
#'   \code{\link{data.proc}}, a text file, named "summary_nReads.txt" is saved
#'   in the same location where the raw data are, summarising the number of
#'   reads retained in each step of the analysis
#' @export
#' @import data.table
raw2data.proc <- function(fn, nRead=1e8, EndAdapter="P7_last10", 
                          adapter.mismatch=0, info.file, sample.IDs="Sample_IDs", 
                          Fprimer="F_Primer", Rprimer="R_Primer", 
                          primer.mismatch=0, Find="F_ind", Rind="R_ind", 
                          index.mismatch=0, gene="Gene",  
                          amplic.size="Amplicon", truncQ=2, qrep=FALSE,
                       dada=TRUE, pool=FALSE, plot.err=FALSE, chim=TRUE,
                       orderBy="abundance") {
  #----------------------------------------------------------------------------#
  # Helper functions
  #----------------------------------------------------------------------------#
extract.sums <- function(ldproc, el)  {
  l <- lapply(ldproc, extr <- function(dproc, el) {
                                      if(el == "nFiltered") i <- 1
                                      if(el == "nDerep") i <- 2
                                      if(el == "nSeq") i <- length(dproc)
                                      r <- if(is.data.table(dproc$lsummary[[i]])) {
                                        sum(dproc$lsummary[[i]][, el, with=FALSE], 
                                            na.rm=TRUE)
                                      } else {
                                        sum(dproc$lsummary[[i]][el], na.rm=TRUE)
                                      }
                                      return(r)
                                    },
              el)
  s <- sum(unlist(l))
  return(s)
}
  #----------------------------------------------------------------------------#
  info_table <- read.csv(info.file)
  primers <- unique(info_table[, Fprimer])
  genes <- unique(info_table[, gene])
  if(length(primers) != length(genes)) 
    stop(paste("Detected a different number of unique forward primers and",
               gene))
  
  rme <- rmEndAdapter(fn, nRead, EndAdapter, adapter.mismatch)
  
  fn_Endrm <- paste0(substr(fn, start=1, stop=regexpr(".fastq", fn)[1] - 1), 
                        "_EndAdRm.fastq.gz")
  dec <- deconv(fn_Endrm, nRead, info.file, 
          sample.IDs, Fprimer, Rprimer, primer.mismatch,
          Find, Rind, index.mismatch,
          gene)
  
  path.results <- paste(dirname(fn), genes, "Final", sep="/")
  names(path.results) <- genes
  ldproc <- list()
  for(g in genes) {
    sel <- info_table[, gene] == g
    bp <- info_table[sel, amplic.size][1]
  
      txt <- capture.output(
      ldproc[[g]] <- data.proc(dir.in=path.results[g], bp=bp, truncQ=truncQ, 
                            qrep=qrep, dada=dada, pool=pool, plot.err=plot.err, 
                            chim=chim, orderBy=orderBy, verbose=FALSE)
    )
  }
  
  writeLines(c(paste("The number of reads found in", fn, "was", rme[[1]]), 
               paste("The end adapter was found and removed in", rme[[2]], "reads"),
               paste("The index(es) were found and removed in", dec[[2]], "reads"),
               paste("Primers were found and removed in ", dec[[3]], "reads"),
               paste("The number of reads retained after applying the quality filter was",
                     extract.sums(ldproc, "nFiltered")),
               paste("The number of unique reads retained across all samples was", 
                     extract.sums(ldproc, "nDerep")),
               paste("The number of unique reads retained across all samples at completion of data.proc() was",
                     extract.sums(ldproc, "nSeq")),
               "More details are in the files 'data.proc.summary.csv' in each folder"),
             con=paste(dirname(fn), "summary_nReads.txt", sep="/"))
  return(ldproc)
}

