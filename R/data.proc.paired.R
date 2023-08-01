#' Data processing
#' 
#' \code{data.proc.paired} is a function to process (quality checking, error and 
#' chimeras filtering) data from a NGS run after these have been deconvoluted.
#' 
#' \code{data.proc.paired} locates the .fastq files (can be compressed) in the 
#' directory indicated in \code{dir.in}. If the directory path is not provided, 
#' this will be selected using an interactive window.
#' 
#' This is the first implementation to deal with paired-end reads and assumes 
#' that adapters, primers 
#' and indexes have been already removed and that each file represents a sample.
#' 
#' The \code{data.proc.paired} pipeline is as follows: fastq files are read in. A 
#' filter is applied to truncate reads at the first instance of a quality score 
#' less than \code{truncQ}, remove reads  that are of low quality (currently the
#' threshold is hard-coded and reads are discarded if the expected errors is
#' higher than 3 - from documentation in the R package \code{dada2}, the
#' expected errors are calculated from the nominal definition of the quality
#' score: EE = sum(10^(-Q/10)) - and remove reads that (after truncation) do not
#' match the target length. A quality report can be (optionally) generated with the
#' R package \code{ShortReads} to verify the quality of the reads retained
#' after this step. Reads are then dereplicated. Optionally, the dada (Callahan
#' et al 2015) algoritm is applied and bimeras are searched and removed with
#' default settings of the relative functions in the package \code{dada2}. The
#' sequences that were retained at completion of \code{data.proc.paired} are saved in
#' fasta files in the subfolder "Final_seqs" and a .csv with a summary of the
#' number of reads that have been retained in each step is also written. These
#' two outputs are also returned at the end of the function.
#' 
#' The sequence data handling is done by using functionalities from the packages
#' \code{dada2} and \code{ShortRead}, so make sure to cite them (in addition to 
#' \code{amplicR} of course!) if you report your results in a paper or report.
#' 
#' @param dir.in The directory where the fastq files are located. It needs a 
#'   vector of length=2 If NULL 
#'   (default) an interactive window is used to select a folder
#' @param dir.out The directory where to save the results. If NULL (default) 
#'   then \code{dir.out <- dir.in}
#' @param bp An integer indicating the expected length (base-pairs) of the
#'   reads. If zero (default) no truncation is applied
#' @param truncQ Truncate reads at the first instance of a quality score less 
#'   than or equal to truncQ when conducting quality filtering. See
#'   \code{\link[dada2]{fastqFilter}} for details
#' @param qrep Logical. Should the quality report be generated? (default 
#'   \code{FALSE})
#' @param dada Logical. Should the dada analysis be conducted? (default 
#'   \code{TRUE})
#' @param pool Logical. Should samples be pooled together prior to sample 
#'   inference? (default \code{FALSE}). See \code{\link[dada2]{dada}} for
#'   details
#' @param plot.err Logical. Whether error rates obtained from \code{dada} should
#'   be plotted
#' @param chim Logical. Should the bimera search and removal be performed? 
#'   (default \code{TRUE})
#' @param minOverlap When merging paired ends, what is the minimum of basis that 
#'   should overlap
#' @param maxMismatch When merging paired ends, the maximum number of allowed mismatches
#' @param {trimOverhang} {When merging paired ends, whether the overhanging 
#'   extreme of the sequence should be trimmed off}
#' @param orderBy Character vector specifying how the returned sequence table 
#'   should be sorted. Default "abundance". See 
#'   \code{\link[dada2]{makeSequenceTable}} for details
#' @param verbose Logical. Whether information on progress should be outputted 
#'   (default: TRUE)
#' @return Return a list with several elements:
#'   
#'   \itemize{ \item $luniseqsFinal: A list with unique sequences (names) that
#'   were retained at completion of \code{data.proc} and their abundance
#'   (values). \item $lsummary: A list where each element is a summary of the
#'   number of reads that were retained in each step. This can be converted in a
#'   \code{data.frame} by running the following
#'   
#'   \code{summary <- plyr::join_all(lsummary, by="Sample", type="left")} \item
#'   $stable: The sequence table \item $seq_list: The sequences and matching 
#'   sequence IDs \item $call: The function call}
#'   
#'   Several files are also returned, these include: \itemize{ \item
#'   Seq_table.csv Sequence table \item Seq_list.csv List of sequences and their
#'   matching IDs \item data.proc.summary.csv A summary of the number of reads
#'   that were retained in each step \item data.proc.rda R data file containing
#'   the list returned by data.proc (see above) }
#'   
#' @references Benjamin J Callahan, Paul J McMurdie, Michael J Rosen, Andrew W 
#'   Han, Amy J Johnson, Susan P Holmes (2015). DADA2: High resolution sample 
#'   inference from amplicon data.
#' @import ggplot2
#' @import data.table
#' @seealso \code{\link[dada2]{dada}}, \code{\link[dada2]{makeSequenceTable}}
#' @export
#' @examples 
#' # Select the directory where the example data are stored
#' example.data <- system.file("extdata", "HTJ", package="amplicR")
#' # Select a temporary directory where to store the outputs
#' out <- tempdir()
#' 
#' HTJ.test <- data.proc(example.data, out, bp=140)
#' 
#' # To clean up the temp directory
#' unlink(out, recursive=TRUE)

data.proc.paired <- function(dir.in=NULL, dir.out=NULL, bp=0, truncQ=2, qrep=FALSE,
                      dada=TRUE, pool=FALSE, plot.err=FALSE, chim=TRUE, 
                      minOverlap = 12,
                      maxMismatch = 0,
                      trimOverhang = FALSE,
                      orderBy="abundance", paired = TRUE, verbose=TRUE) {
  
#----------------------------------------------------------------------------#
  if (!requireNamespace("dada2", quietly = TRUE)) {
    stop("Package 'dada2' needed for this function to work. Please install it 
         either manually or using the function amplicR::setup().",
         call. = FALSE)
  }
#----------------------------------------------------------------------------#
# Helper functions
#----------------------------------------------------------------------------#


# pass dada list, extracts unique reads with dada2::getUniques() and returns
# only those that are not chimeras.

# NOTE: the same process can be applied to dada elements (for example data.frame)
# to maintain additional info if needed later on
rm.chim<-function(i, dada_el, chim_el) {
  no_chim <- dada2::getUniques(dada_el[[i]])[!chim_el[[i]]]
  return(no_chim)
}

nChim <- function(chim_el) {
  nSeqs <- sum(chim_el)
  return(nSeqs)
}

dp_extract <- function(nm, derep) {
  seqs <- dada2::getUniques(derep[[nm]])
  return(seqs)
}

#----------------------------------------------------------------------------#
call <- sys.call(1)

# if(is.null(dir.in)) {
#   dir.in <- choose.dir(caption="Please, select the directory where the fastq
#                        files are located")
# }

#if(paired == TRUE) 
  if(length(dir.in) !=2) 
  stop("It seems you are trying to process paired end reads but you didn't provide two paths for dir.in")

if(is.null(dir.out)) {
  #if(paired == TRUE) {
    stop("dir.out must be provided when prosessing paired end reads")
  #} else {
  #  dir.out <- dir.in
  }
    
dir.create(dir.out, showWarnings=FALSE, recursive=TRUE)
filt_fold <- "Filtered_seqs"
dir.create(paste(dir.out, filt_fold, sep="/"), showWarnings=FALSE, recursive=TRUE)

fns <- vector("list", length = length(dir.in))
fastqs <- vector("list", length = length(dir.in))
filtRs <- vector("list", length = length(dir.in))
for(d in seq_along(dir.in)) {
  fns[[d]] <- list.files(path=dir.in[d])
  fastqs[[d]] <- fns[[d]][grepl(".fastq.{,3}$", fns[[d]])]
  if(length(fastqs[[d]]) == 0) stop(paste("There are no files in", dir.in[d],
                                     "with either fastq or fastq.gz extension"))
  filtRs[[d]] <- paste(dir.out,
                filt_fold,
                 sapply(fastqs[[d]],
                        sub,
                        pattern="\\.fastq.{,3}$",
                        replacement=paste0("R", d, "_filt.fastq.gz")),
                sep="/"
                 )
}

if(length(fastqs[[1]]) == length(fastqs[[2]])) {
  n <- length(fastqs[[1]]) 
  message(paste("Found", n, "fasta files in each input dir.in folders"))
} else {
  stop("Different number of fasta files in the two dir.in folders")
}

#### Filter ####
sample_names <- sub(".fastq.{,3}$", "", fastqs[[1]])
fastqsPairs <- mapply(c, file.path(dir.in[1], fastqs[[1]]), 
                      file.path(dir.in[2], fastqs[[2]]))
filtRsPairs <- mapply(c, filtRs[[1]], filtRs[[2]])

for(i in seq_len(n)) {
  suppressWarnings(dada2::fastqPairedFilter(fastqsPairs[, i], 
                                       filtRsPairs[, i], 
                                       maxN=0, 
                                       maxEE=3,
                                       truncQ=truncQ, 
                                       truncLen=bp, 
                                       compress=TRUE, 
                                       verbose=FALSE))
}

filtRs <- list.files(path=file.path(dir.out, filt_fold), full.names=TRUE)
sample_names_fil <- sub("R[1-2]_filt.fastq.gz", "", 
                        sapply(filtRs, basename, USE.NAMES=FALSE))

if(qrep == TRUE) browseURL(report(qa(filtRs)))

#### Dereplicate ####
# Fwd
derepReadsF <- lapply(filtRs[grep("R1", x = filtRs)], dada2::derepFastq, verbose=FALSE)
names(derepReadsF) <- sample_names_fil
fnSeqsF <- unlist(lapply(derepReadsF, getnFiltered))
if(length(fnSeqsF) == 0) stop(paste("There are no sequences that passed the filter in", 
                                   file.path(dir.out, filt_fold)))
unSeqsF <- unlist(lapply(derepReadsF, getnUniques))

# Rev
derepReadsR <- lapply(filtRs[grep("R2", x = filtRs)], dada2::derepFastq, verbose=FALSE)
names(derepReadsR) <- sample_names_fil
fnSeqsR <- unlist(lapply(derepReadsR, getnFiltered))
if(length(fnSeqsR) == 0) stop(paste("There are no sequences that passed the filter in", 
                                    file.path(dir.out, filt_fold)))
unSeqsR <- unlist(lapply(derepReadsR, getnUniques))

# summary
lsummary <- list()
lsummary[[1]] <- merge(data.table(Sample=sample_names), 
                                   data.table(Sample=sample_names_fil, 
                                              nFilteredF=fnSeqsF,
                                              nFilteredR=fnSeqsR), 
                                   all.x=TRUE,
                                   by="Sample")

derepReads <- list(derepReadsF, derepReadsR)


lsummary[[2]] <- data.frame(Sample=sample_names_fil, nDerepF=unSeqsF, nDerepR=unSeqsR)
el <- 2

#### dada ####
if(dada == TRUE) {
  for(d in 1:2) {
  dadaReads[[1]] <- dada2::dada(derepReads[[d]], err=dada2::inflateErr(dada2::tperr1,3),
                    errorEstimationFunction=dada2::loessErrfun,
                    selfConsist=TRUE, pool=pool)
  if(plot.err == TRUE) {
    pdf(file = file.path(dir.out, paste0("R", d, "_Plot_ErrorRates.pdf")))
    if(length(derepReads[[d]]) > 1) {
      for (i in seq_along(dadaReads[[d]])) {
        p <- dada2::plotErrors(dadaReads[[d]][[i]], nominalQ=TRUE)
        print(p + ggtitle(names(dadaReads[[d]][i])))
      }
    } else {
      p <- dada2::plotErrors(dadaReads[[d]], nominalQ=TRUE)
      print(p + ggtitle(names(derepReads[[d]])))
    }
    dev.off()
  }
  
  if(length(derepReads[[d]]) > 1) {
    nDenoised[[d]] <- sapply(dadaReads[[d]], ndada)
    lda[[d]] <- lapply(dadaReads[[d]], dada2::getUniques)
  } else {
    nDenoised[[d]] <- length(dada2::getUniques(dadaReads[[d]]))
    names(nDenoised[[d]]) <- names(derepReads[[d]])
    lda[[d]] <- list(dada2::getUniques(dadaReads[[d]]))
    names(lda[[d]]) <- names(derepReads[[d]])
  }
  }
  
el <- el + 1
  lsummary[[el]] <- data.frame(Sample=names(nDenoised[[1]]), 
                               nDenoisedF=nDenoised[[1]],
                               nDenoisedR=nDenoised[[2]])

  #### Merge ####
  mergePairedEnds <- mergePairs(dadaF = dadaReads[[1]], derepF = derepReads[[1]],
             dadaR = dadaReads[[2]], derepR = derepReads[[2]],
             minOverlap = 12,
             maxMismatch = 0,
             returnRejects = FALSE,
             propagateCol = character(0),
             justConcatenate = FALSE,
             trimOverhang = FALSE,
             verbose = verbose)
  
  el <- el + 1
  lsummary[[el]] <- mergePairedEnds
  }

  if(verbose) {
    message("Number of sequenced retained after sample inference with dada2
          (Callahan et al 2015):")
    message(cat(nDenoised, "\n"))
  }
  
  
  
#### Bimera ####
if(chim == TRUE)  {
  if(dada == TRUE) {
    uniqSeqs <- mergePairedEnds$abundance
    names(uniqSeqs) <- mergePairedEnds$sequence
    if(length(uniqSeqs) > 1) {
      single <- FALSE
      bimReads <- lapply(uniqSeqs, dada2::isBimeraDenovo, verbose=FALSE)
      names(bimReads) <- names(uniqSeqs)
    } else {
      single <- TRUE
      bimReads <- dada2::isBimeraDenovo(uniqSeqs, verbose=FALSE)
    }
  } else {
    if(length(derepReadsF) > 1) {
      single <- FALSE
    } else {
      single <- TRUE
    }
    warning("The data.proc.paired function is not fully developed to summarise results 
            when dada == FALSE and results are reported only for the Forward reads")
    bimReads <- lapply(derepReadsF, dada2::isBimeraDenovo, verbose=FALSE)
    names(bimReads) <- names(derepReadsF)
  }
  
  if(verbose) {
    message("Proportion of bimeras found in each sample")
      message(cat(round(sapply(bimReads, mean), digits=2), "\n"))
  }
  
    nChimeras <- sapply(bimReads, nChim)

  if(verbose) {
    message("Number of bimeras found in each sample:")
    message(cat(nChimeras, "\n"))
  }
  
  el <- el + 1
  lsummary[[el]] <- data.frame(Sample=names(nChimeras), nChimeras)
  if(dada == TRUE) {
    if(single) {
      no_chim <- list(dada2::getUniques(dadaReads)[!bimReads])
      names(no_chim) <- names(derepReads)
    } else {
      no_chim <- lapply(seq_along(dadaReads), rm.chim, dadaReads, bimReads)
      names(no_chim) <- names(bimReads)
    }
  } else {
    no_chim <- lapply(seq_along(derepReads), rm.chim, derepReads, bimReads)
    names(no_chim) <- names(bimReads)
  }
}

  #### Reporting ####
if(chim == TRUE) {
  luniseqsFinal <- no_chim
  nSeq <- sapply(luniseqsFinal, length)
} else {
  if(dada == TRUE) {
    luniseqsFinal <- lda
    nSeq <- sapply(luniseqsFinal, length)
  } else {
    luniseqsFinal <- derepReads
    nSeq <- unSeqs
  }
}

stable <- dada2::makeSequenceTable(luniseqsFinal, orderBy=orderBy)
seqs <- colnames(stable)
seq_names <- paste0("seq", 1:dim(stable)[2])
colnames(stable) <- seq_names
seq_list <- data.frame(seq_names, sequence=seqs)
write.csv(stable, file=paste(dir.out, "Seq_table.csv", sep="/"))
write.csv(seq_list, file=paste(dir.out, "Seq_list.csv", sep="/"), row.names=FALSE)

el <- el + 1
lsummary[[el]] <- data.frame(Sample=names(nSeq), nSeq)
summary <- plyr::join_all(lsummary, by="Sample", type="left")
write.csv(summary, file=paste(dir.out, "data.proc.summary.csv", sep="/"), 
          row.names=FALSE)

fasta_fold <- "Final_seqs"
fasta_dir <- paste(dir.out, fasta_fold, sep="/")
table2fasta(stable, seq.list=seq_list, dir.out=fasta_dir, verbose)
dproc <- list(luniseqsFinal=luniseqsFinal, 
              lsummary=lsummary, 
              stable=stable, 
              seq_list=seq_list,
              call=call)
save(dproc, file=paste(dir.out, "data.proc.rda", sep="/"))
return(dproc)
}
