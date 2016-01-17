#' Data processing
#' 
#' \code{data.proc} is a function to process (quality checking, error and 
#' chimeras filtering) data from a NGS run after these have been deconvoluted.
#' 
#' \code{data.proc} locates the .fastq files in the directory indicated in 
#' \code{dir.in}. If the directory path is not provided, this will be selected 
#' using an interactive window.
#' 
#' It is currently limited to single-reads and assumes that adapters, primers 
#' and indexes have been already removed and that each file represents a sample.
#' 
#' The  \code{data.proc} pipeline is as follows: fastq files are read in. A 
#' filter is applied to truncate reads at the first instance of a quality score 
#' less than 2, remove reads  that are of low quality (currently the threshold 
#' is hard-coded and reads are discarded if the expected errors is higher than 3
#' - from documentation in the R package \code{dada2}, the expected errors are 
#' calculated from the nominal definition of the quality score: EE = 
#' sum(10^(-Q/10)) - and remove reads that (after truncation) do not match the 
#' target length. A quality report can be (optionally) generated with R package 
#' \code{shortReads} the to verify the quality of the reads retained after this 
#' step. Reads are then dereplicated. Optionally, the dada (Callahan et al 2015)
#' algoritm is applied and bimeras are searched and removed with default 
#' settings of the relative functions in the package \code{dada2}. The sequences
#' that were retained at completion of \code{data.proc} are saved in fasta files
#' in the subfolder "Final_seqs" and a .csv with a summary of the number of 
#' reads that have been retained in each step is also written. These two outputs
#' are also returned at the end of the function.
#' 
#' The sequence data handling is done by using functionalities from the packages
#' \code{dada2} and \code{ShortRead}, so make sure to cite them (in addition to
#' \code{amplicR} of course!) if you report your results in a paper or report.
#' 
#' @param dir.in The directory where the fastq files are located. If NULL 
#'   (default) an interactive window is used to select a folder
#' @param dir.out The directory where to save the results. If NULL (default) 
#'   then \code{dir.out <- dir.in}
#' @param bp An integer indicating the expected length (base-pairs) of the reads
#' @param qrep Logical. Should the quality report be generated? (default 
#'   \code{FALSE})
#' @param dada Logical. Should the dada analysis be conducted? (default 
#'   \code{TRUE})
#' @param chim Logical. Should the bimera search and removal be performed? 
#'   (default \code{TRUE})
#' @return Return a list with two elements and several files (see details). The 
#'   first is a \code{DNAString} object with the sequences that were retained at
#'   completion of \code{data.proc}. The second element is a a list where each 
#'   element is a summary of the number of reads that have been retained in each
#'   step. This can be converted in a \code{data.frame} by running the following
#'   
#'   \code{summary <- plyr::join_all(lsummary, by="Sample", type="left")}
#' @references Benjamin J Callahan, Paul J McMurdie, Michael J Rosen, Andrew W 
#'   Han, Amy J Johnson, Susan P Holmes (2015). DADA2: High resolution sample 
#'   inference from amplicon data.
#' @import ggplot2
#' @import data.table
#' @export

data.proc <- function(dir.in=NULL, dir.out=NULL, bp, qrep=FALSE,
                      dada=TRUE, chim=TRUE) {
#----------------------------------------------------------------------------#
library(dada2)
library(ShortRead)
#----------------------------------------------------------------------------#
# Helper functions
#----------------------------------------------------------------------------#

getnUniques <- function (derep) {
  nUniques <- length(derep$uniques)
  return(nUniques)
}

ndada <- function(dada_el) {
  nSeqs <- length(getUniques(dada_el))
  return(nSeqs)
}

# pass dada list, extracts unique reads with dada2::getUniques() and returns
# only those that are not chimeras.

# NOTE: the same process can be applied to dada elements (for example data.frame)
# to maintain additional info if needed later on
rm.chim<-function(i, dada_el, chim_el) {
  no_chim <- getUniques(dada_el[[i]])[!chim_el[[i]]]
  return(no_chim)
}

nChim <- function(chim_el) {
  nSeqs <- sum(chim_el)
  return(nSeqs)
}

dp_extract <- function(nm, derep) {
  seqs <- getUniques(derep[[nm]])
  return(seqs)
}

make.DNAString <- function(luniseqs_el) {
  DNAstr <- DNAStringSet(names(luniseqs_el))
  return(DNAstr)
}

w.fasta <- function(i, lseqs, sample_names, dir.out, fasta_fold) {
  writeFasta(lseqs[[i]], paste(dir.out, fasta_fold, 
                               paste0(sample_names[[i]], ".fasta"), 
                               sep="/"))
}

#----------------------------------------------------------------------------#
if(is.null(dir.in)) {
  dir.in <- choose.dir(caption="Please, select the directory where the fastq
                       files are located")
}

if(is.null(dir.out)) dir.out <- dir.in
dir.create(dir.out, showWarnings=FALSE, recursive=TRUE)

fns <- list.files(path=dir.in)
fastqs <- fns[grepl(".fastq$", fns)]

#### Filter ####
filtRs <- paste(dir.out,
                 sapply(fastqs,
                        sub,
                        pattern="\\.fastq$",
                        replacement="_filt.fastq.gz"),
                sep="/"
                 )

sample_names <- sub(".fastq$", "", fastqs)


for(i in seq_along(fastqs)) {
  suppressWarnings(fastqFilter(paste(dir.in, fastqs[i], sep="/"), filtRs[i], maxN=0, maxEE=3, 
              truncQ=2, truncLen=bp, compress=TRUE, verbose=TRUE))
}

filtered <- lapply(filtRs, readFastq)
fnSeqs <- unlist(lapply(filtered, length))
retain <- fnSeqs > 0
if (0 %in% fnSeqs) {
  message("NOTE: no reads were retained for one or more samples.
Samples with no reads will be removed from the the filtered set.")
  message("List of sample(s) with zero reads:")
  message(cat(sample_names[!retain], sep="\n"))
}

lsummary <- list()
lsummary[[1]] <- data.frame(Sample= sample_names, nFiltered=fnSeqs)

filtRs <- filtRs[retain]
if(qrep == TRUE) browseURL(report(qa(filtRs)))

#### Dereplicate ####
derepReads <- lapply(filtRs, derepFastq, verbose=TRUE)
names(derepReads) <- sample_names[retain]

unSeqs <- unlist(lapply(derepReads, getnUniques))
lsummary[[2]] <- data.frame(Sample=sample_names[retain], nDerep=unSeqs)
el <- 2
retain <- unSeqs > 1
#### dada ####
if(dada == TRUE) {
  if (1 %in% unSeqs) {
    message("NOTE: Some samples had only one unique sequence.
These samples will be removed before dada analysis.")
    message("List of sample(s) with one unique sequence:")
    message(cat(names(derepReads)[!retain], sep="\n"))
  }
  derepReadsdada <- derepReads[retain]
  dadaReads <- dada(derepReadsdada , err=inflateErr(tperr1,3),
                  errorEstimationFunction=loessErrfun,
                                   selfConsist = TRUE)
  pdf(file = paste(dir.out, "Plot_ErrorRates.pdf", sep="/"))
  for (i in seq_along(dadaReads)) {
    p <- plotErrors(dadaReads[[i]], nominalQ=TRUE)
    print(p + ggtitle(names(dadaReads[i])))
  }
  dev.off()

  nInference <- sapply(dadaReads, ndada)
  el <-el + 1
  lsummary[[el]] <- data.frame(Sample=names(nInference), nInference)

  message("Number of sequenced retained after sample inference with dada2
          (Callahan et al 2015):")
  message(cat(nInference, "\n"))

  lda <- lapply(dadaReads, getUniques)
  nms <- names(derepReads)[!retain]
  ldp <- lapply(nms, dp_extract, derepReads)
  names(ldp) <- nms
  luniseqsFinal <- c(lda, ldp)
}

if(chim == TRUE)  {
  #### Bimera ####
  if(dada == T) {
    bimReads <- sapply(dadaReads, isBimeraDenovo, verbose=TRUE)
  } else {
    bimReads <- sapply(derepReads, isBimeraDenovo, verbose=TRUE)
  }

  message("Proportion of bimeras found in each sample")
  message(cat(round(sapply(bimReads, mean), digits=2), "\n"))

  nChimeras <- sapply(bimReads, nChim)
  message("Number of bimeras found in each sample:")
  message(cat(nChimeras, "\n"))
  el <- el + 1
  lsummary[[el]] <- data.frame(Sample=names(nChimeras), nChimeras)
  if(dada == T) {
    dada_no_chim <- lapply(seq_along(dadaReads), rm.chim, dadaReads, bimReads)
  } else {
    dada_no_chim <- lapply(seq_along(derepReads), rm.chim, derepReads, bimReads)
  }
  names(dada_no_chim) <- names(bimReads)

  #### Reporting ####

  lda <- lapply(dada_no_chim, getUniques)
  nms <- names(derepReads)[!retain]
  ldp <- lapply(nms, dp_extract, derepReads)
  names(ldp) <- nms
  luniseqsFinal <- c(lda, ldp)
} else {
  if(dada == FALSE) {
    nms <- names(derepReads)
    luniseqsFinal <- lapply(nms, dp_extract, derepReads)
    names(luniseqsFinal) <- nms
  }
}

lDNAstr <- lapply(luniseqsFinal, make.DNAString)

nSeq <- sapply(lDNAstr, length)
el <- el + 1
lsummary[[el]] <- data.frame(Sample=names(nSeq), nSeq)
summary <- plyr::join_all(lsummary, by="Sample", type="left")
write.csv(summary, file=paste(dir.out, "Summary.csv", sep="/"), row.names=FALSE)

fasta_fold <- "Final_seqs"
dir.create(paste(dir.out, fasta_fold, sep="/"), showWarnings=FALSE, recursive=TRUE)
lapply(seq_along(lDNAstr), w.fasta, lDNAstr, names(lDNAstr), dir.out, fasta_fold)

return(list(lDNAstr, lsummary))

}
