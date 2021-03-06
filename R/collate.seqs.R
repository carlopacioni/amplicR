


#' Collate sequences
#' 
#' \code{collate.seqs} can be used to collate together outputs from 
#' \code{data.proc}.
#' 
#' This function was developed because, if, for example, sequencing runs are 
#' carried out for the same samples, different  \code{data.proc} outputs will 
#' also exists and you may want to collate together all the sequences from the 
#' same sample. \code{collate.seqs} does exactly this. By passing either the 
#' path to, or a list of, \code{data.proc} outputs, \code{collate.seqs} collate 
#' together the sequences.
#' 
#' \code{collate.seqs} will generate an output that it is very similar to
#' \code{data.proc} but with the combined seqs for the same samples.
#' 
#' @param ldproc If the data are already in R memory, a list whose element are
#'   \code{data.proc} outputs.
#' @param rdas.in If \code{ldproc=NULL}, a character vector with the path to the 
#'   rda files where the outputs from data.proc that have to be collated were 
#'   saved. 
#' @param dir.out The directory where to save the results. If NULL (default) 
#'   then it will be set to the working directory
#' @return Return a list with several elements:
#'   
#'   \itemize{ 
#'   \item $lbysamples: A list whose elements are vectors with the unique 
#'         sequences (names) for each sample and their abundance (values). 
#'   \item $stable: The sequence table 
#'   \item $seq_list: The sequences and matching
#'   sequences IDs 
#'   \item $call: The function call}
#'   
#'   Several files are also saved to disk, these include:
#'   \itemize{ 
#'   \item Seq_table.csv Sequence table
#'   \item Seq_list.csv List of sequences and their matching IDs
#'   \item collated.rda R data file containing the list returned by collate.seqs 
#'      (see above)
#'      }
#' @import data.table
#' @export
#' @examples 
#' # Select the directory where the example data are stored
#' example.data <- system.file("extdata", "HTJ", package="amplicR")
#' # Select a temporary directory where to store the outputs
#' out <- tempdir()
#' 
#' # Run data.proc with bp=140 and save in a sub folder
#' out140 <- paste(out, "out140", sep="/")
#' HTJ.140 <- data.proc(example.data, out140, bp=140)
#' 
#' # Repeat the process in a different sub folder
#' out140bis <- paste(out, "out140bis", sep="/")
#' HTJ.140bis <- data.proc(example.data, out140bis, bp=140)
#'
#' 
#' # Collate seqs of the two objects HTJ.140 and HTJ.140bis 
#' # Note that the sequences are recognised to be identical and size is summed
#' collate.seqs(ldproc=list(HTJ.140, HTJ.140bis), dir.out=out)
#' 
#' # As above but using the .rda files
#' rdas <- paste(out, c("out140", "out140bis"), "data.proc.rda", sep="/")
#' collate.seqs(rdas.in=rdas, dir.out=out)
#' 
#' # Clean up the temp directory
#' unlink(out, recursive=TRUE)


collate.seqs <- function(ldproc=NULL, rdas.in=NULL, dir.out=NULL) {
  if (!requireNamespace("dada2", quietly = TRUE)) {
    stop("Package 'dada2' needed for this function to work. Please install it 
         either manually or using the function amplicR::setup().",
         call. = FALSE)
  }
  #----------------------------------------------------------------------------#
  # Helper functions
  #----------------------------------------------------------------------------#
  make.df <- function(el) {
    df <- data.frame(sequence=names(el), abundance=el, row.names=NULL)
    return(df)
  }
  
  extr.seqs <- function(l) {
    ldf <- lapply(l[[1]], make.df)
    return(ldf)
  }
  
  extr.names <- function(l) {
    snames <- names(l[[1]])
    return(snames)
  }
  
  combine.sample <- function(sample, l) {
    ldf <- l[names(l) == sample]
    dt <- rbindlist(ldf)
    return(dt)
  }
  
  #----------------------------------------------------------------------------#
  call <- sys.call(1)
  if(is.null(ldproc)) {
    if(is.null(rdas.in)) stop("Please, pass a list of data.proc outputs or their path")
    temp.space <- new.env()
    ldproc <- list()
    for (i in seq_along(rdas.in)) {
      dproc <- load(rdas.in[[i]], temp.space)
      ldproc[[i]] <- get(dproc, temp.space)
    }
    rm(temp.space)
  } else {
      if(length(ldproc) < 2) stop("ldproc must be a list of length > 1")
  }
  
  if(is.null(dir.out)) dir.out <- getwd()
  ##################
  lseqs <- list()
  z <- 0
  for (i in seq_along(ldproc)) {
    nsamples <- length(ldproc[[i]][[1]])
    for (j in 1:nsamples) {
      z <- z + 1
      lseqs[[z]] <- make.df(ldproc[[i]][[1]][[j]])
    }
  }
  ##################
  snames_all <- unlist(lapply(ldproc, extr.names))
  names(lseqs) <- snames_all
  
  samples <- unique(snames_all)
  lbysamples <- lapply(samples, combine.sample, lseqs)
  names(lbysamples) <- samples
  
  stable <- suppressWarnings(dada2::makeSequenceTable(lbysamples))
  seqs <- colnames(stable)
  seq_names <- paste0("seq", 1:dim(stable)[2])
  colnames(stable) <- seq_names
  seq_list <- data.frame(seq_names, sequence=seqs)
  
  collate.out <- paste(dir.out, "Collated_seqs", sep="/")
  dir.create(collate.out, recursive=TRUE, showWarnings=FALSE)
  write.csv(stable, file=paste(collate.out, "Seq_table.csv", sep="/"))
  write.csv(seq_list, file=paste(collate.out, "Seq_list.csv", sep="/"), 
            row.names=FALSE)
  
  fasta_fold <- "Final_seqs"
  fasta_dir <- paste(collate.out, fasta_fold, sep="/")
  table2fasta(stable, seq.list=seq_list, dir.out=fasta_dir)
  
  collated <- list(lbysamples=lbysamples, 
                   stable=stable, 
                   seq_list=seq_list, 
                   call=call)
  save(collated, file=paste(collate.out, "collated.rda", sep="/"))
  return(collated)
}


