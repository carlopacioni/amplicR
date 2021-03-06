#' Detect a reference sequence
#' 
#' \code{detect} takes in either the output of \code{data.proc}, or load it up 
#' from a .rda file, and compare the sequences with a reference sequence 
#' reporting the number of mismatch using \code{srdistance} from the package 
#' \code{\link[ShortRead]{ShortRead}}.
#' 
#' The output from \code{data.proc} can be passed with \code{data}. If no data 
#' is passed to \code{detect}, then it will load the .rda file passed with 
#' \code{rda.in}. If \code{rda.in=NULL}, then an interactive window will open to
#' select the location of the file.
#' 
#' If both \code{dir.out=NULL} and \code{rda.in=NULL}, then the path where to 
#' save the results will be asked with an interactive window. If the .rda file 
#' path is provided, then the folder where .rda is located will be selected as
#' output folder.
#' 
#' A summary of the number of sequences found and the minimum number of mismatch
#' within each sample is returned as \code{matrix}, for each reference sequence, 
#' with the same layout as the 
#' sequence table. There will be as many tables as the length of the character 
#' vector passed with \code{ref_seqs}. These results are also written to disk as
#' text files along with the alignments of the sequences provided with the 
#' reference sequence (in the folder "Final_alns"). The alignments are built using
#' \code{PairwiseAlignments} from the package \code{Biostrings}.
#' 
#' Lastly a detect_table is returned (and written to disk) where each row is a 
#'   sequence with the number of mismatch with each reference sequence (columns).
#'   The first column ("nSeq_tot") is the total number of reads for each sequence.
#'   
#' 
#' @param data The output from \code{data.proc}
#' @param rda.in The fully qualified (i.e. including the path) name of the .rda 
#'   file where the output from \code{data.proc} is saved
#' @param ref_seqs A \strong{named} character vector with the reference 
#'   sequence(s)
#' @param dir.out The path where to save the results. If NULL and data is also 
#'   NULL, the directory where the .rda file is located is used. If no file path
#'   is provided, an interactive windows is used to select the folder
#' @return A \code{list} with four elements:
#'   \itemize{
#'     \item $detect_results A list with a result table, for each reference 
#'            sequence, with the minimum mismatch count
#'     \item $alns A list with an alignment, for each reference sequence as 
#'            elements, of the sequences in the sequence table with the 
#'            reference sequence
#'     \item $detect_table A \code{data.table} with the sequence IDs as rows and 
#'            a column with total sequence abundance. All the other columns are 
#'            reference sequences. Values are the minimum number of 
#'            differences with the reference sequence 
#'     \item $call: The function call      
#'     }
#' 
#'   These results are also written to text files
#' @seealso \code{\link{data.proc}}, \code{\link[ShortRead]{srdistance}}, 
#'   \code{\link[Biostrings]{PairwiseAlignments}}
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
#' # Use 'detect' to verify the presence of Mycobacteriumavium subspecies 
#' # paratuberculosis
#' det <- detect(HTJ.test, dir.out=out, ref_seqs=HTJ)
#' 
#' # Clean up the temp directory
#' unlink(out, recursive=TRUE)

detect <- function(data=NULL, rda.in=NULL, dir.out=NULL, ref_seqs) {
  #----------------------------------------------------------------------------#
  if (!requireNamespace("dada2", quietly = TRUE)) {
    stop("Package 'dada2' needed for this function to work. Please install it 
         either manually or using the function amplicR::setup().",
         call. = FALSE)
  }
  #----------------------------------------------------------------------------#
  call <- sys.call(1)
  
  if(is.null(data)) {
    if(is.null(rda.in)) {
      message("Please, select the directory where data.proc output is saved (.rda file)")
      rda.in <- file.choose()
    }
  }
  
  if(is.null(data)) {
    temp.space <- new.env()
    data <- load(rda.in, temp.space)
    data <- get(data, temp.space)
    rm(temp.space)
  } 
  seq_list <- data[[4]]
  stable <- data[[3]]
  DNAstr <- Biostrings::DNAStringSet(seq_list[, "sequence"])
  names(DNAstr) <- seq_list[, "seq_names"]

    if(is.null(dir.out)) {
    if(is.null(rda.in)) {
      dir.out <- choose.dir(caption="Please, select the directory where to save 
                            results") 
    } else {
      dir.out <- basename(rda.in)
    }
  }
      
  #### Detect ####
  lres <- list()
  lalns <- list()
  nSeq_tot <- apply(stable, 2, sum)
  detect_table <- data.table(seq_names=seq_list[, "seq_names"], nSeq_tot)
  aln.out <- paste(dir.out, "Final_alns", sep="/")
  dir.create(aln.out, showWarnings=FALSE, recursive=TRUE)
  
  for (i in seq_along(ref_seqs)) {
    nDiff <- unlist(ShortRead::srdistance(DNAstr, ref_seqs[i]))
    MDiff <- matrix(nDiff, nrow=dim(stable)[1], ncol=length(nDiff), byrow=TRUE)
    stable[stable == 0] <- NA
    stable <- stable > 0
    MDiff <- MDiff * stable
    lres[[names(ref_seqs)[i]]] <- MDiff
    write.csv(MDiff, 
              file=paste(dir.out,  
                         paste0("Results_detect_", names(ref_seqs)[i], ".csv"), 
                         sep="/"))
    aln <- Biostrings::pairwiseAlignment(DNAstr, ref_seqs[i])
    Biostrings::writePairwiseAlignments(aln, 
                            paste(aln.out, 
                            paste0("aln_", names(ref_seqs)[i], ".fasta"), sep="/"), 
                            block.width=2000)
    lalns[[names(ref_seqs)[i]]] <- aln
    detect_table[, (names(ref_seqs)[i]) := unname(nDiff)]
  }
 
  write.csv(detect_table, 
            file=paste(dir.out, "detect_table.csv", sep="/"), row.names=FALSE)
  return(list(detect_results=lres, 
              alns=lalns, 
              detect_table=detect_table,
              call=call))
  }
