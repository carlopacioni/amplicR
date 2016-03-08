#' Detect a reference sequence
#' 
#' \code{detect} takes in either the output of \code{data.proc}, or load it up 
#' from a .rda file, and compare the sequences with a reference sequence 
#' reporting the number of mismatch using \code{srdistance} from the package 
#' \code{shortRead}.
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
#' within each sample is returned as \code{matrix} with the same layout as the 
#' sequence table and a \code{data.frame} with the sequence IDs and the 
#' number of mismatch. There will be as many tables as the length of the character 
#' vector passed with \code{ref_seqs}. These results are also written to disk as
#' text files along with the alignments of the sequences provided with the 
#' reference sequence (in the folder "Final_alns"). The alignments are built using
#' \code{PairwiseAlignments} from the package \code{Biostrings}.
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
#'     \item $lseq_list A list where elements (one for each reference 
#'            sequence) are a \code{data.frame} with the sequence IDs and the 
#'            minimum number of differences with the reference sequence 
#'     \item $call: The function call      
#'     }
#' 
#'   These results are also written to text files
#' @seealso \code{\link{data.proc}}, \code{\link[shortRead]{srdistance}}, 
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
  library(dada2)
  library(ShortRead)
 
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
  DNAstr <- DNAStringSet(seq_list[, "sequence"])
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
  lseq_list <- list()
  aln.out <- paste(dir.out, "Final_alns", sep="/")
  dir.create(aln.out, showWarnings=FALSE, recursive=TRUE)
  
  for (i in seq_along(ref_seqs)) {
    nDiff <- unlist(srdistance(DNAstr, ref_seqs[i]))
    MDiff <- matrix(nDiff, nrow=dim(stable)[1], ncol=length(nDiff), byrow=TRUE)
    stable[stable == 0] <- NA
    stable <- stable > 0
    MDiff <- MDiff * stable
    lres[[names(ref_seqs)[i]]] <- MDiff
    write.csv(MDiff, 
              file=paste(dir.out,  
                         paste0("Results_detect_", names(ref_seqs)[i], ".csv"), 
                         sep="/"))
    aln <- pairwiseAlignment(DNAstr, ref_seqs[i])
    writePairwiseAlignments(aln, 
                            paste(aln.out, 
                            paste0("aln_", names(ref_seqs)[i], ".fasta"), sep="/"), 
                            block.width=2000)
    lalns[[names(ref_seqs)[i]]] <- aln
    detect.seq_list <- data.frame(seq_names=seq_list[, "seq_names"], 
                                  nDiff=unname(nDiff))
    lseq_list[[names(ref_seqs)[i]]] <- detect.seq_list
    write.csv(detect.seq_list, 
              file=paste(dir.out,  
                         paste0("detect.seq_list_", names(ref_seqs)[i], ".csv"), 
                         sep="/"),
              row.names=FALSE)
  }
  return(list(detect_results=lres, 
              alns=lalns, 
              detect_seq_list=lseq_list, 
              call=call))
  }
