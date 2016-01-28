


#' Collate sequences
#' 
#' \code{collate.seqs} can be used to collate together sequences contained in 
#' different fasta files with the same name (where, most likely, the file name 
#' is the sample name where these sequences are coming from).
#' 
#' This function was developed because, if \code{data.proc} is looped over a 
#' number of different \code{bp} values for the same sample files, and the 
#' results stored in several sub folders, you will end up with several fasta 
#' files in each sub folder. Many of these files are generated from the same 
#' samples. However, the sequences contained in each file will be different 
#' because they have been processed using a different \code{bp} values. In such 
#' cases you may want to collate together all the sequences from the same 
#' sample. \code{collate.seqs} does exactly this. By passing the path to 
#' \code{data.proc} outputs, \code{collate.seqs} will search for the sub folder 
#' "Final_seqs", read the fasta files and collate together the ones that have 
#' matching names.
#' 
#' \code{collate.seqs} assumes that the fasta files have extension .fasta.
#' 
#' \code{collate.seqs} assumes that the fasta files were generated with 
#' \code{data.proc}, so it will look for a "Final_seqs" sub folders within the 
#' \code{dirs} path. As long as the folder structure is respected (i.e. the 
#' fasta files are in a sub folder "Final_seqs"), the files from the same 
#' samples have matching names and .fasta is their extension, this function can
#' still be used .
#' 
#' @param dirs A character vector with the path to the folders where the fasta 
#'   files are contained (e.g. where the outputs from data.proc were saved. See 
#'   details)
#' @param dir.out The directory where to save the results. If NULL (default) 
#'   then it will be set to the working directory
#' @return Saves sequences to file in fasta format
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
#' # Run data.proc with bp=140 and save in a sub folder
#' out141 <- paste(out, "out141", sep="/")
#' HTJ.141 <- data.proc(example.data, out141, bp=141)
#' 
#' # Collate seqs of matching file names across the two sub folders 
#' # out140 & out141
#' collate.seqs(dirs=c(out140, out141), out)
#' 
#' # Clean up the temp directory
#' unlink(out, recursive=TRUE)


collate.seqs <- function(dirs, dir.out=NULL) {
  library(ShortRead)
  #----------------------------------------------------------------------------#
  # Helper functions
  #----------------------------------------------------------------------------#
  list.fasta <- function(dir) {
    list.files(paste(dir, "Final_seqs", sep="/"), "*.fasta$")
  }
  
  collate <- function(ffile, dirs, collate.out) {
    for (i in seq_along(dirs)) {
      seqs <- readFasta(paste(dirs[i], "Final_seqs", sep="/"), pattern=ffile)
      writeFasta(seqs, file=paste(collate.out, ffile, sep="/"), mode="a")
    }
  }
  #----------------------------------------------------------------------------#
  if(is.null(dir.out)) dir.out <- getwd()
  ffiles <- lapply(dirs, list.fasta)
  ffiles <- Reduce(intersect, ffiles)
  collate.out <- paste(dir.out, "Collated_seqs", sep="/")
  dir.create(collate.out, recursive=TRUE, showWarnings=FALSE)
  
  lapply(ffiles, collate, dirs, collate.out)
}


