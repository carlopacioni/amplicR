library(amplicR, quietly=TRUE)
context("Collate sequences")

test_that("test detect: check sum of number of seqs", {
  sink(file="NUL")
  #----------------------------------------------------------------------------#
  # Helper functions
  check.length <- function(ffile, collated.out) {
    seqs <- readFasta(collated.out, pattern=ffile)
    nSeqs <- length(seqs)
    return(nSeqs)
  }
  
  nSingle <- function(sample_name, data1, data2) {
    nSeqs <- length(data1[[1]][sample_name][[1]]) + 
      length(data2[[1]][sample_name][[1]])
    return(nSeqs)
  }
  #----------------------------------------------------------------------------#
  example.data <- system.file("extdata", "HTJ", package="amplicR")
  out <- tempdir()
  
  out140 <- paste(out, "out140", sep="/")
  HTJ.140 <- data.proc(example.data, out140, bp=140)
  
  out141 <- paste(out, "out141", sep="/")
  HTJ.141 <- data.proc(example.data, out141, bp=141)
  
  collate.seqs(dirs=c(out140, out141), out)
  
  list.fasta <- list.files(paste(out, "Collated_seqs", sep="/"), "*.fasta$")
  collated.out <- paste(out, "Collated_seqs", sep="/")
  
  nCollated <- sapply(list.fasta, check.length, collated.out)
  sample_names <- sub(".fasta$", "", list.fasta)
  names(nCollated) <- sample_names
  
  sum_nSeqs_single<-sapply(sample_names, nSingle, HTJ.140, HTJ.141)
  sink()
  expect_equal(nCollated, sum_nSeqs_single)
  
  # Clean up the temp directory
  unlink(out, recursive=TRUE)
  
})
