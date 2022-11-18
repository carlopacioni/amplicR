library(amplicR, quietly=TRUE)
context("Sequence detection")

test_that("test detect: object, read files", {
  #sink(file="NUL")
  seqs <- c("AAATTTT", "GAATTCT")
  names(seqs) <- c("seq1", "seq2")
  x <- Biostrings::DNAStringSet(seqs)
  locusNames <- c("locus1", "locus2")
  locusLen <- c(3, 4)
  tmpDir <- tempdir(check = TRUE)
  write.nexus(x, aln=TRUE, dir.out=tmpDir, fn="aln.nex", charset=TRUE, 
              locusIDs=locusNames, locusLength=locusLen)
  
  #sink()
  
  expect_true(file.exists(file.path(tmpDir, "aln.nex")))
  unlink(tmpDir)
})
