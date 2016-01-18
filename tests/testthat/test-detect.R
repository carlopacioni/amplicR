library(amplicR, quietly=TRUE)
context("Sequence detection")

test_that("test detect: object, read files", {
  test.data <- system.file("extdata", "HTJ", package="amplicR")
  out <- tempdir()
  HTJ <- "CTGCGCGCCGGCGATGACATCGCAGTCGAGCTGCGCATCCTGACCAGCCGACGTTCCGATCTGGTGGCTGATCGGACCCGGGCGATCGAACCGAATGCGCGCCCAGCTGCTGGAATACTTTCGGCGCTGGAACGCGCCTT"
  bp<-140
  
  sink(file="NUL")
  HTJ.test <- data.proc(test.data, out, bp)
  det1 <- detect(HTJ.test, dir.out=out, ref_seq=HTJ)
  det2 <- detect(dir.in=paste(out, "Final_seqs", sep="/"), 
                 dir.out=out, 
                 ref_seq=HTJ)
  sink()
  
  expect_equal(det1[, "nDiff"], c(73, 2))
  expect_equal(det2[, "nDiff"], c(73, 2))
  unlink(out)
})
