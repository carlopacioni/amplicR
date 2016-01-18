library(amplicR, quietly=TRUE)
context("Sequence detection")

test_that("test detect: object, read files", {
  test.data <- system.file("extdata", "HTJ", package="amplicR")
  out <- tempdir()
  HTJ <- "CTGCGCGCCGGCGATGACATCGCAGTCGAGCTGCGCATCCTGACCAGCCGACGTTCCGATCTGGTGGCTGATCGGACCCGGGCGATCGAACCGAATGCGCGCCCAGCTGCTGGAATACTTTCGGCGCTGGAACGCGCCTT"
  names(HTJ)<- "HTJ"
  bp<-140
  
  sink(file="NUL")
  HTJ.test <- data.proc(test.data, out, bp)
  det1 <- detect(HTJ.test, dir.out=out, ref_seqs=HTJ)
  det2 <- detect(dir.in=paste(out, "Final_seqs", sep="/"), 
                 dir.out=out, 
                 ref_seqs=HTJ)
  sink()
  
  expect_equal(det1[, "HTJ"], c(73, 2))
  expect_equal(det2[, "HTJ"], c(73, 2))
  unlink(out)
})
