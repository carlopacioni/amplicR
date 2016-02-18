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
  det2 <- detect(rda.in=paste(out, "data.proc.rda", sep="/"), 
                 dir.out=out, 
                 ref_seqs=HTJ)
  sink()
  
  expect_equal(unname(det1[[1]][["HTJ"]][1,]), c(NA, 75, 73, 74, 75, 74, 75, 74, 76))
  expect_equal(unname(det1[[1]][["HTJ"]][2,]), c(2, NA, NA, NA, NA, NA, NA, NA, NA))
  expect_equal(unname(det2[[1]][["HTJ"]][1,]), c(NA, 75, 73, 74, 75, 74, 75, 74, 76))
  expect_equal(unname(det2[[1]][["HTJ"]][2,]), c(2, NA, NA, NA, NA, NA, NA, NA, NA))
  unlink(out)
})
