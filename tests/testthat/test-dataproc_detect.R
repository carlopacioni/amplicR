library(amplicR)
context("Data processing and detection")

test_that("test data.proc: defaults, no dada, no chim", {
  test.data <- system.file("extdata", "MAPIS900", package="amplicR")
  out <- tempdir()
  
  data(MAPIS900.ref)
  sink(file="NUL")
  suppressWarnings(MAPIS900.test <- data.proc(test.data, out, bp=35))
  MAPIS900.test2 <- data.proc(test.data, out, bp=35, dada=FALSE)
  suppressWarnings(MAPIS900.test3 <- data.proc(test.data, out, bp=35, chim=FALSE))
  sink()
  
  expect_equal(MAPIS900.test , MAPIS900.ref)
  expect_equal(MAPIS900.test2[[2]][[4]][, "nSeq"], c(9, 18))
  expect_equal(MAPIS900.test3[[2]][[4]][, "nSeq"], c(5, 3))
  unlink(out)
})

test_that("test detect: object, read files", {
  test.data <- system.file("extdata", "MAPIS900", package="amplicR")
  out <- tempdir()
  MAPWA <- DNAString("GTGGCACAACCTGTCTGGGCGGGCGTGGACGCCGG")
  
  sink(file="NUL")
  suppressWarnings(MAPIS900.test <- data.proc(test.data, out, bp=35))
  det1 <- detect(MAPIS900.test, dir.out=out, ref_seq=MAPWA)
  det2 <- detect(dir.in=paste(out, "Final_seqs", sep="/"), 
                 dir.out=out, 
                 ref_seq=MAPWA)
  sink()
  
  expect_equal(det1[, "nDiff"], c(19, 19))
  expect_equal(det2[, "nDiff"], c(19, 19))
  unlink(out)
})

