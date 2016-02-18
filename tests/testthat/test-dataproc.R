library(amplicR, quietly=TRUE)
context("Data processing")

test_that("test data.proc: defaults, no dada, no chim", {
  test.data <- system.file("extdata", "HTJ", package="amplicR")
  out <- tempdir()
  bp<-140
  
  data(HTJ.ref)
  sink(file="NUL")
  HTJ.test <- data.proc(test.data, out, bp)
  HTJ.test2 <- data.proc(test.data, out, bp, dada=FALSE)
  HTJ.test3 <- data.proc(test.data, out, bp, chim=FALSE)
  sink()
  
  expect_equal(HTJ.test , HTJ.ref)
  expect_equal(HTJ.test2[[2]][[4]][, "nSeq"], c(122, 374))
  expect_equal(HTJ.test3[[2]][[4]][, "nSeq"], c(8, 1))
  unlink(out)
})



