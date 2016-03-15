library(amplicR, quietly=TRUE)
context("Collate sequences")

test_that("test detect: check sum of number of seqs", {
  sink(file="NUL")
 
  example.data <- system.file("extdata", "HTJ", package="amplicR")
  out <- tempdir()
  
  out140 <- paste(out, "out140", sep="/")
  HTJ.140 <- data.proc(example.data, out140, bp=140)
  
  out141 <- paste(out, "out141", sep="/")
  HTJ.141 <- data.proc(example.data, out141, bp=141)
  
  dirs <-c("out140", "out141")
  rdas <- paste(out, dirs, "data.proc.rda", sep="/")
  col_data <- collate.seqs(ldproc=list(HTJ.140, HTJ.141), dir.out=out)

  col_read <- collate.seqs(rdas.in=rdas, dir.out=out)
  checkSum <- sum(col_data$stable[1,])
  sink()
  expect_equal(col_data, col_read)
  expect_equal(checkSum, 666)
  
  # Clean up the temp directory
  unlink(out, recursive=TRUE)
  
})
