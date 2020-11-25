library(amplicR, quietly=TRUE)
library(ShortRead, quietly=TRUE)
context("Collate sequences")

test_that("test collate.seqs: check sum of number of seqs", {
  sink(file="NUL")
 
  example.data <- system.file("extdata", "HTJ", package="amplicR")
  out <- tempdir()
  
  out140 <- paste(out, "out140", sep="/")
  HTJ.140 <- data.proc(example.data, out140, bp=140)
  
  dir.create(path=paste(out, "exDatMod", sep="/"), 
             showWarnings=FALSE, 
             recursive=TRUE)
  
  fastqList <- list.files(example.data, full.names = TRUE)
  
  file.copy(from = fastqList, 
            to = paste(out, "exDatMod", basename(fastqList), sep="/"))
  
  writeFastq(
    ShortReadQ(sread=DNAStringSet(rep(paste(rep("A", 140), collapse=""), 100)),
               quality=BStringSet(rep(paste(rep("F", 140), collapse=""), 100)),
               id=BStringSet(paste("fake-seq", 1:100, sep="_"))), 
             mode="a", 
             file=paste(out, "exDatMod", "DAFWA#11.fastq.gz", sep="/"))
  
  out140mod <- paste(out, "out140mod", sep="/")
  HTJ.140mod <- data.proc(paste(out, "exDatMod", sep="/"), out140mod, bp=140)
  
  dirs <-c("out140", "out140mod")
  rdas <- paste(out, dirs, "data.proc.rda", sep="/")
  col_data <- collate.seqs(ldproc=list(HTJ.140, HTJ.140mod), dir.out=out)

  col_read <- collate.seqs(rdas.in=rdas, dir.out=out)
  checkSum10 <- sum(col_data$stable[1,])
  checkSum11 <- sum(col_data$stable[2,])
  sink()
  expect_equal(col_data[1:3], col_read[1:3])
  expect_equal(checkSum10, 656)
  expect_equal(checkSum11, 7264)
  
  # Clean up the temp directory
  unlink(out, recursive=TRUE)
  
})
