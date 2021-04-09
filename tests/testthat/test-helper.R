library(amplicR, quietly=TRUE)
context("Test helper")

test_that("Extract rows from matrix", {
  m <- matrix(LETTERS[1:21], ncol = 7, byrow = TRUE)
  oneRow <- getRow(m, 2)
  allRows <- lapply(1:3, getRow, m=m)
  rowsAsVec <- unlist(allRows)
  expect_equal(oneRow, "HIJKLMN")
  expect_length(allRows, 3)
  expect_type(allRows, "list")
  expect_type(rowsAsVec, "character")
  expect_length(rowsAsVec, 3)
})

test_that("Extract rows from matrix", {
  seqs <- c("AAATTTT", "GAATTCT")
 names(seqs) <- c("seq1", "seq2")
  x <- Biostrings::DNAStringSet(seqs)
       locusNames <- c("locus1", "locus2")
       locusLen <- c(3, 4)
       tmpDir <- tempdir()
       tmpFile <- basename(tempfile(tmpdir = tmpDir))
       write.nexus(x, aln=TRUE, dir.out=tmpDir, fn=tmpFile, charset=TRUE, 
                 locusIDs=locusNames, locusLength=locusLen)
       expect_true(file.exists(file.path(tmpDir, tmpFile)))
       rl <- readLines(file.path(tmpDir, tmpFile))
       expect_equal(rl[4], "2")
       expect_equal(rl[7], "seq1")
       expect_equal(rl[22], "charset locus1 = 1 - 3 ;")
       
       write.nexus(x, aln=FALSE, dir.out=tmpDir, fn=tmpFile, charset=TRUE, 
                   locusIDs=locusNames, locusLength=locusLen)
       expect_true(file.exists(file.path(tmpDir, tmpFile)))
})
