library(amplicR, quietly=TRUE)
library(adegenet)
context("Test fastsimcoal")

test_that("gl2sfs", {
  
  temp <- tempdir(check = TRUE)
  #temp<-"E:/amplicR/TestFsc"
  suppressMessages(
    amplicR::setup()
  )
  test.data <- system.file("extdata", "fastsimcoal", package="amplicR")
  file.copy(from = test.data, to = temp, overwrite = TRUE, recursive = TRUE)
  dir.in <- file.path(temp, "fastsimcoal")
  load(file.path(dir.in, "glindAmplic.rda"))
  
  sfs <- gl2sfs(glindAmplic, outfile_root ="glindAmplic",  outpath=dir.in)
  res_DAF <- readLines(file.path(dir.in, "Results", "glindAmplic_DAFpop0.obs"))
  res_MAF <- readLines(file.path(dir.in, "Results", "glindAmplic_MAFpop0.obs"))
  
  test_DAF <- readLines(file.path(dir.in, "glindAmplic_DAFpop0.obs"))
  test_MAF <- readLines(file.path(dir.in, "glindAmplic_MAFpop0.obs"))
  
  expect_equal(res_DAF, test_DAF)
  expect_equal(res_MAF, test_MAF)
  unlink(temp, recursive=TRUE)
})

  # fsc.estimate(file.path(dir.in, "Bottleneck"), n=500000, L=100, maf=TRUE, ncpu=4, nBatches=NULL, 
  #              fsc.cmd="fsc2702", fsc.path=dir.in) 
  # AICtable <- AIC.comp(dir.in)
  
test_that("fastsimcoal_pipeline", {
  
  temp <- tempdir(check = TRUE)
  #temp<-"E:/amplicR/TestFsc"
  suppressMessages(
    amplicR::setup()
  )
  test.data <- system.file("extdata", "fastsimcoal", package="amplicR")
  file.copy(from = test.data, to = temp, overwrite = TRUE, recursive = TRUE)
  dir.in <- file.path(temp, "fastsimcoal")
 expect_invisible(
  fsc.multiple.estimate(file.path(dir.in, "Test"), n=500, L=30, fsc.path=dir.in)
  )
  
  AICtable <- AIC.comp(file.path(dir.in, "Test"))
  expect_s3_class(AICtable, "data.frame")
  
  bootstr <- fsc.bootstraps(file.path(dir.in,  "Test", AICtable[AICtable$Rank == 1, "Model"]), 
                            nLoci=10000, nSim=5, 
                            maf=TRUE, ncpu=0, 
                            nBatches=NULL, fsc.cmd="fsc2702", fsc.path=dir.in,
                            par.indLoci="2000000 0", par.nBlocks=1, 
                            par.data="DNA 100 0 2.5e-8 0.33",
                            n=10, L=5,
                            nBoot=100, conf=0.95, boot.type="perc")
  expect_type(bootstr, "list")
  expect_length(bootstr, 4)
  unlink(temp, recursive=TRUE)
})

