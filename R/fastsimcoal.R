#' Convert a genlight object into a sfs input file
#' 
#' The output of this function is suitable for analysis in fastsimcoal2 or dada.
#' 
#' It saves a derived sfs, assuming that the reference allele is the ancestral, 
#' and a MAF sfs.
#' 
#' At this stage this function caters only for diploid organisms, 
#' for samples from one population only, and for genotypes without missing data.
#' 
#' It expects a dartR formatted genlight object, but will work with other gl as 
#' long as the have a \code{data.frame} in \cod{other} called loc.metrics with 
#' a column named \code{TrimmedSequence} that contains a character string whose 
#' length corrensponds to the number of sites of the allele. This information is 
#' used to obtained the number of non-polymorphic sites.
#' 
#' TO-DO: deal with missing data. See http://evomics.org/learning/population-and-speciation-genomics/2018-population-and-speciation-genomics/fastsimcoal2-activity/
#' 
#' @param gl the genlight object
#' @param outfile_root The root of the name of the output file
#' @param outpath The output path where to save the file
#' @import dartR
#' @export 
gl2sfs <- function(gl, outfile_root="gl2sfs", outpath=tempdir()) {
#---------- Helper ---------------#
count.freq.ref <- function (x) {
  if(x == 0) return(2) else if(x == 1) return(1) else return(0)
}           
#---------------------------------#

if(!all(gl$ploidy == 2)) 
  stop("This function currently caters only for diploid organisms")
if(length(unique(gl$pop))) 
  stop("This function currently caters only for samples from one population only")
if(sum(is.na(glm))>0) 
  stop("This function  currently caters only for genotypes without missing data")
glm <- as.matrix(gl)
daf <- table(colSums(glm))

names(daf) <- paste("d0", seq(0, length(daf) - 1), sep="_")
inv <- sum(nchar(as.character(gl$other$loc.metrics$TrimmedSequence))) - gl@n.loc
daf[1] <- daf[1] + inv
writeLines(c(paste(paste("d0", seq(0, length(daf) - 1), sep="_"), collapse=" "),
           paste(as.character(daf), collapse=" ")), 
           con=file.path(outpath, outfile_root, "_DAFpop0.obs"))


freq.ref <- colSums(apply(glm, FUN=count.freq.ref, MARGIN=c(1,2)))
freq.alt <- colSums(glm)
maf <- table(c(freq.ref, freq.alt))
maf <- maf[which(as.numeric(names(maf)) <= nrow(glm))]
maf[1] <- maf[1] + inv
if(length(maf) == nrow(glm)) maf[length(maf)] <- maf[length(maf)]/2
writeLines(c(paste(paste("d0", seq(0, length(maf) - 1), sep="_"), collapse=" "),
             paste(as.character(maf), collapse=" ")), 
           con=file.path(outpath, outfile_root, "_MAFpop0.obs"))

return(list(DAF=daf, MAF=maf))
}

#' Run fastsimcoal to estimate parameters' values
#'
#' Run fastsimcoal assuming that in \code{dir.in} there are the tpl, est and obs
#' files properly named and formatted.It writes a fsc_run.txt file and execute
#' fsc in \code{dir.in}
#'
#' If \code{ncpu=0} and \code{nBatches=NULL}, \code{nBatches} is set to 12 (fsc
#' default value). If \code{ncpu}>0 and \code{nBatches=NULL}, \code{nBatches} is
#' set to twice ]code{nspu}.
#'
#'
#' @param dir.in The path where to run fsc
#' @param n The number of coalescent simulations to approximate the expected SFS
#'   (-n option). This should be larger than 100,000.
#' @param L The number of optimization cycles (-L option). It should be >50
#' @param maf Whether a MAF SFS (default) or a derived SFS is provided (if
#'   \code{FALSE})
#' @param ncpu The number of CPU (threads) to use in the analysis. Automatically
#'   handle if \code{ncpu=0} (default)
#' @param nBatches The number of batches (-B option)
#' @fsc.cmd The command to use to call fsc (that may be different depending on
#' the version installed)
#' @param fsc.path The path where fsc is installed or \code{"path"} if it is in
#'   the PATH (that means that it can be called regardless of the working
#'   directory)
#' @export
fsc.estimate <- function(dir.in, n=500000, L=100, maf=TRUE, ncpu=0, nBatches=NULL, 
                         fsc.cmd="fsc2702", fsc.path="path") {
  if(!is.integer(as.integer(ncpu))) stop("The number of CPUs has to be an integer")
  if(ncpu == 0 & is.null(nBatches)) nBatches <- 12
  if(ncpu > 0 & is.null(nBatches)) nBatches <- 2 * ncpu
  tpl <- list.files(dir.in, pattern="*.tpl$", full.names=FALSE)
  est <- list.files(dir.in, pattern="*.est$", full.names=FALSE) 
  cmd <- paste("-t", tpl, "-n", as.integer(n), if(maf == TRUE) "-m" else "-d", "-e", est, 
               "-M", "-L", as.integer(L), "-q", "-c", ncpu, "-B", nBatches, collapse=" ")
  writeLines(c(dir.in, cmd), con=file.path(dir.in, "fsc_run.txt"))
  old.wd <- getwd()
  on.exit(setwd(old.wd))
  setwd(dir.in)
  system2(if(fsc.path == "path") fsc.cmd else file.path(fsc.path, fsc.cmd))
}


#' Run fsc over multiple models
#'
#' This function expects one or more folders within \code{dir.in}, where in each
#' directory there are all the files needed to run \code{fsc.estimate}, which is run iteratively 
#' @inheritParams fsc.estimate
#' @export
fsc.multiple.estimate <- function(dir.in, n=500000, L=100, maf=TRUE, ncpu=0, 
                                  nBatches=NULL, fsc.cmd="fsc2702", fsc.path="path") {
  ld <- list.dirs(dir.in,full.names=TRUE, recursive=FALSE)
  for(d in ld) {
    fsc.estimate(dir.in=d, n=n, L=L, maf=maf, ncpu=ncpu, nBatches=nBatches, 
                 fsc.cmd=fsc.cmd, fsc.path=fsc.path)
  }
  
}


#' AIC
#' 
#' compute the AIC given the likelihood in log10 scale (default fsc output) 
#' and number of parameters, which is automatically read by the function if 
#' \code{k=NULL).
#' 
#' It is important to note that this function will probably fail if there are 
#' more than one .est within each directory as it will read these to evaluate 
#' the number of parameters in the model. If there is more than one file, there 
#' is no way to know which one is correct.
#' 
#' @param dir.in The directory where the analysis was conducted
#' @param k The number of parameters in the model. If \code{NULL} is read 
#' automatically by the function
#' @export
AIC <- function(dir.in, k=NULL) {
  log10toln<-function(l10) {
    rlns=l10/log10(exp(1))
  }
  bLhood <- list.files(dir.in, pattern="*.bestlhoods$", recursive=TRUE, full.names=TRUE)
  path.length <- sapply(bLhood, nchar)
  bLhood <- bLhood[which.min(path.length)]
  lbest <- read.table(bLhood, h=TRUE)[1, "MaxEstLhood"]
  if(is.null(k)) {
  est <- list.files(dir.in, pattern="*.est$", recursive=TRUE, full.names=TRUE)
  est <- est[which.min(sapply(est, nchar))]
  suppressWarnings(
  rlns <- readLines(est)
  )
  start <- max(
    grep("\\[PARAMETERS|//#isInt|//all Ns are in number of haploid individuals", rlns)
  ) + 1
  end <- min(grep("COMPLEX PARAMETERS|RULES", rlns)) - 1
  k <- length(grep("^[0-1]", rlns[start:end]))
  }
  ln=log10toln(lbest)
  return(2*k-2*ln)
}

#' Compute and compare AIC for different models
#' 
#' This functions expects that different models were run in different directories 
#' within the same parent directory (as for \code{fsc.multiple.estimate}). Assuming
#' that \code{dir.in} is the parent directory, it will extract and compute AIC 
#' from all of them.
#' 
#' It is important to note that this function will probably fail if there are 
#' more than one .est within each directory as it will read these to evaluate 
#' the number of parameters in the model. If there is more than one file, there 
#' is no way to know which one is correct. Also, that will also create an 
#' inconsistency between the number of models and the number of files, which may 
#' cause problems in reporting the results even if the function completes. 
#' 
#' @param dir.in The parent directory within which the models have been run
#' @export
AIC.comp <- function(dir.in) {
  mod.nms <- list.dirs(dir.in, full.names=FALSE, recursive=FALSE)
  ld <- list.dirs(dir.in,full.names=TRUE, recursive=FALSE)
  lAICs <- sapply(ld, amplicR::AIC, USE.NAMES=FALSE)
  mod.rank <- rank(lAICs)
  return(data.frame(Model=mod.nms, AIC=lAICs, Delta=lAICs - min(lAICs), Rank=mod.rank))
}

#'Compute bootstrap confidence intervals for estimated parameters with fsc
#'
#'Once fsc estimated parameter values, this function uses the *_maxL.par to
#'simulate data with the model parameters. From there, it estimates the
#'parameter values and confidence intervals are finally returned using the R
#'package \code{boot},
#'
#'It is important that enough sites are simulated to ensure that sufficient
#'polymorphic loci are present in the simulated data. It is better to simulate
#'an excess of sites and retained those needed using \code{nLoci}.
#'
#'Initial values are passed using the .pv so that a reduced numner of replicates can be passed
#' @param dir.in The directory where the analysis was conducted
#'@param nLoci The number of polymorphic loci to retain
#'@param nSim The number of datasets that need to be simulated
#'@param par.indLoci The two integers value that need to be used for the number
#'  of independent loci in the .par that it is used to run the simulations.
#' @param par.nBlocks The number of linkage blocks in the .par that it is used 
#' to run the simulations.
#' @param par.data The string to be used in the 'per Block: data type, num loci, 
#' rec. rate and mut rate + optional parameters' line of the .par that it is used 
#' to run the simulations. The default value is 
#' \code{par.data="DNA 100 0 2.5e-8 0.33"}
#' @param nBoot The number of bootstrap replicates
#' @param conf The confidence level to compute the confidence intervals
#' @param boot.type The method to be used to compute the confidence intervals. 
#'   See \code{?boot.ci} for details. By default, percentile intervals are computed.
#' @inheritParams fsc.estimate
#' @import boot
#' @import data.table
#' @export
fsc.bootstraps <- function(dir.in, nLoci=10000, nSim=100, maf=TRUE, ncpu=0, 
                           nBatches=NULL, fsc.cmd="fsc2702", fsc.path="path",
                           par.indLoci="200000 0", par.nBlocks=1, 
                           par.data="DNA 100 0 2.5e-8 0.33",
                           n=100000, L=50, nBoot=1000, conf=0.95, boot.type="perc") {
  #---------- Helper ---------------#
  med.i <- function(x, i) median(x[i])
  
  extractCI <- function(x, nBoot, conf, boot.type) {
    bootstr <- boot(x, statistic=med.i, R=nBoot)
    b.ci <- boot.ci(bootstr, conf=conf, type=boot.type)
    t0 <- b.ci$t0
    ci <- b.ci[[4]][1, tail(seq_len(ncol(b.ci[[4]])), 2)]
    return(c(t0, ci))
  }
  #----------------------------------------------#
  
  if(!is.integer(as.integer(ncpu))) stop("The number of CPUs has to be an integer")
  if(ncpu == 0 & is.null(nBatches)) nBatches <- 12
  if(ncpu > 0 & is.null(nBatches)) nBatches <- 2 * ncpu
tpl <- list.files(dir.in, pattern="*.tpl$", full.names=TRUE)
est <- list.files(dir.in, pattern="*.est$", full.names=TRUE)
res.dir <- gsub(".tpl$", replacement = "", x = basename(tpl))
par.file <- list.files(file.path(dir.in, res.dir), pattern="*_maxL.par$", recursive=FALSE, full.names=TRUE)
dir.create(file.path(dirname(par.file), "Bootstraps"), showWarnings=FALSE)
new.par.nm <- strsplit(basename(par.file), "_")[[1]][1]
new.par.full <- file.path(dirname(par.file), "Bootstraps", paste0(new.par.nm, "_boot.par"))
file.copy(from=par.file, to=new.par.full, overwrite=TRUE)
# Modify the par file
rlns <- readLines(new.par.full)
par.indLoci.pos <- grep("independent", rlns)
rlns[par.indLoci.pos + 1] <- par.indLoci
par.nBlocks.pos <- grep("linkage blocks$", rlns)
rlns[par.nBlocks.pos + 1] <- par.nBlocks
par.data.pos <- grep("data type", rlns)
rlns[par.data.pos + 1] <- par.data
writeLines(rlns, new.par.full)
# generate cmd and run fsc
cmd <- paste("-i", basename(new.par.full), "-n", as.integer(nSim), if(maf == TRUE) "-m" else "-d", 
             "-j", "-q", "-s", as.integer(nLoci), "-x", "-c", ncpu, "-B", nBatches, collapse=" ")
writeLines(c(dirname(new.par.full), cmd), con=file.path(dirname(new.par.full), "fsc_run.txt"))
old.wd <- getwd()
on.exit(setwd(old.wd))
setwd(dirname(new.par.full))
system2(if(fsc.path == "path") fsc.cmd else file.path(fsc.path, fsc.cmd))
d <- list.dirs(dirname(new.par.full), recursive=FALSE)
ld <- list.dirs(d, recursive=FALSE)


new.tpl.nm <- sub(".tpl", "_boot.tpl", x=basename(tpl))
new.est.nm <- sub(".est", "_boot.est", x=basename(est))
pv <- list.files(dirname(par.file), pattern="*.pv$", full.names=TRUE, recursive=FALSE)
for(dsim in ld) {
  file.copy(from=tpl, to=file.path(dsim, new.tpl.nm), overwrite=TRUE)
  file.copy(from=est, to=file.path(dsim, new.est.nm), overwrite=TRUE)
  cmd <- paste("-t", new.tpl.nm, "-n", as.integer(n), if(maf == TRUE) "-m" else "-d", "-e", 
               new.est.nm, "-M", "-L", as.integer(L), "--initialValues", pv, "-q", "-c", ncpu, 
               "-B", nBatches, collapse=" ")
  writeLines(c(dsim, cmd), con=file.path(dsim, "fsc_run.txt"))
  
  setwd(dsim)
  system2(if(fsc.path == "path") fsc.cmd else file.path(fsc.path, fsc.cmd))
}
bLhood <- list.files(d, pattern="*.bestlhoods$", recursive=TRUE, full.names=TRUE)
lbLhoods <- lapply(bLhood, read.table, header=TRUE)
est.sim <- do.call(rbind, args=lbLhoods)

est.dt <- data.table(est.sim)
bootstr.ci <- est.dt[, lapply(.SD, extractCI, nBoot, conf, boot.type)]
nms <- names(bootstr.ci)
bootstr.ci[, Param:=c("Median", "Lower", "Upper")]
setcolorder(bootstr.ci, "Param")
return(bootstr.ci)
}

