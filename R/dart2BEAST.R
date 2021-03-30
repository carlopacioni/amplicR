library(amplicR)
library(data.table)
oldwd <- getwd()
setwd("C:/Users/Carlo/Dropbox/BEASTly things/Data_handlingTest_Dec2020")

#' @param LocMetrics Fully qualified path to the LocMetrics file
#' @param sampleIDs Character vector with the sample labels to retain
#' @param fastq.dir.in Character vector wtih the path to the directorh where the
#'   fastq files are located
#' @param min.nSNPs Integer indicating the minimum number of SNPs that a locus
#'   has to have to be retained
#' @param dir.out Character vector with the name of the directory where to save
#'   the results
#' @param singleAllele Whether only one random allele should be selected for
#'   each sample (TRUE), or both (FALSE)
#' @param nCPUs Integer for the number of CPUs to use (for parallel computation)
#'   or "auto" to automatically call all available CPUs. If 1, no parallel
#'   computation.
#'   @inheritParams data.proc
#'  @import data.table
#'  @import parallel 
#'   
dart2nexus <- function(LocMetrics, samplesIDs, fastq.dir.in=NULL, min.nSNPs=3, truncQ=20, minQ=25,
                       dir.out="Processed_data", singleAllele=TRUE, dada=TRUE, nCPUs="auto") {
  amplicR::setup()
  
  #### Reading and getting basic info from LocMetrics ####
  proc.data <- fread(LocMetrics)
  proc.data[, lenSeq := nchar(AlleleSequence)]
  proc.data[, lenTrimSeq := nchar(TrimmedSequence)]
  loci <- proc.data[, .N, by=c("CloneID")]
  message(paste(nrow(loci), "loci were detected, with a total of", loci[, sum(N)], "SNPs. "))
  message("The range of SNPs within each locus is ", paste(range(loci[, N]), collapse=" - "))
  message("The frequencies of the number of SNPs within each locus are\n")
  print(table(loci[, N]))
  max.nSNPs <- max(loci[, N]) # max n SNPs within one locus.
  setkey(loci, N)
  target.loci <- loci[J(min.nSNPs:max.nSNPs), CloneID]
  setkey(proc.data, CloneID)
  
  sub.proc.data <- proc.data[J(target.loci), mult="first"]
  setkey(sub.proc.data, CloneID)

  #### Prep to read fastq files ####
  if(is.null(fastq.dir.in)) {
    fastq.dir.in <- choose.dir(caption="Please, select the directory where the fastq
                       files are located")
  }
  
  dir.create(file.path(fastq.dir.in, dir.out), showWarnings=FALSE, recursive=TRUE)
  
  fns <- list.files(path=fastq.dir.in)
  fastqs <- fns[grepl(".fastq|FASTQ.{,3}$", fns)]
  if(length(fastqs) == 0) stop(paste("There are no files in", fastq.dir.in,
                                     "with either fastq or fastq.gz extension"))
  
  #### Read csv files and identify samples needed ####
  csvs <- list.files(fastq.dir.in, ".csv$")
  readInfo <- lapply(csvs, fread)
  readInfo <- rbindlist(readInfo)
  keep.these <- readInfo[genotype %in% samplesIDs, targetid]
  fastqs <- fastqs[grep(paste0(paste0("^", keep.these), collapse = "|"), fastqs)]
  if(length(fastqs) == length(samplesIDs)) 
    message("fastq files were identified for all samples provided") else
      warning(paste(length(samplesIDs), "sample labels were provided, but", 
                    length(fastqs), "were found"))
  # set up a cluster 
  if(nCPUs != 1) {
    if(nCPUs == "auto") nCPUs <- parallel::detectCores()
    if(length(fastqs)<nCPUs) nCPUs <- length(fastqs)
    cl <- parallel::makeCluster(nCPUs)
    on.exit(expr = parallel::stopCluster(cl))
    catch <- parallel::clusterEvalQ(cl, library("dada2"))
  }
    
  
  #### Filter ####
  filt_fold <- "Filtered_seqs"
  dir.create(file.path(fastq.dir.in, dir.out, filt_fold), showWarnings=FALSE, recursive=TRUE)
  filtRs <- paste(fastq.dir.in, dir.out, filt_fold,
                  sapply(fastqs,
                         sub,
                         pattern="\\.fastq.{,3}$|\\.FASTQ.{,3}$",
                         replacement="_filt.fastq.gz"),
                  sep="/"
  )
  message("Applying filter to reads...")
  sys_time <- system.time(
  if(nCPUs == 1) {
    
    for(i in seq_along(fastqs)) {
      
        # suppressWarnings(
        dada2::fastqFilter(fn = file.path(fastq.dir.in, fastqs[i]),
                           fout = filtRs[i],
                           maxN=0,
                           maxEE=Inf,
                           truncQ=truncQ,
                           minLen=77,
                           compress=TRUE,
                           OMP = FALSE,
                           verbose=FALSE)
        #)
    }
    
  } else {
    # mapply(dada2::fastqFilter, fn = file.path(fastq.dir.in, fastqs), fout=filtRs, 
    #        MoreArgs = list(
    #          maxN=0,
    #          maxEE=Inf,
    #          truncQ=0,
    #          minLen=77,
    #          compress=TRUE,
    #          OMP = FALSE,
    #          verbose=FALSE))
    
    clusterExport(cl, varlist=c("fastq.dir.in", "fastqs", "filtRs"), envir=.GlobalEnv) 
      clusterMap(cl = cl, dada2::fastqFilter, fn = file.path(fastq.dir.in, fastqs), fout=filtRs, 
             MoreArgs = list(
               maxN=0,
               minQ=minQ,
               maxEE=Inf,
               truncQ=truncQ,
               minLen=77,
               compress=TRUE,
               OMP = FALSE,
               verbose=FALSE),
               .scheduling="dynamic")
    
  }
  )
  message("Done!")
  message(paste("Time needed (in seconds)", round(sys_time[3]), sep = "\n"))
  
  filtRs <- list.files(path=file.path(dir.out, filt_fold), full.names=TRUE)
  sample_names_fil <- as.integer(sub("_filt.fastq.gz", "", 
                          sapply(filtRs, basename, USE.NAMES=FALSE)))
  
  #### Dereplicate ####
  message("Dereplicating the reads...")
  sys_time <- system.time(
  if(nCPUs == 1) {
  derepReads <- dada2::derepFastq(filtRs, verbose=FALSE)
  } else {
    clusterExport(cl, varlist=c("filtRs"), envir=.GlobalEnv) 
    #derepReads <- lapply(filtRs, dada2::derepFastq, verbose=FALSE)
    derepReads <- parLapply(cl, filtRs, fun = dada2::derepFastq, verbose=FALSE)
  }
  )
  message("Done!")
  message(paste("Time needed (in seconds)", round(sys_time[3]), sep = "\n"))
  names(derepReads) <- sample_names_fil
  
  #lsummary <- list()
  fnSeqs <- unlist(lapply(derepReads, getnFiltered))
  if(length(fnSeqs) == 0) stop(paste("There are no sequences that passed the filter in", 
                                     fastq.dir.in))
  
  #### dada ####
  if(dada == TRUE) {
    message("Applying denoising algorythm...")
    sys_time <- system.time(
      if(nCPUs == 1) {
        dadaReads <- dada2::dada(derepReads, err=dada2::inflateErr(dada2::tperr1,3),
                                 errorEstimationFunction=dada2::loessErrfun, multithread=nCPUs)
      } else {
        clusterExport(cl, varlist=c("derepReads"), envir=.GlobalEnv)
        dadaReads <- parLapply(cl, derepReads, err=dada2::inflateErr(dada2::tperr1,3),
                               errorEstimationFunction=dada2::loessErrfun, multithread=FALSE)
      }
      
      )
    message("Done!")
    message(paste("Time needed (in seconds)", round(sys_time[3]), sep = "\n"))
    
  } else {
    dadaReads <- lapply(derepReads, getSeqFromDerep)
  }
  
  #### Search for alleles and separate by loci ####
  setkey(readInfo, targetid)
  
  #dir.out <- file.path(dir.out, filt_fold)
  
  #cloneIDs <- sub.proc.data[, CloneID]
  #info_table <- read.csv(info.file)
  #alleles <- sub.proc.data[, TrimmedSequence]
  concatAllele1 <- vector("list", length = length(dadaReads))
  names(concatAllele1) <- names(dadaReads)
  concatAllele2 <- vector("list", length = length(dadaReads))
  names(concatAllele2) <- names(dadaReads)
  
  for(i in seq_along(dadaReads)) {
    sampleID <- readInfo[J(sample_names_fil[i]), genotype]
    #dir.create(file.path(dir.out, sampleID))
    barcode <- readInfo[J(sample_names_fil[i]), barcode]
    
    #illqual <- lapply(seq_len(nrow(dadaReads[[i]]$quality)), getRow, 
     #                 m=dadaReads[[i]]$quality)
    #qual <- BStringSet(do.call(c, illqual))
    
    #names(qual) <- nms
    seqs <- dadaReads[[i]]$sequence
    nms <- paste0("seq", seq_along(seqs))
    names(seqs) <- nms
    #ids <- nms
    bar_time <- system.time(
      # Search the sample barcode in the reads
    barcode_hits <- Biostrings::vmatchPattern(
      pattern=Biostrings::DNAString(barcode), 
      subject=seqs,
      max.mismatch=0, 
      min.mismatch=0,
      with.indels=FALSE, fixed=TRUE,
      algorithm="auto")
    )
    
    #if(verbose) {
    # info(patt_hits=ind_hits, fn=fn, table=sub_info_table, row=row, 
    #     gene=gene, name_patt=Find, mismatch=index.mismatch)
    #}
    
    starts <- where2trim(mismatch=0, 
                         subjectSet=seqs, 
                         patt_hits=barcode_hits, 
                         patt=barcode, 
                         type="F")
    # remomve barcode
    seqs <- Biostrings::DNAStringSet(seqs, start=unlist(starts) + 1)
    #qual <- Biostrings::BStringSet(qual, start=unlist(starts) + 1)
    retain <- as.logical(S4Vectors::elementNROWS(barcode_hits))
    seqs <- seqs[retain] # keep seqs where the barcode was found
    #qual <- qual_rm[retain]
    #ids <- ids[retain]
    lenSeqTable <- length(table(width(seqs))) # This computes the length in bp of the seqs
    lenSeq <- as.integer(names(table(width(seqs))[which.max(table(width(seqs)))]))
    if(lenSeqTable > 1 ) { # if the seqs have variable length
      dir.create(file.path(dir.out, sampleID))
      fn_warning <- file.path(dir.out, sampleID, "inconsistent_length.csv")
      warning(paste("After removing the barcode", barcode, "for sample", sampleID,
                    "not all the reads have the same length/n", "See", fn_warning, 
                    "for details"))
      write.csv(data.frame(Seq=names(seqs)[width(seqs) != lenSeq], 
                           Length=width(seqs)[width(seqs) != lenSeq]),
                file = fn_warning, row.names = FALSE)
      seqs <- seqs[width(seqs) == lenSeq]
      retain <- retain[seq_along(seqs)]
      warning(paste("All reads of length different from the mode -", lenSeq, "- were removed."))
    }
    
    setkey(loci, CloneID)
    allele1 <- Biostrings::DNAStringSet()
    allele2 <- Biostrings::DNAStringSet()
    loci_time <- system.time(
    for(locus in target.loci) {
      print(locus)
     if(lenSeq < sub.proc.data[J(locus), lenTrimSeq]) {
        warning(paste("The length of the reads -", lenSeq, 
                      "- is less than the length of the alleles for the locus",
                      locus, "so all alleles for this locus will be Ns "))
        allele1[as.character(locus)] <-  Biostrings::DNAStringSet(paste0(rep("N", 
                                              sub.proc.data[J(locus), lenTrimSeq]), 
                                              collapse = ""))
        allele2[as.character(locus)] <- Biostrings:: DNAStringSet(paste0(rep("N", 
                                              sub.proc.data[J(locus), lenTrimSeq]), 
                                              collapse = ""))
        next
     }
      # Trim seqs for the length of that allele in a temp DNAStringSet
      temp.seqs <- DNAStringSet(seqs, start=1, 
                                end=sub.proc.data[J(locus), lenTrimSeq])
      max.dist <- loci[J(locus), N] # Max n of mismatches based on nSNPs
      
      # compute distance
      system.time(
      sr <- ShortRead::srdistance(temp.seqs, DNAString(sub.proc.data[J(locus), TrimmedSequence]))
      )
      sel.matches <- which(sr[[1]] <= max.dist) # sr is subset [[1]] because 
         # it is a list of length=1 because only had one pattern searched (the target allele)
      sel.matches # This is to speed up testing. Rm once done
      if(length(sel.matches) == 0) {
        warning(paste("No suitable reads found for sample", sampleID, 
                      "and the locus",
                      locus, "so alleles for this locus will be Ns "))
        allele1[as.character(locus)] <- Biostrings::DNAStringSet(paste0(rep("N", 
                                            sub.proc.data[J(locus), lenTrimSeq]), 
                                                            collapse = ""))
        allele2[as.character(locus)] <- Biostrings:: DNAStringSet(paste0(rep("N", 
                                            sub.proc.data[J(locus), lenTrimSeq]), 
                                                            collapse = ""))
        next
      }
      # Check if there are gaps
      aln <- pairwiseAlignment(temp.seqs[sel.matches], DNAString(sub.proc.data[J(locus), TrimmedSequence]))
      noGaps <- grep("-", alignedPattern(aln), invert = TRUE)
      sel.matches <- sel.matches[noGaps] # remove the seqs that align with gaps
      if(length(sel.matches)>2) {
        n.matches <- length(sel.matches)
          temp.dist <- sr[[1]]
          names(temp.dist) <- names(temp.seqs)
          temp.dist <- sort(temp.dist)
          sel.matches <- which(names(temp.seqs) %in% names(temp.dist)[1:2])
          warning(c(paste("In the sample", sampleID, "there were", 
                          n.matches, 
                          "reads that were \npossibly compatible alleles for locus:",
                          locus),
                    #paste(seqs[sel.matches], sep = "\n"), 
                    "\nretaining only the two with the lowest distance from the reference allele\n"))
        } else {
          sel.matches <- which(names(seqs) %in% names(temp.seqs)[sel.matches])
        }
      
      if(length(sel.matches) == 1) sel.matches <- c(sel.matches, sel.matches)
      
      allele1[as.character(locus)] <- DNAStringSet(seqs[sel.matches[1]], start=1, 
                                           end=sub.proc.data[J(locus), lenTrimSeq])
      allele2[as.character(locus)] <-  DNAStringSet(seqs[sel.matches[2]], start=1, 
                                           end=sub.proc.data[J(locus), lenTrimSeq])
      retain[sel.matches] <- FALSE
      seqs <- seqs[retain]
      retain <- rep(TRUE, length(seqs))
    }
    )
    concatAllele1[[i]] <- do.call(xscat, allele1) # concatenate alleles
    concatAllele2[[i]] <- do.call(xscat, allele2)
    
  }
  alnAllele1 <- DNAStringSet(concatAllele1)
  alnAllele2 <- DNAStringSet(concatAllele2)
  names(alnAllele1) <- names(dadaReads)
  names(alnAllele2) <- names(dadaReads)
  
  write.nexus(if(singleAllele == FALSE) c(alnAllele1, alnAllele2) else alnAllele1, 
              dir.out=dir.out, fn="phasedAln.nex", charset=TRUE, 
              locusIDs=sub.proc.data[, as.character(CloneID)], locusLength==sub.proc.data[, lenTrimSeq])
}







