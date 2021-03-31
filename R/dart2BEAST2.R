library(amplicR)
library(data.table)
oldwd <- getwd()
setwd("C:/Users/Carlo/Dropbox/BEASTly things/Data_handlingTest_Dec2020")

#' @param gl The ginlight object with the processed data
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
dart2nexus <- function(gl, fastq.dir.in=NULL, min.nSNPs=3, truncQ=20, minQ=25,
                       dir.out="Processed_data", singleAllele=TRUE, dada=TRUE, nCPUs="auto") {
  amplicR::setup()
  
  #### Reading and getting basic info from LocMetrics ####
  proc.data <- gl$other$loc.metrics
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
  
  sub.proc.data <- proc.data[J(target.loci), mult="all"]
  setkey(sub.proc.data, CloneID)

  #### Prep to read fastq files ####
  glm <- as.matrix(gl)
  locNamesgl <- names(glm[1,])
  sampleIDs <- row.names(glm)
  sampleNeeded <- rep(FALSE, length(sampleIDs))
  names(sampleNeeded) <- sampleIDs
  
  for(locus in target.loci) {
    genotypes <- glm[, grep(locus, x = locNamesgl)]
    isHet <- genotypes == 1
    res <- apply(isHet, MARGIN = 1, sum, na.rm=TRUE)
    sampleNeeded[res>1] <- TRUE
  }
  sampleNeeded <- names(sampleNeeded[sampleNeeded])
  if(length(sampleNeeded)>0) {
    warning(c("Raw sequence data are needed for the following samples:\n", paste(sampleNeeded, collapse = "; ")))
  }
  
  
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
  if(!is.null(minAbund)) derepReads <- lapply(derepReads, subsetDerep, minAbund)
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
    finReads <- dadaReads
  } else {
    finReads <- lapply(derepReads, getSeqFromDerep)
  }
  
  #### Search for alleles and separate by loci ####
  setkey(readInfo, targetid)
  
  #dir.out <- file.path(dir.out, filt_fold)
  
  #cloneIDs <- sub.proc.data[, CloneID]
  #info_table <- read.csv(info.file)
  #alleles <- sub.proc.data[, TrimmedSequence]
  concatAllele1 <- vector("list", length = length(finReads))
  names(concatAllele1) <- names(finReads)
  concatAllele2 <- vector("list", length = length(finReads))
  names(concatAllele2) <- names(finReads)
  
  for(i in seq_along(finReads)) {
    sampleID <- readInfo[J(sample_names_fil[i]), genotype]
    #dir.create(file.path(dir.out, sampleID))
    barcode <- readInfo[J(sample_names_fil[i]), barcode]
    
    #illqual <- lapply(seq_len(nrow(finReads[[i]]$quality)), getRow, 
     #                 m=finReads[[i]]$quality)
    #qual <- BStringSet(do.call(c, illqual))
    
    #names(qual) <- nms
    seqs <- finReads[[i]]$sequence
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
    
    #### new approach based on genotypes ####
    IUPAC <- c("AC", "AG", "AT", "CG", "CT", "GT", "CA", "GA", "TA", "GC", "TC", "TG")
    names(IUPAC) <- c("M", "R", "W", "S", "Y", "K", "M", "R", "W", "S", "Y", "K")
    glm <- as.matrix(gl)
    locNamesgl <- names(glm[1,])
    loci_time <- system.time(
    for(locus in target.loci) {
      print(locus)
      baseAllele <- sub.proc.data[J(locus), TrimmedSequence, mult="first"]
      genotypes <- glm[sampleID, grep(locus, x = locNamesgl)]
      # remember that SNP position are one behind because position 1 is indexed as 0
      SNPpositions <- sub.proc.data[J(locus), SnpPosition, mult="all"] 
      breaks <- c(SNPpositions, 
                  # if the last SNP position is at the end of the allele sequence, 
                  # this need to be dropped as it is already being taken care of by the next line
                  if((max(SNPpositions) + 1) == sub.proc.data[J(locus), lenTrimSeq, mult="first"]) {
                    head(SNPpositions, -1) + 1
                  } else {
                    SNPpositions + 1  
                    }
                  , sub.proc.data[J(locus), lenTrimSeq, mult="first"])
      breaks <- sort(breaks)
      nsections <- length(SNPpositions) * 2 + 
        if((max(SNPpositions) + 1) == sub.proc.data[J(locus), lenTrimSeq, mult="first"]) 0 else 1
      sections <- vector("list", length=nsections)
      seqAlleles <- vector("character")
      seqAlleles <- ""
      s <- 1
      #section <- 1
      for(section in (seq_len(nsections))) {
        sections[[section]] <- substr(baseAllele, start = s, stop = breaks[section])
        if(breaks[section] %in% (SNPpositions + 1)) {
          whichGen <- which(SNPpositions + 1 == breaks[section])
          baseSNP <- substr(names(genotypes)[whichGen], 
                 start = nchar(names(genotypes)[whichGen]) - 2, 
                 stop = nchar(names(genotypes)[whichGen]) - 2)
          altSNP <- substr(names(genotypes)[whichGen], 
                           start = nchar(names(genotypes)[whichGen]), 
                           stop = nchar(names(genotypes)[whichGen]))
          if(genotypes[whichGen] == 0) {
            sections[[section]] <- baseSNP
          } else {
            if(genotypes[whichGen] == 2) {
              sections[[section]] <- altSNP
            } else {
              if(is.na(genotypes[whichGen])) {
                sections[[section]] <- names(which(IUPAC == 
                                      paste(c(baseSNP, altSNP), collapse = "")))
            } else {
              sections[[section]] <- c(baseSNP, altSNP)
              }
          }
          }
        }
          seqAlleles <- unlist(lapply(seqAlleles, paste0, sections[[section]]))
          s <- breaks[section] + 1
          #section <- section + 1
      }
      if(length(seqAlleles) == 1) {
        seqAlleles <- c(seqAlleles, seqAlleles)
      } else {
        if(length(seqAlleles)>2) {
          
          #### handle here where there are more than two alleles based on genotypes ####
          if(lenSeq < sub.proc.data[J(locus), lenTrimSeq]) {
            warning(paste("Warning code: 1. Sample:", sampleID, "Locus:", locus,
            "Number of possible alleles:", length(seqAlleles), 
            "No reads of sufficient length. IUPAC ambiguity codes used"))
            
            hetGen <- which(genotypes == 1)
            
            replaceSNP <- function(genotypes, seqAlleles, SNPpositions) {
              seq <- seqAlleles[1]
              for(whichGen in seq_along(genotypes)) {
                baseSNP <- substr(names(genotypes)[whichGen], 
                                  start = nchar(names(genotypes)[whichGen]) - 2, 
                                  stop = nchar(names(genotypes)[whichGen]) - 2)
                altSNP <- substr(names(genotypes)[whichGen], 
                                 start = nchar(names(genotypes)[whichGen]), 
                                 stop = nchar(names(genotypes)[whichGen]))
                substr(seq, start = SNPpositions[whichGen] + 1, stop = SNPpositions[whichGen] + 1) <-
                  names(which(IUPAC == paste(c(baseSNP, altSNP), collapse = "")))
              }
              return(seq)
            }
            
            seqAlleles <- replaceSNP(genotypes[hetGen], seqAlleles, SNPpositions[hetGen])
            seqAlleles <- c(seqAlleles, seqAlleles)
          }
          # Trim seqs for the length of that allele in a temp DNAStringSet
          temp.seqs <- DNAStringSet(seqs, start=1, 
                                    end=sub.proc.data[J(locus), lenTrimSeq, mult="first"])
          
          # compute distance
          system.time(
            sr <- ShortRead::srdistance(temp.seqs, DNAStringSet(seqAlleles))
          )
          sel.matches <- lapply(sr, function(x) which(x == 0)) # sr is a list of 
          # length=length(subject)
          anything <- sapply(sel.matches, length)
          if(sum(anything>0)) {
            if(sum(anything>0) == 2) {
              seqAlleles <- seqAlleles[anything>0]
            } else {
              if(sum(anything>0) == 1) {
                warning(paste("Warning code: 2. Sample:", sampleID, "Locus:", locus,
                              "Number of possible alleles:", length(seqAlleles), 
                              "Only one sequence found and used"))
                seqAlleles <- seqAlleles[anything>0]
                seqAlleles <- c(seqAlleles, seqAlleles)
              } else {
                newsel.matches <- sel.matches[anything>0]
                seqsIndex <- as.integer(
                  substring(names(temp.seqs[newsel.matches]), first = 4)
                  )
                abund <- if(dada == TRUE) {
                  dadaReads[[i]]$denoised[seqsIndex]
                  } else {
                    derepReads[[i]]$uniques[seqsIndex]
                  }
                warning(paste("Warning code: 3. Sample:", sampleID, "Locus:", locus,
                              "Number of possible alleles:", length(seqAlleles), 
                              "but n suitable sequences:", sum(anything>0),
                              "with respective abundance:", abund))
                seqAlleles <- seqAlleles[anything>0]
                ordAbund <- order(abund)
                seqAlleles <- seqAlleles[ordAbund[1:2]]
              } 
            }
          } else {
            warning(paste("Warning code: 4. Sample:", sampleID, "Locus:", locus,
                          "Number of possible alleles:", length(seqAlleles), 
                          "No suitable reads found. IUPAC ambiguity codes used"))
            
            hetGen <- which(genotypes == 1)
            seqAlleles <- replaceSNP(genotypes[hetGen], seqAlleles, SNPpositions[hetGen])
            seqAlleles <- c(seqAlleles, seqAlleles)
          }
        }
      }
      allele1[as.character(locus)] <- DNAStringSet(seqAlleles[1])
      allele2[as.character(locus)] <-  DNAStringSet(seqAlleles[2])
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
              locusIDs=sub.proc.data[, as.character(CloneID), mult="first"], 
              locusLength==sub.proc.data[, lenTrimSeq], mult="first")
}







