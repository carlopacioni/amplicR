
#' Convert SNPs data from a genlight object (from Dart sequencing) to phased
#' alleles
#'
#'
#' @param gl The genlight object with the processed data
#' @param fastq.dir.in Character vector with the path to the directory where the
#'   fastq files are located. If \code{NULL} and interactive pop up windows will
#'   be used to select the directory. If \code{NA} no sequences are used.
#' @param minLen Minimum length of reads to keep when applying the filter
#' @param min.nSNPs Integer indicating the minimum number of SNPs that a locus
#'   has to have to be retained
#' @param minAbund Either NULL (default) or the minimum number of identical
#'   reads to retain an alleles from the raw sequences
#' @param dir.out Character vector with the name of the directory where to save
#'   the results
#' @param singleAllele Whether only one random allele should be selected for
#'   each sample (TRUE), or both (FALSE)
#' @param nCPUs Integer for the number of CPUs to use (for parallel computation)
#'   or "auto" to automatically call all available CPUs. If 1, no parallel
#'   computation.
#' @inheritParams data.proc
#' @import data.table
#' @import parallel
#' @export
#' 
dart2nexus <- function(gl, fastq.dir.in=NULL, min.nSNPs=3, minAbund=NULL, 
                       minLen=77, truncQ=20, minQ=25,
                       dir.out="Processed_data", singleAllele=TRUE, dada=TRUE, 
                       nCPUs="auto") {
  amplicR::setup()
  
  #### Getting basic info from LocMetrics ####
  proc.data <- data.table(gl$other$loc.metrics)
  proc.data[, lenSeq := nchar(as.character(AlleleSequence))]
  proc.data[, lenTrimSeq := nchar(as.character(TrimmedSequence))]
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
  samplesIDs <- row.names(glm)
  sampleNeeded <- rep(FALSE, length(samplesIDs))
  countLoci <- rep(0, length(samplesIDs))
  names(sampleNeeded) <- samplesIDs
  
  # Check whether samples have more than one heterozygous SNP, which need the raw 
  # sequences to be read to resolve the phase
  for(locus in target.loci) {
    genotypes <- glm[, grep(locus, x=locNamesgl)]
    isHet <- genotypes == 1
    res <- apply(as.matrix(isHet), MARGIN=1, sum, na.rm=TRUE)
    countLoci <- countLoci + (res>1)
    sampleNeeded[res>1] <- TRUE
  }
  
  # These are the samples that will need the sequences 
  countLoci <- countLoci[sampleNeeded]
  
  if(length(countLoci)>0) {
    message("The following samples have this number of loci with multiple 
              possible alleles and sequences would be needed to resolve their phase:\n")
    print(countLoci)
  }
  
  if(is.null(fastq.dir.in)) {
    fastq.dir.in <- choose.dir(caption="Please, select the directory where the fastq
                       files are located")
  } else {
    if(is.na(fastq.dir.in)&length(countLoci>0)) {
      message("No raw sequences were provided but they were needed for at least some samples\n
              For samples with multiple possible alleles, the SNPs will be replace with IUPAC ambiguities")
    }
  }
  #### Filter ####
  if(!is.na(fastq.dir.in)&length(countLoci>0)) {
    
  
  # Read csv files and identify samples needed 
  csvs <- list.files(fastq.dir.in, ".csv$", full.names = TRUE)
  readInfo <- lapply(csvs, fread)
  readInfo <- rbindlist(readInfo)
  keep.these <- readInfo[genotype %in% samplesIDs, targetid]
  
  dir.create(file.path(fastq.dir.in, dir.out), showWarnings=FALSE, recursive=TRUE)
    fastqs <- list.files(path=fastq.dir.in, pattern = ".fastq|FASTQ.{,3}$")
  #fastqs <- fns[grepl(".fastq|FASTQ.{,3}$", fns)]
  if(length(fastqs) == 0) stop(paste("There are no files in", fastq.dir.in,
                                     "with either fastq or fastq.gz extension"))
  fastqs <- fastqs[grep(paste0(paste0("^", keep.these), collapse = "|"), fastqs)]
  if(length(fastqs) == length(keep.these)) 
    message("fastq files were identified for all needed samples") else
      warning(paste(length(keep.these), "fastq files were needed, but", 
                    length(fastqs), "were found"))
  # set up a cluster 
  if(nCPUs != 1) {
    if(nCPUs == "auto") nCPUs <- parallel::detectCores()
    if(length(fastqs)<nCPUs) nCPUs <- length(fastqs)
    cl <- parallel::makeCluster(nCPUs)
    on.exit(expr=parallel::stopCluster(cl))
    catch <- parallel::clusterEvalQ(cl, library("dada2"))
  }
  
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
        dada2::fastqFilter(fn=file.path(fastq.dir.in, fastqs[i]),
                           fout=filtRs[i],
                           maxN=0,
                           minQ=minQ,
                           maxEE=Inf,
                           truncQ=truncQ,
                           minLen=minLen,
                           compress=TRUE,
                           OMP=FALSE,
                           verbose=FALSE)
        #)
    }
    
  } else {
    # mapply(dada2::fastqFilter, fn=file.path(fastq.dir.in, fastqs), fout=filtRs, 
    #        MoreArgs=list(
    #          maxN=0,
    #          maxEE=Inf,
    #          truncQ=0,
    #          minLen=77,
    #          compress=TRUE,
    #          OMP=FALSE,
    #          verbose=FALSE))
    
    clusterExport(cl, varlist=c("fastq.dir.in", "fastqs", "filtRs"), envir=.GlobalEnv) 
      clusterMap(cl=cl, dada2::fastqFilter, fn=file.path(fastq.dir.in, fastqs), fout=filtRs, 
             MoreArgs=list(
               maxN=0,
               minQ=minQ,
               maxEE=Inf,
               truncQ=truncQ,
               minLen=minLen,
               compress=TRUE,
               OMP=FALSE,
               verbose=FALSE),
               .scheduling="dynamic")
    
  }
  )
  message("Done!")
  message(paste("Time needed (in seconds)", round(sys_time[3]), sep = "\n"))
  
  filtRs <- list.files(path=file.path(fastq.dir.in, dir.out, filt_fold), full.names=TRUE)
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
    derepReads <- parLapply(cl, filtRs, fun=dada2::derepFastq, verbose=FALSE)
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
  } # close  if(is.na(fastq.dir.in)&length(countLoci>0))
  
  #### new approach based on genotypes ####
  setkey(readInfo, genotype)
  
  #dir.out <- file.path(dir.out, filt_fold)
  
  #cloneIDs <- sub.proc.data[, CloneID]
  #info_table <- read.csv(info.file)
  #alleles <- sub.proc.data[, TrimmedSequence]
  concatAllele1 <- vector("list", length=length(samplesIDs))
  names(concatAllele1) <- samplesIDs
  concatAllele2 <- vector("list", length=length(samplesIDs))
  names(concatAllele2) <- samplesIDs
  
  setkey(loci, CloneID)
  for(sampleID in samplesIDs) {
message(paste("Processing samples" , sampleID))
    
    allele1 <- Biostrings::DNAStringSet()
    allele2 <- Biostrings::DNAStringSet()
    
    loci_time <- system.time(
    for(locus in target.loci) {
      print(locus)
      baseAllele <- sub.proc.data[J(locus), TrimmedSequence, mult="first"]
      genotypes <- glm[sampleID, grep(locus, x=locNamesgl)]
      SNPpositions <- sub.proc.data[J(locus), SnpPosition, mult="all"]
      lenAllele <- sub.proc.data[J(locus), lenTrimSeq, mult="first"]
      seqAlleles <- make.alleles(baseAllele=baseAllele, 
                                 genotypes=genotypes, 
                                 SNPpositions=SNPpositions, 
                                 lenAllele=lenAllele)
      if(length(seqAlleles) == 1) {
        seqAlleles <- c(seqAlleles, seqAlleles)
      } else {
        if(length(seqAlleles)>2) {
          if(!is.na(fastq.dir.in)&length(countLoci>0)) {
          targets <- readInfo[sampleID, targetid, mult="all"]
          seqs <- vector("list", length(targets))
          for(target in targets) {
          barcode <- readInfo[targetid == target, barcode]
          
          seqs[[which(targets == target)]] <- finReads[[as.character(target)]]$sequence
          nms <- paste(target, seq_along(seqs[[which(targets == target)]]), sep = "_")
          names(seqs[[which(targets == target)]]) <- nms
          
          bar_time <- system.time(
            # Search the sample barcode in the reads
            barcode_hits <- Biostrings::vmatchPattern(
              pattern=Biostrings::DNAString(barcode), 
              subject=seqs[[which(targets == target)]],
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
                               subjectSet=seqs[[which(targets == target)]], 
                               patt_hits=barcode_hits, 
                               patt=barcode, 
                               type="F")
          # remomve barcode
          seqs[[which(targets == target)]] <- Biostrings::DNAStringSet(
            seqs[[which(targets == target)]], start=unlist(starts) + 1)
          retain <- as.logical(S4Vectors::elementNROWS(barcode_hits))
          seqs[[which(targets == target)]] <- seqs[[which(targets == target)]][retain] # keep seqs where the barcode was found
          lenSeqTable <- length(table(width(seqs[[which(targets == target)]]))) # This computes the length in bp of the seqs
          lenSeq <- as.integer(names(table(width(seqs[[which(targets == target)]]))[
            which.max(table(width(seqs[[which(targets == target)]])))]))
          if(lenSeqTable > 1 ) { # if the seqs have variable length
            dir.create(file.path(dir.out, sampleID))
            fn_warning <- file.path(dir.out, sampleID, "inconsistent_length.csv")
            warning(paste("After removing the barcode", barcode, "for sample", 
                          sampleID, "and target", target,
                          "not all the reads have the same length/n", "See", 
                          fn_warning,  "for details"))
            write.csv(data.frame(Seq=names(seqs[[which(targets == target)]])[width(seqs[[which(targets == target)]]) != lenSeq], 
                                 Length=width(seqs[[which(targets == target)]])[width(seqs[[which(targets == target)]]) != lenSeq]),
                      file=fn_warning, row.names=FALSE)
            seqs[[which(targets == target)]] <- seqs[[which(targets == target)]][width(seqs[[which(targets == target)]]) == lenSeq]
            retain <- retain[seq_along(seqs[[which(targets == target)]])]
            warning(paste("All reads of length different from the mode -", lenSeq, 
                          "- were removed."))
          }
          if(lenSeq < sub.proc.data[J(locus), lenTrimSeq, mult="first"]) {
            warning(paste("Warning code: 1. Sample:", sampleID, "Locus:", locus,
                          "target:", target,
                          "Number of possible alleles:", length(seqAlleles), 
                          "No reads of sufficient length. All sequences from this file removed"))
            seqs[[which(targets == target)]] <- DNAStringSet(character(0))
          }
          } # end for target in targets
          # Check if files were dumped
          w <-sapply(seqs, function(x) as.integer(names(table(width(x)))))
          if(class(w) == "list") w <- sapply(w, length)
          
          #### handle here where there are more than two alleles based on genotypes ####
          if(all(w == 0)) { # IF there are no seqs left
            warning(paste("Warning code: 2. Sample:", sampleID, "Locus:", locus,
                          "target: all",
                          "Number of possible alleles:", length(seqAlleles), 
                          "After scanning all target files, no reads of sufficient length found. IUPAC ambiguity codes used"))
            
            hetGen <- which(genotypes == 1)
            seqAlleles <- replaceSNP(genotypes[hetGen], seqAlleles, SNPpositions[hetGen])
            seqAlleles <- c(seqAlleles, seqAlleles)
            next # next locus
          } else {
            seqs <- unlist(do.call(DNAStringSetList, seqs[w>0]))
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
                warning(paste("Warning code: 3. Sample:", sampleID, "Locus:", locus,
                              "Number of possible alleles:", length(seqAlleles), 
                              "Only one sequence found and used"))
                seqAlleles <- seqAlleles[anything>0]
                seqAlleles <- c(seqAlleles, seqAlleles)
              } else {
                newsel.matches <- sel.matches[anything>0] # more than 2 matches
                sourceSeqs <- strsplit(names(temp.seqs[newsel.matches]), "_")
                targets <- sapply(sourceSeqs, '[[', 1)
                
                seqsIndex <- as.integer(
                  sapply(sourceSeqs, '[[', 2)
                  )
                abund <- vector("integer", length=length(targets))
                for(i in seq_along(targets)) {
                  abund <- if(dada == TRUE) {
                    dadaReads[[targets[i]]]$denoised[seqsIndex]
                  } else {
                    derepReads[[targets[i]]]$uniques[seqsIndex]
                  }
                }
                
                warning(paste("Warning code: 4. Sample:", sampleID, "Locus:", locus,
                              "Number of possible alleles:", length(seqAlleles), 
                              "but n suitable sequences:", sum(anything>0),
                              "with respective abundance:", abund))
                seqAlleles <- seqAlleles[anything>0]
                ordAbund <- order(abund)
                seqAlleles <- seqAlleles[ordAbund[1:2]]
              } 
            }
          } else { # Close if(sum(anything>0))
            # if no matches are found
            warning(paste("Warning code: 5. Sample:", sampleID, "Locus:", locus,
                          "Number of possible alleles:", length(seqAlleles), 
                          "No suitable reads found. IUPAC ambiguity codes used"))
            
            hetGen <- which(genotypes == 1)
            seqAlleles <- replaceSNP(genotypes[hetGen], seqAlleles, SNPpositions[hetGen])
            seqAlleles <- c(seqAlleles, seqAlleles)
          }

          } else { # close if seqs were provided
            # if no seqs were provided
            hetGen <- which(genotypes == 1)
            seqAlleles <- replaceSNP(genotypes[hetGen], seqAlleles, SNPpositions[hetGen])
            seqAlleles <- c(seqAlleles, seqAlleles)
          } 
        } # close if(length(seqAlleles)>2)
      }
      seqAlleles <- sample(seqAlleles, 2, replace = FALSE) # Shuffle the alleles so that the base allele is not always the first
      allele1[as.character(locus)] <- DNAStringSet(seqAlleles[1])
      allele2[as.character(locus)] <-  DNAStringSet(seqAlleles[2])
    } # close for(locus)
    )
    loci_time
    concatAllele1[[sampleID]] <- do.call(xscat, allele1) # concatenate alleles
    concatAllele2[[sampleID]] <- do.call(xscat, allele2)
    
  } # close samples
  alnAllele1 <- DNAStringSet(concatAllele1)
  alnAllele2 <- DNAStringSet(concatAllele2)
  names(alnAllele1) <- samplesIDs
  names(alnAllele2) <- samplesIDs
  
  write.nexus(if(singleAllele == FALSE) c(alnAllele1, alnAllele2) else alnAllele1, 
              dir.out=file.path(fastq.dir.in, dir.out), fn="phasedAln.nex", 
              charset=TRUE, aln = TRUE,
              locusIDs=sub.proc.data[, unique(as.character(CloneID))], 
              locusLength=sub.proc.data[J(as.numeric(unique(as.character(CloneID)))), 
                                         lenTrimSeq, mult="first"])
}







