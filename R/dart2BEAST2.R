
#'Convert SNPs data from a genlight object (from Dart sequencing) to phased
#'alleles
#'
#'This function was developed to convert Dart sequencing data into phased allele
#'sequences. As such, the formatting of the data is expected to be what is
#'provided by the dartR package. That is, the data used as input are generally
#'generated using \code{dartR::gl.read.dart()} and processed (filter) in dartR.
#'Custom data can be used with this function as long as they are formatted in a
#'compatible way. At its mimimum, the last three characters of
#'\code{locNames(gl)} has to have the base and alternative SNP separate by a
#'forward slash '/', and have a loc.metrics element with at least the following
#'headings: \itemize{ \item CloneID \item AlleleSequence \item TrimmedSequence
#'\item SnpPosition }
#'
#'The output of this function is a nexus alignment with the concatenated
#'sequences of the alleles. Once generated, these data can be used for
#'phylogenetic analyses in software such as BEAST (ref).
#'
#'Here, I define 'loci' as one segment (read) from the sequencing data. These
#'are generally identified as CloneID in Dart data. One locus can contain
#'multiple SNPs and often, for phylogenetic analyses, only loci with multiple
#'SNPs are of interest because they contain more information (see Trucchi et al.
#'2014).
#'
#'\code{dart2nexus} initially uses information from the genotypes to create all
#'possible allele sequences. If within any given read, there is no more than one
#'SNP that is heterozygous, there are at the most two possible alleles. If
#'instead, there are >1 heterozygous SNPs, raw sequence data need to be provided
#'to resolve the phase of the allele. If raw sequences are not available, IUPAC
#'ambiguities will be used.
#'
#'Where raw sequence data are available, these are read and processed using a
#'combination of R packages. Most importantly, filtering of the sequences is
#'done using \code{dada::fastqFilter}. Even when the sequecens are provided,
#'there might be situations where multiple reads are possible candidates for the
#'alleles. in such cases also IUPAC ambiguities codes are used. This may happen
#'when spurious reads are retained after filtering. There is not clear way to
#'replicate Dart sequence processing, and playing around with the settings may
#'help in removing these spurious reads. In the limited testing I have
#'conducted, the default settings seems to work adequately in most cases, but
#'there is no guarantee.
#'
#'After filtering the sequences with \code{dada::fastqFilter} it is possible to
#'remove all the (unique) reads that do not have a minimum abundance threshold
#'with \code{minAbund}.
#'
#'It is also possible to include and additional step using dada::dada (see
#'\code{?dada::dada} for more information) after filtering. This is achieved by
#'setting \code{dada=TRUE}. In the testing done, this doesn't seem to help much,
#'while it may help resolving a few loci, it seems to leave many more loci with
#'missing data and it should be consider somewhat experimental. I found that
#'often, using \code{minAbund}, gives the same improvement without leaving other
#'loci with missing data.
#'
#'If sequences are provided, \code{dart2nexus} is also expecting to find in the
#'same location the .csv files the map the dastq files of the sequences with the
#'sample name that are listed in \code{gl$other$loc.metrics}. In the .csv, the
#'fastq file names are generally identified as targetid. The sample ID are
#'typically in a column named genotype. These files also contain teh barcode
#'sequence. The same samples may have sequences in multiple targetid.
#'
#'@param gl The genlight object with the processed data
#'@param dir.in Character vector with the path to the directory where the fastq
#'  and targets.csv files are located. If \code{NULL} and interactive pop up
#'  windows will be used to select the directory. If \code{NA} no sequences are
#'  used.
#'@param minLen Minimum length of reads to keep when applying the filter
#'@param min.nSNPs Integer indicating the minimum number of SNPs that a locus
#'  has to have to be retained
#'@param minAbund Either NULL (default) or the minimum number of identical reads
#'  that an alleles from the raw sequences needs to have to be retained
#'@param dir.out Character vector with the name of the directory where to save
#'  the results
#'@param singleAllele Whether only one random allele should be selected for each
#'  sample (TRUE), or both (FALSE)
#'@param nCPUs Integer for the number of CPUs to use (for parallel computation)
#'  or "auto" to automatically call all available CPUs. If 1, no parallel
#'  computation.
#'@inheritParams data.proc
#'@import data.table
#'@import parallel
#'@export
#'@references Trucchi, E., P. Gratton, J. D. Whittington, R. Cristofari, Y. Le
#'  Maho, N. C. Stenseth and C. Le Bohec, 2014: King penguin demography since
#'  the last glaciation inferred from genome-wide data. Proceedings of the Royal
#'  Society B: Biological Sciences, 281, 20140528.
#'
#'  
dart2nexus <- function(gl, dir.in=NULL, min.nSNPs=3, minAbund=NULL, 
                       minLen=77, truncQ=20, minQ=25,
                       dir.out=NULL, singleAllele=TRUE, dada=TRUE, 
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
  
  if(is.null(dir.in)) {
    dir.in <- choose.dir(caption="Please, select the directory where the fastq
                       files are located")
  } else {
    if(is.na(dir.in)&length(countLoci>0)) {
      message("No raw sequences were provided but they were needed for at least some samples\n
              For samples with multiple possible alleles, the SNPs will be replace with IUPAC ambiguities")
    }
  }
  if(is.null(dir.out)) dir.out <- dir.in

  #### Filter ####
  if(!is.na(dir.in)&length(countLoci>0)) {
    
  
  # Read csv files and identify samples needed 
  csvs <- list.files(dir.in, ".csv$", full.names = TRUE)
  readInfo <- lapply(csvs, fread)
  readInfo <- rbindlist(readInfo)
  keep.these <- readInfo[genotype %in% names(countLoci), targetid]
  dir.out <- file.path(dir.out, "Processed_data")
  dir.create(dir.out, showWarnings=FALSE, recursive=TRUE)
    fastqs <- list.files(path=dir.in, pattern = ".fastq|FASTQ.{,3}$")
  #fastqs <- fns[grepl(".fastq|FASTQ.{,3}$", fns)]
  if(length(fastqs) == 0) stop(paste("There are no files in", dir.in,
                                     "with either fastq or fastq.gz extension"))
  fastqs <- fastqs[grep(paste0(paste0("^", keep.these), collapse = "|"), fastqs)]
  if(length(fastqs) == length(keep.these)) 
    message("fastq files were identified for all needed samples") else
      warning(paste(length(keep.these), "fastq files were needed, but", 
                    length(fastqs), "were found"))
  r <- regexec("\\.fastq|\\.FASTQ.{,3}$", fastqs)
  m <- regmatches(fastqs, r, invert = TRUE)
  targetidsExist <- sapply(m, function(x) x[1])
  sampleExist <- readInfo[targetid %in% targetidsExist, unique(genotype), mult="all"]
  
    # set up a cluster 
  if(nCPUs != 1) {
    if(nCPUs == "auto") nCPUs <- parallel::detectCores()
    if(length(fastqs)<nCPUs) nCPUs <- length(fastqs)
    cl <- parallel::makeCluster(nCPUs)
    on.exit(expr=parallel::stopCluster(cl))
    catch <- parallel::clusterEvalQ(cl, library("dada2"))
  }
  
  filt_fold <- "Filtered_seqs"
  dir.create(file.path(dir.out, filt_fold), showWarnings=FALSE, recursive=TRUE)
  filtRs <- paste(dir.out, filt_fold,
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
        dada2::fastqFilter(fn=file.path(dir.in, fastqs[i]),
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
    # mapply(dada2::fastqFilter, fn=file.path(dir.in, fastqs), fout=filtRs, 
    #        MoreArgs=list(
    #          maxN=0,
    #          maxEE=Inf,
    #          truncQ=0,
    #          minLen=77,
    #          compress=TRUE,
    #          OMP=FALSE,
    #          verbose=FALSE))
    
    clusterExport(cl, varlist=c("dir.in", "fastqs", "filtRs"), envir=environment()) 
      clusterMap(cl=cl, dada2::fastqFilter, fn=file.path(dir.in, fastqs), fout=filtRs, 
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
  
  filtRs <- list.files(path=file.path(dir.out, filt_fold), full.names=TRUE)
  sample_names_fil <- as.integer(sub("_filt.fastq.gz", "", 
                          sapply(filtRs, basename, USE.NAMES=FALSE)))
  
  #### Dereplicate ####
  message("Dereplicating the reads...")
  sys_time <- system.time(
  if(nCPUs == 1) {
  derepReads <- dada2::derepFastq(filtRs, verbose=FALSE)
  } else {
    clusterExport(cl, varlist=c("filtRs"), envir=environment()) 
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
                                     dir.in))
  if(!is.null(minAbund)) derepReads <- lapply(derepReads, subsetDerep, minAbund)
  #### dada ####
  if(dada == TRUE) {
    message("Applying denoising algorythm...")
    sys_time <- system.time(
      if(nCPUs == 1) {
        dadaReads <- dada2::dada(derepReads, err=dada2::inflateErr(dada2::tperr1,3),
                                 errorEstimationFunction=dada2::loessErrfun, multithread=nCPUs)
      } else {
        clusterExport(cl, varlist=c("derepReads"), envir=environment())
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
  } # close  if(is.na(dir.in)&length(countLoci>0))
  
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
          if(!is.na(dir.in)&length(countLoci>0)&sampleID %in% sampleExist) {
          targets <- readInfo[sampleID, targetid, mult="all"]
          seqs <- vector("list", length(targets))
          for(target in targets) {
            if(target %in% targetidsExist) {
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
            } else {
              next
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
            suppressWarnings(
            sr <- ShortRead::srdistance(temp.seqs, DNAStringSet(seqAlleles))
            )
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
    concatAllele1[[sampleID]] <- do.call(xscat, as.list(allele1)) # concatenate alleles
    concatAllele2[[sampleID]] <- do.call(xscat, as.list(allele2))
    
  } # close samples
  alnAllele1 <- DNAStringSet(concatAllele1)
  alnAllele2 <- DNAStringSet(concatAllele2)
  names(alnAllele1) <- samplesIDs
  names(alnAllele2) <- samplesIDs
  
  write.nexus(if(singleAllele == FALSE) c(alnAllele1, alnAllele2) else alnAllele1, 
              dir.out=dir.out, fn="phasedAln.nex", 
              charset=TRUE, aln = TRUE,
              locusIDs=sub.proc.data[, unique(as.character(CloneID))], 
              locusLength=sub.proc.data[J(as.numeric(unique(as.character(CloneID)))), 
                                         lenTrimSeq, mult="first"])
  return(list(Allele1=alnAllele1, Allele2= alnAllele2))
}







