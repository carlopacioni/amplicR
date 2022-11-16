#' Replace null values with zero or x with x - 1
#' 
#' Replace null values with zero or x with x - 1
#' 
#' This function is used internally by amplicR
#' 
#' @param x the element to evaluate
#' @return An obj of same class and length as the input
#' @export
putzero <- function(x){
  if(is.null(x)) x <- 0 else x <- x - 1
  return(x)
}

#' Replace null values with NA
#' 
#' Replace null values with NA
#' 
#' This function is used internally by amplicR
#' 
#' @param x the element to evaluate
#' @return An obj of same class and length as the input
#' @export
putna <- function(x){
  if(is.null(x)) {
    x <- NA 
  } else {
    if(length(x) == 0) x <- NA
  }
  return(x)
}

#' Select position where a pattern has been found
#'
#' Takes a set of sequences and the output from
#' \code{Biostrings::vmatchPattern()} and identify the first positions after the
#' match with a pattern
#' 
#' This function is used internally by amplicR
#' 
#' @param mismatch The number of allowed mismatches
#' @param subjectSet The set of sequences that were used for the search
#' @param patt_hits The output from Biostrings::vmatchPattern()
#' @param patt The patter used for the search
#' @param type The type of match, either forward "F", or reverse "R"
#' @return A list with start positions
#' @export
where2trim <- function(mismatch, subjectSet, patt_hits, patt, type="F") {
  if(length(patt_hits) == 0) return(list())
  
  if(mismatch == 0) {
    if(type == "F") {
      starts <- Biostrings::endIndex(patt_hits)
      starts <- lapply(starts, "[", 1)
      starts <- lapply(starts, putna)
    } else {
      starts <- Biostrings::startIndex(patt_hits)
      starts <- lapply(starts, putna)
      starts <- lapply(starts, max)
    }
    
  } else {
    
    sel_hits <- lapply(1:length(subjectSet), sel.match, 
                       subjectSet=subjectSet, 
                       patt_hits=patt_hits, 
                       patt=patt,
                       type=type)
    
    if(type == "F") {
      starts <- lapply(sel_hits, BiocGenerics::end)
    } else {
      starts <- lapply(sel_hits, BiocGenerics::start)
    }
    
    starts <- lapply(starts, putna)
  }
  
  return(starts) 
}

#' Select appropriate hits when searching a patter
#' 
#' This function is used internally by amplicR
#' 
#' @param i The index to subset patt_hitts
#' @inheritParams where2trim
#' @export
sel.match <- function(i, subjectSet, patt_hits, patt, type="F") {
  if (!requireNamespace("ShortRead", quietly = TRUE)) {
    stop("Package 'ShortRead' needed for this function to work. Please install it 
         either manually or using the function amplicR::setup().",
         call. = FALSE)
  }
  if(length(patt_hits[[i]]) > 1) {
    if(type == "F") {
      patt <- Biostrings::DNAString(patt)
    } else {
      patt <- Biostrings::reverseComplement(Biostrings::DNAString(patt))
    }
    
    v <- IRanges::Views(subjectSet[[i]], patt_hits[[i]])
    hits <- suppressWarnings(as.character(v, use.names=FALSE, check.limits=FALSE))
    sr <- ShortRead::srdistance(Biostrings::DNAStringSet(hits), patt)
    nDiff2 <- unlist(sr)
    nDiff <- unlist(sr, use.names=FALSE)
    if(type == "F") {
      m <- patt_hits[[i]][which.min(nDiff)]
    } else {
      if(length(nDiff) == 0) {
        m <- patt_hits[[i]][integer()]
      } else {
        m <- patt_hits[[i]][max(which(nDiff == min(nDiff)))]
      }
    }
  } else {
    return(patt_hits[[i]])
  }
  
  return(m)
}


#' Return the number of filtered sequences from a derep obj
#' 
#' Return the number of filtered sequences from a derep obj
#' 
#' This function is used internally by amplicR
#' 
#' @param derep 
#' @return A vector
#' @export
getnFiltered <- function (derep) {
  nFiltered <- length(derep$map)
  return(nFiltered)
}

#' Return the number of unique sequences in a derep obj
#' 
#' Return the number of unique sequences in a derep obj
#' 
#' This function is used internally by amplicR
#' 
#' @param derep 
#' @return A vector
#' @export
getnUniques <- function (derep) {
  nUniques <- length(derep$uniques)
  return(nUniques)
}

#' Return the number of unique sequences in a dada obj
#' 
#' Return the number of unique sequences in a dada obj
#' 
#' This function is used internally by amplicR
#' 
#' @param dada_el the element of a dada object
#' @return A vector
#' @export
ndada <- function(dada_el) {
  nSeqs <- length(dada2::getUniques(dada_el))
  return(nSeqs)
}

#' Extracts rows from a matrix into a vector
#' 
#' The row rn becomes an element of a vector. This function is used internally by amplicR
#' @param m A matrix
#' @param rn The row number to extract
#' @return A vector of same type as the matrix
#' @export
#' 
getRow <- function(m, rn) {
  r <- as.integer(round(m[rn,], 0))
  qual <- Biostrings::quality(NumericQuality(r))
  qual <- IlluminaQuality(qual)
  return(qual)
}

#' Write a nexus file from a DNAStringSet object
#'
#' Takes a \code{character} vector or a \code{DNAStringSet} and write a nexus
#' file to disk.
#'
#' If \code{x} is a \code{character} vector, it will convert it to a
#' \code{DNAStringSet}. if \code{aln=FALSE}, \code{write.nexus} will align the
#' sequences first and then write the nexus file. If \code{aln=TRUE}, it will
#' assume that all sequences have the same length and are aligned.
#'
#' If \code{charset=TRUE}, \code{write.nexus} will append a character block at
#' the end of the file. This is sometimes used to partition the alignment when
#' it is imported in software like BEAST2
#'
#' @param x Either a \code{character} vector or a \code{DNAStringSet}
#' @param aln Logical. Whether x is an alignment
#' @param gapOening Integer (negative). Penalty for opening a gap in the
#'   alignment
#' @param gapExtension Integer (negative). Penalty for extending a gap in the
#'   alignment
#' @param dir.out The path where to save the output
#' @param fn The output file name
#' @param charset Whether a charset block should be added at the end of the
#'   nexus file
#' @param locusIDs \code{character} vector with the names of the loci (if
#'   \code{charset=TRUE}) in the same order as they appear in the sequences
#' @param locusLength Integer vector with the length of each locus (if
#'   \code{charset=TRUE}) in the same order as they appear in the sequences
#' @return Writes a nexus file to disk
#' @export
#' @examples seqs <- c("AAATTTT", "GAATTCT")
#'         names(seqs) <- c("seq1", "seq2")
#'         x <- Biostrings::DNAStringSet(seqs)
#'         locusNames <- c("locus1", "locus2")
#'         locusLen <- c(3, 4)
#'         tmpDir <- tempdir()
#'         tmpFile <- tempfile(tmpdir = tmpDir)
#'         write.nexus(x, aln=TRUE, dir.out=tmpDir, fn=basename(tmpFile), charset=TRUE, 
#'                 locusIDs=locusNames, locusLength=locusLen)

write.nexus <- function(x, aln=FALSE, gapOpening=c(-18, -16), gapExtension=c(-2, -1), 
                        dir.out, fn, charset=FALSE, locusIDs, locusLength, verbose=FALSE) {
  #----------------------------------------------------------------------------#
  if (!requireNamespace("DECIPHER", quietly = TRUE)) {
    stop("Package 'DECIPHER' needed for this function to work. Please install it 
         either manually or using the function amplicR::setup().",
         call. = FALSE)
  }  
  #----------------------------------------------------------------------------#
  
  whatsx <- class(x)
  if(whatsx == "character") {
    x <- Biostrings::DNAStringSet(x)
    names(x) <- paste0("seq", seq_along(x))
    }
  if(whatsx == "DNAStringSet" & aln == FALSE) {
    x <- DECIPHER::AlignSeqs(x, gapOpening=gapOpening, gapExtension=gapExtension) 
  } else {
    if(whatsx == "DNAStringSet" & aln == TRUE & verbose == TRUE) 
      message("'aln' is set to TRUE, assuming 'x' is an alignment")
  }
  
  nex <- c("#NEXUS", 
           "begin taxa;",
           paste0("dimensions ntax=", length(x), ";"),
           "taxlabels", make.names(names(x), unique = TRUE), ";",
           "End;",
           "begin characters;",
           paste0("dimensions nchar=", IRanges::width(x)[1], ";"),
           "Format datatype=dna missing=? gap=- matchchar=.;",
           "Matrix",
           paste(make.names(names(x), unique = TRUE), x),
           ";",
           "End;")
 
  charsets <- NULL
  if(charset == TRUE) {
    pos <- 0
    charsets <- character(length(locusIDs))
    for(i in seq_along(charsets)) {
      charsets[i] <- paste("charset", locusIDs[i], "=", 1 + pos, "-", 
                           pos + locusLength[i], ";")
      pos <- pos + locusLength[i]
      charblock <- c("begin assumptions;", charsets, "end;")
    }
  }
  writeLines(c(nex, if(charset == TRUE) charblock), con=file.path(dir.out, fn))
}

#' Extract sequences from a derep object
#' 
#' This function is used internally to batch-extract sequences from a derep object
#' 
#' @param derepObj The derep-element (that is one element from the list that constitutes the derep list)
#' @return A character vector with the sequences
#' @export
#' 
getSeqFromDerep <- function(derepObj) {
  reads <- list(sequence=names(derepObj$uniques))
  return(reads)
}

#' Subset a derep-class based on abundance of the reads
#'
#' This function retains the reads that have their abundance >= \code{minAbund}.
#' It returns an object similar to a derep-class in its structure. Indeed it has
#' a \code{$uniques} and \code{quals} elemnt for the retained sequences, but,
#' rather than having \code{map} as third element in the list, it has
#' \code{rmed}, which is a vector with the index of the positions of the
#' \code{uniques} from the original object that were removed.
#' 
#' @param derep The derep object to subset
#' @param minAbund The threshold below which reads are removed
#' @return A list where the first two elements are as a derep-class, and the third is rmed. See details.
#' @export
subsetDerep <- function(derep, minAbund=2){
  sel <- derep$uniques >= minAbund
  uniques <- derep$uniques[sel]
  quals <- derep$quals[sel, ] 
  rmed <- seq_along(sel)[sel]
  return(list(uniques=uniques, quals=quals, rmed=rmed))
}

#' Build (possible) allele sequences based on genotypes
#'
#' Given the sequence of the reference allele, the genotypes, the SNPs (as the
#' substitute nucleotide) and their positions, this function builds the possible
#' allele sequences.
#'
#' The genotypes are to be codified as 0 for the reference allele, 2 for the
#' alternative allele, and 1 for the hetorozygous. The genotype vector has to be
#' a named integer vector where the names follow the convention.
#' "something-p-nb/na" something is usually the name of the locus, p is the
#' position of the SNP (as integer), nb is the nucleotide in the base allele and
#' na is the nucleotide of the alternative allele. For examples:
#' "100614668-2-A/C". Of all these, the critical elements are that nb and na has
#' to be in the second to last and penultimate position in the string.
#' @param baseAllele (Character) The sequence of the reference allele
#' @param genotypes (named integer) The genotypes (0,1,2) named with SNPs (see
#'   Details)
#' @param SNPpositions (integer) The position where the SNPs occur
#' @param lenAllele (integer) The length of the allele (as number of
#'   nucleotides)
#' @return A character vector with all the possible allele sequences
#' @export
#' @examples
#' SNPpositions <- list(
#' c(0,2,7),
#' c(0,1,2),
#' c(0,1,3),
#' c(1,2,7),
#' c(1,2,3),
#' c(2,4,6),
#' c(5,6,7)
#' )
#'
#'
#' baseAllele <- "AAAAAAAA"
#' genotypes <- c(2,1,1)
#' names(genotypes) <- paste0("something-p-", c("A/G", "A/C", "A/T"))
#' seqAlleles <- lapply(SNPpositions, make.alleles, baseAllele=baseAllele,
#'                      genotypes=genotypes,lenAllele=8)
#' seqAlleles                      
make.alleles <- function(baseAllele, genotypes, SNPpositions, lenAllele) {
  IUPAC <- c("AC", "AG", "AT", "CG", "CT", "GT", "CA", "GA", "TA", "GC", "TC", "TG")
  names(IUPAC) <- c("M", "R", "W", "S", "Y", "K", "M", "R", "W", "S", "Y", "K")
  # remember that SNP position are one behind because position 1 is indexed as 0
  breaks <- c(SNPpositions, 
              # if the last SNP position is at the end of the allele sequence, 
              # this needs to be dropped as it is already being taken care of by the next line
              if((max(SNPpositions) + 1) == lenAllele) {
                head(SNPpositions, -1) + 1
              } else {
                SNPpositions + 1  
              }
              , lenAllele)
  breaks <- sort(breaks)
  breaks <- breaks[!duplicated(breaks)] # if SNP are consecutive
  SNP1stPos <- which(SNPpositions == 0) # Check if one SNP is at the start
  if(length(SNP1stPos) == 1) breaks <- breaks[-which(breaks == 0)] # if so rm 0
  nsections <- length(breaks) 
  sections <- vector("list", length=nsections)
  seqAlleles <- ""
  s <- 1
  #section <- 1
  for(section in (seq_len(nsections))) {
    sections[[section]] <- substr(baseAllele, start=s, stop=breaks[section])
    if(breaks[section] %in% (SNPpositions + 1)) {
      whichGen <- which(SNPpositions + 1 == breaks[section])
      baseSNP <- substr(names(genotypes)[whichGen], 
                        start=nchar(names(genotypes)[whichGen]) - 2, 
                        stop=nchar(names(genotypes)[whichGen]) - 2)
      altSNP <- substr(names(genotypes)[whichGen], 
                       start=nchar(names(genotypes)[whichGen]), 
                       stop=nchar(names(genotypes)[whichGen]))
      if(is.na(genotypes[whichGen])) {
        sections[[section]] <- names(which(IUPAC == 
                                             paste(c(baseSNP, altSNP), collapse = "")))
      } else {
        if(genotypes[whichGen] == 2) {
          sections[[section]] <- altSNP
        } else {
          if(genotypes[whichGen] == 0) {
            sections[[section]] <- baseSNP
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
  return(seqAlleles)
}


#' Replace SNPs with IUPAC ambiguities
#' 
#' Replace SNPs with IUPAC ambiguities
#' 
#' See \code{make.alleles} for details on how \code{genotypes} should be formatted
#' 
#' @param seqAlleles possible sequences of the alleles or baseAllele (see \code{make.alleles} for details)
#' @inheritParams make.alleles
# Move function at the top
replaceSNP <- function(genotypes, seqAlleles, SNPpositions) {
  IUPAC <- c("AC", "AG", "AT", "CG", "CT", "GT", "CA", "GA", "TA", "GC", "TC", "TG")
  names(IUPAC) <- c("M", "R", "W", "S", "Y", "K", "M", "R", "W", "S", "Y", "K")
  seq <- seqAlleles[1]
  for(whichGen in seq_along(genotypes)) {
    baseSNP <- substr(names(genotypes)[whichGen], 
                      start=nchar(names(genotypes)[whichGen]) - 2, 
                      stop=nchar(names(genotypes)[whichGen]) - 2)
    altSNP <- substr(names(genotypes)[whichGen], 
                     start=nchar(names(genotypes)[whichGen]), 
                     stop=nchar(names(genotypes)[whichGen]))
    substr(seq, start=SNPpositions[whichGen] + 1, stop=SNPpositions[whichGen] + 1) <-
      names(which(IUPAC == paste(c(baseSNP, altSNP), collapse = "")))
  }
  return(seq)
}
