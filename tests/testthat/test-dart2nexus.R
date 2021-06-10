library(amplicR, quietly=TRUE)
library(adegenet, quietly = TRUE)
context("Test dart2nexus")

test_that("dart2nexus", {
  
temp <- tempdir(check = TRUE)
suppressMessages(
amplicR::setup()
)

# Create some seqs
# Seqs created for one samples, two locus over two different files
# First locus
barcode1 <- "TCGA"
barcode2 <- "AGCT"
baseRead <- paste(rep("A", 5), collapse = "")
altRead <- "AACAT"
tail <- "CCC"
# Second locus
baseReadbis <- paste(rep("G", 5), collapse = "")
altReadbis <- "GGCGT"
tailbis <- "AAA"
# Thirs locus
baseReadtris <- paste(rep("C", 5), collapse = "")
altReadtris <- "CCGCT"
altReadtrisWrong <- "CCCCT"
tailtris <- "TTT"

toFile1 <- c(
  paste0(barcode1, c(baseRead, altRead), tail), # First locus
  paste0(barcode1, c(baseReadbis, altReadbis), tailbis), # Second locus
  rep(paste0(barcode1, c(baseReadtris, altReadtris), tailtris), # Third locus
      2), # so that there is higher abundance of the correct ones
  paste0(barcode1, c(baseReadtris, altReadtrisWrong), tailtris) # Wrong third locus
  )                                                             # for S1
toFile2 <- c(
  paste0(barcode2, c(baseRead, altRead), tail),
  paste0(barcode2, c(baseReadbis, altReadbis), tailbis),
  paste0(barcode2, c(baseReadtris, altReadtris), tailtris)
)

seq1 <- DNAStringSet(toFile1)
seq2 <- DNAStringSet(toFile2)
qual <- Biostrings::quality(NumericQuality(rep(30L, width(toFile1)[1])))
qual <- IlluminaQuality(qual)
sr1 <- ShortRead::ShortReadQ(sread=DNAStringSet(toFile1), 
                             quality=BStringSet(
                               rep(qual, length(toFile1))), 
                             id=BStringSet(letters[seq_len(length(toFile1))]))

qual <- Biostrings::quality(NumericQuality(rep(30L, width(toFile2)[1])))
qual <- IlluminaQuality(qual)
sr2 <- ShortRead::ShortReadQ(sread=DNAStringSet(toFile2), 
                             quality=BStringSet(
                               rep(qual, length(toFile2))), 
                             id=BStringSet(letters[seq_len(length(toFile2))]))

fns <- c("1518950", "1518951")

# Write to fastq  files
writeFastq(object=sr1, file=file.path(temp, paste0(fns[1],".fastq.gz")))
writeFastq(object=sr2, file=file.path(temp, paste0(fns[2],".fastq.gz")))

# Generate a gl
gl <- new("genlight", list(S1=c(0,1,1, 0,1,1, 0,1,1), S2=c(0,0,1, 0,0,1, 0,0,1)))

SNPname <- c(100614668, 100614670, 100614671)
locNames(gl) <- c(
  paste(SNPname[1], c(0,2,4), c("A/G", "A/C", "A/T"), sep="-"),
  paste(SNPname[2], c(0,2,4), c("A/G", "A/C", "A/T"), sep="-"),
  paste(SNPname[3], c(0,2,4), c("C/G", "C/G", "C/T"), sep="-")
)
#as.matrix(gl)

# loc.metrics
loc.metrics <- data.frame(CloneID=rep(SNPname, each=3),
                          AlleleSequence=c(
                            rep(paste0(baseRead, tail), 3),
                            rep(paste0(baseReadbis, tailbis), 3),
                            rep(paste0(baseReadtris, tailtris), 3)),
                          TrimmedSequence=c(
                            rep(baseRead, 3), rep(baseReadbis, 3),
                            rep(baseReadtris, 3)),
                          SnpPosition=rep(c(0,2,4), length(SNPname)))
other(gl) <- list(loc.metrics=loc.metrics)

# csv
targetid <- as.numeric(fns)
genotype	<- rep("S1", 2)
barcode <- c(barcode1, barcode2)
write.csv(data.frame(targetid, genotype, barcode), file = file.path(temp, "targets_test.csv"))

# run dart2nexus
expect_warning(
res <- dart2nexus(gl, dir.in=temp, min.nSNPs=3, 
                       minLen=4, truncQ=0, minQ=25,
                       singleAllele=TRUE, dada=FALSE, 
                       nCPUs="auto")
)
#debug(dart2nexus)
together <- c(res[["Allele1"]], res[["Allele2"]])
locus1 <- DNAStringSet(together, start = 1, end = 5)
locus2 <- DNAStringSet(together, start = 6, end = 10)
locus3 <- DNAStringSet(together, start = 11, end = 15)
suppressWarnings(
distLocus1 <- ShortRead::srdistance(locus1, DNAStringSet(c("AACAT", "AAAAA", "AAAAT")))
)
suppressWarnings(
  distLocus2 <- ShortRead::srdistance(locus2, DNAStringSet(c("AGMGW", "AGAGA", "AGAGT")))
)
suppressWarnings(
  distLocus3 <- ShortRead::srdistance(locus3, DNAStringSet(c("CCCCC", "CCCCT", "CCGCT")))
)
sumDistLocus1 <- sapply(distLocus1, sum)
sumDistLocus2 <- sapply(distLocus2, sum)
sumDistLocus3 <- sapply(distLocus3, sum)

expect_equal(unname(sumDistLocus1), c(5,3,3))
expect_equal(unname(sumDistLocus2), c(4,3,3))
expect_equal(unname(sumDistLocus3), c(3,3,5))

resNoSeqs <- dart2nexus(gl, dir.in=NA, dir.out=temp, min.nSNPs=3, 
                  minLen=4, truncQ=0, minQ=25,
                  singleAllele=TRUE, dada=FALSE, 
                  nCPUs="auto")
S1 <- "AAMAWAGMGWCCSCY"
names(S1) <- "S1"
expect_equal(as.character(resNoSeqs[["Allele1"]][1]), S1)

S2 <- "AAAAAAGAGACCCCC"
names(S2) <- "S2"
S2bis <- "AAAATAGAGACCCCC"
names(S2bis) <- "S2"
S2tris <- "AAAAAAGAGTCCCCC"
names(S2tris) <- "S2"
S2tetra <- "AAAATAGAGTCCCCC"
names(S2tetra) <- "S2"

S2penta <- "AAAAAAGAGACCCCT"
names(S2penta) <- "S2"
S2esa <- "AAAATAGAGACCCCT"
names(S2esa) <- "S2"
S2hepta <- "AAAAAAGAGTCCCCT"
names(S2hepta) <- "S2"
S2octa <- "AAAATAGAGTCCCCT"
names(S2octa) <- "S2"

expect_true(
  as.character(resNoSeqs[["Allele1"]][2]) == S2 | 
  as.character(resNoSeqs[["Allele1"]][2]) == S2bis |
  as.character(resNoSeqs[["Allele1"]][2]) == S2tris |
  as.character(resNoSeqs[["Allele1"]][2]) == S2tetra |
    as.character(resNoSeqs[["Allele1"]][2]) == S2penta | 
    as.character(resNoSeqs[["Allele1"]][2]) == S2esa |
    as.character(resNoSeqs[["Allele1"]][2]) == S2hepta |
    as.character(resNoSeqs[["Allele1"]][2]) == S2octa
  )

# Test make.alleles
SNPpositions <- list(
  c(0,2,7), 
  c(0,1,2),
  c(0,1,3),
  c(1,2,7),
  c(1,2,3),
  c(2,4,6),
  c(5,6,7)
)


baseAllele <- "AAAAAAAA"
genotypes <- c(2,1,1)
names(genotypes) <- paste0("something-p-", c("A/G", "A/C", "A/T"))
seqAlleles <- lapply(SNPpositions, make.alleles, baseAllele=baseAllele, 
                           genotypes=genotypes,lenAllele=8)
nalleles <- sapply(seqAlleles, length)

expect_equal(nalleles, rep(4, 7))
expect_equal(seqAlleles[[1]], c("GAAAAAAA", "GAAAAAAT", "GACAAAAA", "GACAAAAT"))
expect_equal(seqAlleles[[2]], c("GAAAAAAA", "GATAAAAA", "GCAAAAAA", "GCTAAAAA"))
expect_equal(seqAlleles[[3]], c("GAAAAAAA", "GAATAAAA", "GCAAAAAA", "GCATAAAA"))
expect_equal(seqAlleles[[4]], c("AGAAAAAA", "AGAAAAAT", "AGCAAAAA", "AGCAAAAT"))
expect_equal(seqAlleles[[5]], c("AGAAAAAA", "AGATAAAA", "AGCAAAAA", "AGCTAAAA"))
expect_equal(seqAlleles[[6]], c("AAGAAAAA", "AAGAAATA", "AAGACAAA", "AAGACATA"))
expect_equal(seqAlleles[[7]], c("AAAAAGAA", "AAAAAGAT", "AAAAAGCA", "AAAAAGCT"))

unlink(temp, recursive=TRUE)
})

