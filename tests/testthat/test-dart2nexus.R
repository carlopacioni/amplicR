temp <- tempdir(check = TRUE)
library(amplicR)
library(data.table)
amplicR::setup()
library(adegenet)

# Create some seqs
barcode1 <- "TCGA"
barcode2 <- "AGCT"
baseRead <- paste(rep("A", 5), collapse = "")
altRead <- "AACAT"
tail <- "CCC"
toFile1 <- paste0(barcode1, c(baseRead, altRead), tail)
toFile2 <- paste0(barcode2, c(baseRead, altRead), tail)

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
dat <- 
gl <- new("genlight", list(S1=c(0,1,1), S2=c(0,0,1)))
gl
SNPname <- 100614668
locNames(gl) <- paste(SNPname, c(0,2,4), c("A/G", "A/C", "A/T"), sep="-")
as.matrix(gl)

# loc.metrics
loc.metrics <- data.frame(CloneID=rep(SNPname, 3),
                          AlleleSequence=rep(paste0(baseRead, tail), 3),
                          TrimmedSequence=rep(baseRead, 3),
                          SnpPosition=c(0,2,4))
other(gl) <- list(loc.metrics=loc.metrics)

# csv
targetid <- as.numeric(fns)
genotype	<- rep("S1", 2)
barcode <- c(barcode1, barcode2)
write.csv(data.frame(targetid, genotype, barcode), file = file.path(temp, "targets_test.csv"))

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
all(nalleles == 4)
all.equal(seqAlleles[[1]], c("GAAAAAAA", "GAAAAAAT", "GACAAAAA", "GACAAAAT"))
all.equal(seqAlleles[[2]], c("GAAAAAAA", "GATAAAAA", "GCAAAAAA", "GCTAAAAA"))
all.equal(seqAlleles[[3]], c("GAAAAAAA", "GAATAAAA", "GCAAAAAA", "GCATAAAA"))
all.equal(seqAlleles[[4]], c("AGAAAAAA", "AGAAAAAT", "AGCAAAAA", "AGCAAAAT"))
all.equal(seqAlleles[[5]], c("AGAAAAAA", "AGATAAAA", "AGCAAAAA", "AGCTAAAA"))
all.equal(seqAlleles[[6]], c("AAGAAAAA", "AAGAAATA", "AAGACAAA", "AAGACATA"))
all.equal(seqAlleles[[7]], c("AAAAAGAA", "AAAAAGAT", "AAAAAGCA", "AAAAAGCT"))

# run dart2nexus
oldwd <- getwd()
setwd(temp)

fastq.dir.in <- temp
proc.data.out <- "../data_processing/"
minQ=25
truncQ=0
minAbund=NULL
min.nSNPs=3
dir.out="Processed_data"
singleAllele=TRUE
nCPUs="auto"
i<- 1
locus <- target.loci[1]
dada=FALSE
minLen=4
sampleID<-samplesIDs[1]
target <- targets[1]

#------------#
dada2::fastqFilter(fn = file.path(temp, paste0(fns[1],".fastq.gz")),
                   fout = file.path(temp, paste0("Test",".fastq.gz")),
                   maxN=0,
                   minQ=15,
                   maxEE=Inf,
                   truncQ=10,
                   minLen=3,
                   compress=TRUE,
                   OMP = FALSE,
                   verbose=FALSE)
rf <- readFastq(dirPath = temp, pattern = "^Test")
rf
rf@quality
