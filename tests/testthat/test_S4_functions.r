# Quick test
library(circulaR)
library(magrittr)
library(tidyverse)
library(AnnotationHub)
library(Gviz)
library(GenomicFeatures)
library(BSgenome.Hsapiens.NCBI.GRCh38)
data("supported.organisms")

ah <- AnnotationHub()
ahdb <- query(ah, pattern = c("Homo Sapiens", "EnsDb", 92))[[1]]
source("~/projects/circulaR/R/Experimental code/SetGeneric.R")
source("~/projects/circulaR/R/Experimental code/Sample.R")
source("~/projects/circulaR/R/Experimental code/ChimericExperiment.R")
source("~/projects/circulaR/R/Experimental code/S4Annotation.R")
source("~/projects/circulaR/R/Experimental code/S4Check.R")
source("~/projects/circulaR/R/Experimental code/S4Helper.R")

# Setup sample and experiments for testing
s1 <- circSample()

sample.basedir <- "~/ironside/disk1/Data/RNAseq/SEQ18.004_mapping_jointly_and_separately/aligned/A549/jointly/"

# Setup a sample and an experiment

sample.basedir <- "/disk1/Data/RNAseq/SEQ18.004_mapping_jointly_and_separately/aligned/A549/jointly/"

test_sample <- Sample(
  sample.id = "A549",
  organism = "Homo sapiens",
  gb = "hg38",
  chim.file = "~/ironside/disk1/Data/RNAseq/SEQ18.004_mapping_jointly_and_separately/aligned/A549/jointly/Chimeric.out.junction",
  bam.file = "~/ironside/disk1/Data/RNAseq/SEQ18.004_mapping_jointly_and_separately/aligned/A549/jointly/Aligned.sortedByCoord.out.bam",
  log.file = "~/ironside/disk1/Data/RNAseq/SEQ18.004_mapping_jointly_and_separately/aligned/A549/jointly/Log.final.out",
  sj.file = "~/ironside/disk1/Data/RNAseq/SEQ18.004_mapping_jointly_and_separately/aligned/A549/jointly/SJ.out.tab",
  firstread.firststrand = TRUE,
  paired.end = TRUE
  #, label = "character"
)


s1 <- readBSJdata(s1) # Alternatively: readBSJdata(s1, seqAcrossSj = F)
alignmentStats(s1)
kj <- constructSJDB(annotationDB = ahdb)
s1 <- compareToKnownJunctions(object = s1, known.junctions = kj)
s1 <- readLSJdata(s1)
experiment.basedir <- "~/ironside/disk1/Data/RNAseq/SEQ16.001_circRNA_in_ENCODE_cells/aligned/"

tmp <- bsj.reads(s1)

d.gr <- GRanges(
  seqnames = tmp$X1,
  ranges = IRanges(start = tmp$X2, width = 1),
  strand = tmp$X3
)
a.gr <- GRanges(
  seqnames = tmp$X4,
  ranges = IRanges(start = tmp$X5, width = 1),
  strand = tmp$X6
)

d <- resize(d.gr, width = 2, fix = "start") %>% Views(Hsapiens, .) %>% DNAStringSet %>% as.character
a <- resize(a.gr, width = 2, fix = "end") %>% Views(Hsapiens, .) %>% DNAStringSet %>% as.character
tmp$deduced.junc <- paste(d, a, sep = "/")
table(tmp$deduced.junc == "GT/AG")
table(tmp$deduced.junc == "CT/AC")
subset(tmp, X7 == 1)$deduced.junc %>% table
subset(tmp, X7 == 2)$deduced.junc %>% table
table(tmp$shiftAcceptorToNearestJ == 0 & tmp$shiftDonorToNearestJ == 0)

#tmp <- tmp[checkPEReads(tmp),]
tmp.shifted <- shiftAlignment(tmp)
table(tmp.shifted$deduced.junc == "GT/AG")
table(tmp.shifted$deduced.junc == "CT/AC")
subset(tmp.shifted, X7 == 1)$deduced.junc %>% table
subset(tmp.shifted, X7 == 2)$deduced.junc %>% table
table(tmp.shifted$shiftAcceptorToNearestJ == 0 & tmp.shifted$shiftDonorToNearestJ == 0)




### Experiment-wise ###
exp2 <- circExperiment(
  path = "~/ironside/disk1/Data/RNAseq/SEQ16.001_circRNA_in_ENCODE_cells/aligned/",
  name = "ENCODE data"
)

# Find and populate Experiment class with samples
exp2 <- locateSamples(exp2, organism = "Homo sapiens", genome = "hg38", firstread.firststrand = T)
exp2 <- readBSJdata(exp2, cores = 5)
exp2 <- compareToKnownJunctions(exp2, known.junctions = kj, cores = 5)
exp2 <- readLSJdata(exp2, cores = 5)


tmp <- samples(exp2) %>% lapply(bsj.reads) %>% do.call(rbind, .)

  lapply(shiftAlignment)
tmp2 <-
tmp2$X12 <- unlist(tmp2$X12)
tmp2$X14 <- unlist(tmp2$X14)

tmp3 <- dplyr::filter(tmp2, bsID == "hg38:X:140784661:140783175:+")
tmp3 <- dplyr::filter(tmp2, bsID == "hg38:1:200603225:200614406:-")

fb <- GRanges(
  seqnames = tmp3$X1,
  ranges = IRanges(start = tmp3$X11, width = sapply(tmp3$X12, parseCIGAR)),
  strand = tmp3$X3
)
sb <- GRanges(
  seqnames = tmp3$X1,
  ranges = IRanges(start = tmp3$X13, width = sapply(tmp3$X14, parseCIGAR)),
  strand = tmp3$X3
)




#######
tmp <- samples(exp2) %>% lapply(bsj.reads)
cbs <- do.call(rbind, tmp)

cbs$SHIFTED <- F

## Identify the different kinds shifts ##
shiftSize <- cbs$shiftDonorToNearestJ
dont.touch.these <- cbs$shiftAcceptorToNearestJ == 0 | cbs$shiftDonorToNearestJ == 0

# Group 1: Shifting that will bring both donor and acceptor in alignment
rescue.both <- knownJunctionInRepeatRegion(cbs, junc = "both") &
  cbs$shiftAcceptorToNearestJ == cbs$shiftDonorToNearestJ &
  cbs$donor.closest.type == "donor" &
  cbs$acceptor.closest.type == "acceptor" &
  !dont.touch.these
# Group 2: Shifting that will bring donor into alignment
rescue.donor.only <- knownJunctionInRepeatRegion(cbs, junc = "donor") &
  !knownJunctionInRepeatRegion(cbs, junc = "acceptor") &
  cbs$donor.closest.type == "donor" &
  !dont.touch.these
# group 3: Shifting that will bring acceptor into alignment
rescue.acceptor.only <- knownJunctionInRepeatRegion(cbs, junc = "acceptor") &
  !knownJunctionInRepeatRegion(cbs, junc = "donor") &
  cbs$donor.closest.type == "acceptor" &
  !dont.touch.these
message("Reads that can be shifted: ")
message("- Both OK: ", length(which(rescue.both)))
message("- Donor OK: ", length(which(rescue.donor.only)))
message("- Acceptor OK: ", length(which(rescue.acceptor.only)))

## Check the CIGARS ##
## BOTH DONOR & ACCEPTOR ##
# Plus
cbs$X12[rescue.both & cbs$X3 == "+"] %>% lapply(parseCIGAR, returnLength = F) %>% sapply(function(x)all(names(x[1:2]) == c("S", "M"))) %>% table
cbs$X14[rescue.both & cbs$X3 == "+"] %>% lapply(parseCIGAR, returnLength = F) %>% sapply(function(x)all(names(x[(length(x)-1):(length(x))]) == c("M", "S"))) %>% table

# Minus
cbs$X12[rescue.both & cbs$X3 == "-"] %>% lapply(parseCIGAR, returnLength = F) %>% sapply(function(x)all(names(x[(length(x)-1):(length(x))]) == c("M", "S"))) %>% table
cbs$X14[rescue.both & cbs$X3 == "-"] %>% lapply(parseCIGAR, returnLength = F) %>% sapply(function(x)all(names(x[1:2]) == c("S", "M"))) %>% table

## ONLY DONOR ##
# Plus
cbs$X12[rescue.donor.only & cbs$X3 == "+"] %>% lapply(parseCIGAR, returnLength = F) %>% sapply(function(x)all(names(x[1:2]) == c("S", "M"))) %>% table
cbs$X14[rescue.donor.only & cbs$X3 == "+"] %>% lapply(parseCIGAR, returnLength = F) %>% sapply(function(x)all(names(x[(length(x)-1):(length(x))]) == c("M", "S"))) %>% table

# Minus
cbs$X12[rescue.donor.only & cbs$X3 == "-"] %>% lapply(parseCIGAR, returnLength = F) %>% sapply(function(x)all(names(x[(length(x)-1):(length(x))]) == c("M", "S"))) %>% table
cbs$X14[rescue.donor.only & cbs$X3 == "-"] %>% lapply(parseCIGAR, returnLength = F) %>% sapply(function(x)all(names(x[1:2]) == c("S", "M"))) %>% table

## ONLY ACCEPTOR ##
# Plus
cbs$X12[rescue.acceptor.only & cbs$X3 == "+"] %>% lapply(parseCIGAR, returnLength = F) %>% sapply(function(x)all(names(x[1:2]) == c("S", "M"))) %>% table
cbs$X14[rescue.acceptor.only & cbs$X3 == "+"] %>% lapply(parseCIGAR, returnLength = F) %>% sapply(function(x)all(names(x[(length(x)-1):(length(x))]) == c("M", "S"))) %>% table

# Minus
cbs$X12[rescue.acceptor.only & cbs$X3 == "-"] %>% lapply(parseCIGAR, returnLength = F) %>% sapply(function(x)all(names(x[(length(x)-1):(length(x))]) == c("M", "S"))) %>% table
cbs$X14[rescue.acceptor.only & cbs$X3 == "-"] %>% lapply(parseCIGAR, returnLength = F) %>% sapply(function(x)all(names(x[1:2]) == c("S", "M"))) %>% table

