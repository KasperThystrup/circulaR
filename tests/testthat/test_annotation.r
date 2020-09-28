# Import
context(desc = "circRNA annotation functions")
library(circulaR)
library(GenomicRanges)
library(IRanges)
library(GenomicFeatures)
library(AnnotationHub)
library(ensembldb)

#### Is is a bad idea to construct these things outside test_that???
# Load the txdb and construct SJDB
# txdb <- loadDb("../../misc/test_data/testTxDb.sqlite")
txdb <- loadDb("misc/test_data/testTxDb.sqlite") ### How was paths working again?
txdb.sj <- constructSJDB(db = txdb)

# Load the ahdb and construct SJDB
# ahdb <- EnsDb("../../misc/test_data/testAhDb.sqlite")
ahdb <- EnsDb("misc/test_data/testAhDb.sqlite")
ahdb.sj <- constructSJDB(db = ahdb)

# Define three "artificial" circRNAs, one for each gene in the test databases
circRNAs <- c(
  PVT1 = "hg38:8:127890588:127890999:+",
  TP53 = "hg38:17:7687539:7666085:-",
  GAPDH = "hg38:12:6534511:6538360:+"
)

# Test that the output of constructSJDB is as expected
test_that(desc = "Check constructSJDB",{
  expect_equal(class(txdb.sj), structure("GRanges", package = "GenomicRanges"))
  expect_equal(class(ahdb.sj), structure("GRanges", package = "GenomicRanges"))
  expect_error(constructSJDB(db = "a string"))

  expect_equal(GenomicRanges::seqnames(txdb.sj), Rle(factor(rep(c("8", "12", "17"), c(75, 42, 61)), levels = c("8", "12", "17"))))
  expect_equal(GenomicRanges::seqnames(ahdb.sj), Rle(factor(rep(c("12", "17", "8"), c(42, 61, 75)), levels = c("8", "12", "17"))))

  expect_equal(GenomicRanges::start(txdb.sj)[c(1,length(txdb.sj))], c(127794532L, 7687551L))
  expect_equal(GenomicRanges::start(ahdb.sj)[c(1,length(txdb.sj))], c(6533926L, 128101254L))

  expect_equal(table(txdb.sj$type), structure(c(88L, 90L), .Dim = 2L, .Dimnames = structure(list(c("acceptor", "donor")), .Names = ""), class = "table"))
  expect_equal(table(ahdb.sj$type), structure(c(88L, 90L), .Dim = 2L, .Dimnames = structure(list(c("acceptor", "donor")), .Names = ""), class = "table"))
})
#
# test_that(desc = "Check constructSJDB",{
#   expect_equal(class(txdb.sj), structure("GRanges", package = "GenomicRanges"))
#   expect_equal(class(ahdb.sj), structure("GRanges", package = "GenomicRanges"))
#   expect_error(constructSJDB(db = "a string"))
#
#   expect_equal(GenomicRanges::seqnames(txdb.sj), Rle(factor(rep(c("8", "12", "17"), c(75, 42, 61)), levels = c("8", "12", "17"))))
#   expect_equal(GenomicRanges::seqnames(ahdb.sj), Rle(factor(rep(c("12", "17", "8"), c(42, 61, 75)), levels = c("8", "12", "17"))))
#
#   expect_equal(GenomicRanges::start(txdb.sj)[c(1,length(txdb.sj))], c(127794532L, 7687551L))
#   expect_equal(GenomicRanges::start(ahdb.sj)[c(1,length(txdb.sj))], c(6533926L, 128101254L))
#
#   expect_equal(table(txdb.sj$type), structure(c(88L, 90L), .Dim = 2L, .Dimnames = structure(list(c("acceptor", "donor")), .Names = ""), class = "table"))
#   expect_equal(table(ahdb.sj$type), structure(c(88L, 90L), .Dim = 2L, .Dimnames = structure(list(c("acceptor", "donor")), .Names = ""), class = "table"))
# })

# Test that annotation function returns expected output
annot.ahdb <- annotateByOverlap(bsids = circRNAs, db = ahdb)
annot.txdb <- annotateByOverlap(bsids = circRNAs, db = txdb)

test_that(desc = "Check annotateByOverlap",{
  expect_equal(class(annot.ahdb), c("tbl_df", "tbl", "data.frame"))
  expect_equal(class(annot.txdb), c("tbl_df", "tbl", "data.frame"))

  expect_equal(annot.ahdb$GENENAME, c("PVT1", "TP53", "GAPDH"))
  expect_equal(annot.txdb$GENEID, c("ENSG00000249859", "ENSG00000141510", "ENSG00000111640"))
})

# The addKnownJunctions function
## Generate data that will be used for a couple of tests
l <- c(
  "8\t127890999\t+\t8\t127890588\t+\t2\t1\t2\tDE18INS60503:186:CAJ2RANXX:4:1114:18389:71922\t127890589\t43S57M13p100M\t127890956\t43M57S",
  "17\t7673548\t-\t17\t7674896\t-\t0\t0\t1\tDE18INS60503:186:CAJ2RANXX:5:1207:18411:4951\t7674213\t78M568N22M-24p39M49S\t7673549\t51S49M",
  "17\t7665901\t-\t17\t7668586\t-\t0\t0\t6\tDE18INS60503:186:CAJ2RANXX:5:1104:18038:95844\t7668518\t68M32S\t7665902\t68S32M-1p100M",
  "12\t6537000\t+\t12\t6536715\t+\t0\t0\t6\tDE18INS60503:186:CAJ2RANXX:4:1307:1457:82210\t6536716\t69S31M-37p81M129N12M\t6536931\t69M31S",
  "12\t6538333\t+\t12\t6537979\t+\t0\t0\t0\tDE18INS60503:186:CAJ2RANXX:5:1205:6960:91249\t6537980\t17S17M104N66M34p100M\t6538316\t17M83S",
  "12\t6538350\t+\t12\t6538184\t+\t0\t0\t0\tDE18INS60503:186:CAJ2RANXX:4:2103:3749:54615\t6538185\t41S59M\t6538185\t2S98M26p41M59S",
  "12\t116230532\t-\t12\t116237706\t-\t2\t1\t3\tD74RYQN1:328:C480EACXX:2:2314:7572:36291\t116237660\t46M55S\t116230533\t46S55M6895p101M", # Rescued read, and discarded donor site
  "12\t124398118\t-\t12\t124422556\t-\t2\t2\t2\tD74RYQN1:328:C480EACXX:2:1310:21023:71616\t124402541\t21M17395N80M-6p25M2445N55M21S\t124398119\t80S21M", # Discarded Acceptor
  "8\t55762851\t-\t8\t55773354\t-\t0\t0\t2\tD74RYQN1:328:C480EACXX:2:2306:9429:91923\t55773293\t61M40S\t55762852\t61S40M1054p35M9288N66M", # Unresolved Donor
  "8\t55762851\t-\t8\t55773354\t-\t0\t0\t2\tD74RYQN1:328:C480EACXX:2:1207:2566:94988\t55773293\t61M40S\t55762852\t61S40M1054p35M9288N66M" # Unresolved Acceptor
)

chim <- readChimFile(file = paste0(l, collapse = "\n"))
resolved <- addKnownJunctions(cbs = chim, kj = ahdb.sj)

test_that(desc = "Check addKnownJunctions",{
  expect_equal(resolved$shiftAcceptorToNearestJ, c(0L, 38L, 166L, 32L, -18L, 84L, 0L))
  expect_equal(resolved$shiftDonorToNearestJ, c(0L, 14L, -184L, -11L, -5L, 9L, -6935L))

  expect_equal(resolved$acceptor.closest.type, c("acceptor", "donor", "donor", "acceptor", "donor", "acceptor", "acceptor"))
  expect_equal(resolved$donor.closest.type, c("donor", "donor", "donor", "donor", "donor", "donor", "donor"))
})


