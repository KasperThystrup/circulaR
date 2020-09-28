# Import
context(desc = "Chimeric functions")
library(circulaR)
library(GenomicRanges)
library(IRanges)

test_that(desc = "Check that output of generateCandidateFilters is correct",{
  df <- data.frame(
    X1 = c("A", "B", "A", "C", "C", "A"),
    X2 = c(30, 35, 35, 30, 30, 20),
    X3 = c("-", "-", "+", "-", "+", "+"),
    X4 = c("A", "B", "A", "A", "C", "A"),
    X5 = c(35, 35, 30, 35, 35, 10),
    X6 = c("-", "-", "+", "-", "-", "+"),
    X7 = c(-1, -1, 0, 0, 0, 4)
  )
  expect_equal(generateCandidateFilters(df, backSpliceDist = 10, FALSE), c(T, F, T, F, F, T)) # Check for ALL chim
  expect_equal(generateCandidateFilters(df, backSpliceDist  = 5, FALSE), c(T, F, T, F, F, F)) # Check backSpliceDist
  expect_equal(generateCandidateFilters(df, backSpliceDist  = 0, FALSE), c(F, F, F, F, F, F)) # Check zero dist and filtration
  expect_equal(generateCandidateFilters(df, backSpliceDist = 10,  TRUE), c(F, F, T, F, F, T)) # Check dismissal of encompassing reads
})


test_that(desc = "GRanges object works", {
  expect_equal(junctionsToGRanges(seqnames = 2, start = 5, width = 1, strand = "+"), GRanges(ranges = IRanges(start = 5, end = 5), strand = "+", seqnames = 2))
  expect_equal(junctionsToGRanges(seqnames = 2, start = 5, end = 5, strand = "+"), GRanges(ranges = IRanges(start = 5, width = 1), strand = "+", seqnames = 2))
  expect_equal(junctionsToGRanges(seqnames = 2, start = 5, end = 5, strand = "+", test = "33"), GRanges(ranges = IRanges(start = 5, width = 1), strand = "+", seqnames = 2, test = "33"))
})


# #### THEESE FOLLOWING TESTINGS is not at all complete, and needs further clarification and better tests
# test_that(desc = "Nearest works"{
#   gr1 <- GRanges(ranges = IRanges(start = 5, end = 5), strand = "+", seqnames = 2)
#   gr2 <- GRanges(ranges = IRanges(start = 3, end = 5), strand = "-", seqnames = 1)
#   gr3 <- GRanges(ranges = IRanges(start = 1, end = 2), strand = "+", seqnames = 1)
#   gr4 <- GRanges(ranges = IRanges(start = 1, end = 2), strand = "+", seqnames = 3)
#   kj <- GRanges(
#     seqnames = c(1,1,1,2),
#     strand = c("-", "+", "-", "+"),
#     type = c("acceptor", "donor", "acceptor", "donor"),
#     g_id = c("ENSG007", "ENSG006", "ENSG010", "ENSG200"),
#     tx_id = c("ENST001", "ENST002", "ENST003", "ENST004"),
#     jID = c("J1", "J2", "J3", "J4"),
#     ranges = c(IRanges(
#       start = c(1, 2, 3, 4),
#       width = c(1, 2, 3, 4)
#     )))
#   #expect_equal(as.matrix(findNearest(backsplice_junctions = gr1, junctions_annotation = kj), dimnames = NULL), matrix(c("1", "4", "J4", "donor", "+", "1"), nrow = 1))
#
# })
#
# test_that(desc = "Output of addKnownJunctions return correct data_structure", {
#   gr4 <- GRanges(ranges = IRanges(start = 1, end = 2), strand = "+", seqnames = 3)
#   kj <- GRanges(
#     seqnames = c(1,1,1,2),
#     strand = c("-", "+", "-", "+"),
#     type = c("acceptor", "donor", "acceptor", "donor"),
#     g_id = c("ENSG007", "ENSG006", "ENSG010", "ENSG200"),
#     tx_id = c("ENST001", "ENST002", "ENST003", "ENST004"),
#     jID = c("J1", "J2", "J3", "J4"),
#     ranges = c(IRanges(
#       start = c(1, 2, 3, 4),
#       width = c(1, 2, 3, 4)
#     )))
#   #expect_equal(colnames(addKnownJunctions(cbs = gr4, kj = kj)), c("A", "B", "C"))
# })
