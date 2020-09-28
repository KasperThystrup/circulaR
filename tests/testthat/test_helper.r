# Import
context(desc = "Helper functions")
library(circulaR)
library(GenomicRanges)
library(IRanges)
library(GenomicFeatures)
library(AnnotationHub)
library(ensembldb)

# Test that the output of constructSJDB is as expected
l <- c(
  "2\t38748082\t-\t2\t38750054\t-\t0\t0\t7\tD74RYQN1:328:C480EACXX:2:1103:16101:50753\t38749667\t39M308N40M22S\t38748083\t79S22M1423p101M\tTRUE\t1\tJ495543\tdonor\t41\tJ495531\tdonor\t3",
  "9\t128618048\t+\t9\t128613401\t+\t0\t0\t4\tD74RYQN1:328:C480EACXX:2:1103:8540:51242\t128613402\t35S66M2259p101M\t128618013\t35M66S\tTRUE\t2\tJ891339\tacceptor\t21\tJ891348\tdonor\t-61",
  "6\t75703073\t+\t6\t75702644\t+\t2\t0\t5\tD74RYQN1:328:C480EACXX:2:1103:5721:51585\t75702645\t73S28M\t75702753\t101M146p73M28S\tTRUE\t3\tJ765625\tacceptor\t0\tJ765626\tdonor\t0",
  "19\t14480421\t-\t19\t14480670\t-\t0\t0\t3\tD74RYQN1:328:C480EACXX:2:1103:19372:51552\t14480409\t77M107N24M-135p4M107N77M20S\t14480422\t81S20M\tTRUE\t4\tJ446955\tacceptor\t73\tJ446953\tacceptor\t-65",
  "21\t34385552\t+\t21\t34375547\t+\t0\t0\t5\tD74RYQN1:328:C480EACXX:2:1103:2995:53019\t34375548\t72S25M-19p14M3860N87M\t34385480\t72M29S\tTRUE\t9\tJ580138\tacceptor\t-13\tJ580146\tacceptor\t75"
)
test_data <- readChimFile(file = paste0(l, collapse = "\n"))

test_that(desc = "Check countSumChimFile",{
  expect_error(countSumChimFile(df = "Test"))
  expect_equal(class(countSumChimFile(df = test_data)), c("tbl_df", "tbl", "data.frame"))
  expect_equal(colnames(countSumChimFile(df = test_data)), c("bsID", "count"))
})


### Test various helper functions that mofify the bsID
id.on.plus <- "hg38:11:33287512:33286412:+" # HIPK3
id.on.minus <- "hg38:4:186706562:186709846:-" # FAT1

test_that(desc = "Check bsid2bed",{
  expect_equal(bsid2bed(id.on.plus, UCSCstyle = T), "chr11|33286412-33287511|+")
  expect_equal(bsid2bed(id.on.plus, UCSCstyle = F), "11|33286412-33287511|+")
  expect_equal(bsid2bed(id.on.minus, UCSCstyle = F), "4|186706562-186709845|-")
})

test_that(desc = "Check bsid2gr",{
  expect_equal(class(bsid2gr(c(id.on.plus, id.on.minus))), structure("GRanges", package = "GenomicRanges"))
  expect_equal(width(bsid2gr(c(id.on.plus, id.on.minus))), c(1101L, 3285L))
  expect_equal(start(bsid2gr(c(id.on.plus, id.on.minus))), c(33286412L, 186706562L))
})

test_that(desc = "Check bsid2junc",{
  expect_equal(class(bsid2junc(c(id.on.plus, id.on.minus),junc.type = "acceptor")), structure("GRanges", package = "GenomicRanges"))
  expect_equal(start(bsid2junc(c(id.on.plus, id.on.minus),junc.type = "acceptor")), c(33286412L, 186709846L))
})

test_that(desc = "Check constructBsId",{
  expect_equal(
    constructBsId(df = test_data, gb = "hg38"),
    c("hg38:2:38748082:38750054:-", "hg38:9:128618048:128613401:+",
      "hg38:6:75703073:75702644:+", "hg38:19:14480421:14480670:-",
      "hg38:21:34385552:34375547:+"
    )
  )
})

test_that(desc = "Check parseCIGAR",{
  cs <- c(
    "37S22M36247N33M652N9M-36952p11M36247N33M652N47M22921N10M",
    "43M35495N38M30262N20M-65836p21M35495N38M30262N22M20S",
    "35S66M2259p101M",
    "73S28M",
    "35M66S"
  )
  expect_warning(parseCIGAR(cig = cs))
  expect_equivalent(sapply(cs, parseCIGAR), c(59932, 65860, 2426, 28, 35))
  expect_equal(
    parseCIGAR(cig = cs[1], returnLength = F),
    structure(
      c(37, 22, 36247, 33, 652, 9, -36952, 11, 36247, 33, 652, 47, 22921, 10),
      .Names = c("S", "M", "N", "M", "N", "M", "p", "M", "N", "M", "N", "M", "N", "M")
    )
  )
  expect_equal(
    parseCIGAR(cig = cs[5], returnLength = F),
    structure(
      c(35, 66),
      .Names = c("M", "S")
    )
  )
})

test_that(desc = "junctionsToGRanges", {
  expect_equal(junctionsToGRanges(seqnames = 2, start = 5, width = 1, strand = "+"), GRanges(ranges = IRanges(start = 5, end = 5), strand = "+", seqnames = 2))
  expect_equal(junctionsToGRanges(seqnames = 1, start = 2, end = 4, strand = "-"), GRanges(ranges = IRanges(start = 2, width = 3), strand = "-", seqnames = 1))
  expect_equal(junctionsToGRanges(seqnames = 2, start = 5, end = 5, strand = "+", test = "33"), GRanges(ranges = IRanges(start = 5, width = 1), strand = "+", seqnames = 2, test = "33"))
  expect_error(junctionsToGRanges(seqnames = 2))
  expect_error(junctionsToGRanges(start = 2, end = 4, stand = "-"))
  expect_error(junctionsToGRanges(seqnames = 1, start = 2, end = 4))

})
