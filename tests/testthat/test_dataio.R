# Import
context(desc = "Data I/O functions")
library(circulaR)
library(GenomicRanges)
library(IRanges)

test_that(desc = "Check readSJout",{
  # Generate a small sample of data
  l <- c(
    "1\t19389\t189779\t1\t3\t0\t22\t0\t46", "1\t804223\t804775\t1\t1\t1\t3\t0\t3", "1\t829105\t847653\t1\t1\t1\t1\t0\t36",
    "1\t847807\t850177\t1\t1\t0\t14\t0\t42", "1\t17056\t17914\t2\t2\t0\t0\t6\t25", "1\t17369\t17605\t2\t2\t1\t14\t56\t47",
    "1\t17743\t17914\t2\t2\t1\t23\t597\t50", "1\t18062\t18267\t2\t2\t1\t0\t208\t50", "1\t21550\t191785\t0\t0\t0\t6\t0\t33",
    "1\t1044370\t1045462\t0\t0\t0\t8\t0\t37"
  )
  sampleLinSJ <- textConnection(paste0(l, collapse = "\n"))

  sj.data <- readSJout(fn = sampleLinSJ)
  expect_equal(class(sj.data), structure("GRanges", package = "GenomicRanges"))
  expect_equal(length(sj.data), 10)
  expect_equal(GenomicRanges::start(sj.data), c(19389L, 804223L, 829105L, 847807L, 17056L, 17369L, 17743L, 18062L, 21550L, 1044370L))
  expect_equal(GenomicRanges::end(sj.data), c(189779L, 804775L, 847653L, 850177L, 17914L, 17605L, 17914L, 18267L, 191785L, 1045462L))
  expect_equal(GenomicRanges::strand(sj.data), Rle(factor(c("+", "+", "+", "+", "-", "-", "-", "-", "*", "*"), levels = c("+", "-", "*"))))
  expect_equal(colnames(mcols(sj.data)), c("int.motif", "annotated", "read.count.unique", "read.count.multi", "max.overhang"))
})


test_that(desc = "Check readBSJout",{
  l <- c("bsID\tcount", "hg38:17:8518655:8521151:-\t31", "hg38:7:107159582:107159492:+\t2",
         "hg38:17:42816025:42818856:-\t7", "hg38:19:2067967:2068296:-\t20", "hg38:19:49798871:49782626:+\t9",
         "hg38:10:11869812:11869601:+\t14", "hg38:19:10690351:10689961:+\t14", "hg38:15:74055768:74055408:+\t2",
         "hg38:10:68960250:68959805:+\t21", "hg38:3:127661073:127660597:+\t1"
        )

  sampleBSJ <- textConnection(paste0(l, collapse = "\n"))

  bsj.data <- readBSJout(fn = sampleBSJ)
  expect_equal(class(bsj.data), structure("GRanges", package = "GenomicRanges"))
  expect_equal(length(bsj.data), 10)
  expect_equal(GenomicRanges::start(bsj.data), c(8518655L, 107159492L, 42816025L, 2067967L, 49782626L, 11869601L, 10689961L, 74055408L, 68959805L, 127660597L))
  expect_equal(GenomicRanges::end(bsj.data), c(8521151L, 107159582L, 42818856L, 2068296L, 49798871L, 11869812L, 10690351L, 74055768L, 68960250L, 127661073L))
  expect_equal(GenomicRanges::strand(bsj.data), Rle(factor(c("-", "+", "-", "-", "+", "+", "+", "+", "+", "+"), levels = c("+", "-", "*"))))
  expect_equal(colnames(mcols(bsj.data)), c("bsID", "count", "donor", "acceptor"))
})

## Generate data that will be used for a couple of tests
l <- c(
  "8\t22100667\t-\t1\t150575816\t+\t-1\t0\t0\tD74RYQN1:328:C480EACXX:2:1103:19393:50637\t22100668\t72M104N29M\t150575817\t101M",
  "3\t4831313\t-\t19\t10391631\t+\t0\t0\t0\tD74RYQN1:328:C480EACXX:2:1103:19740:50650\t4831314\t84S17M\t10391632\t17S75M1379N9M284p63M1579N38M",
  "2\t108785689\t-\t2\t108785792\t-\t-1\t0\t0\tD74RYQN1:328:C480EACXX:2:1103:20773:50626\t108785690\t99M2S\t108785691\t101M",
  "3\t128806534\t-\t3\t128807589\t-\t-1\t0\t0\tD74RYQN1:328:C480EACXX:2:1103:2660:50949\t128806535\t56M952N45M\t128806537\t54M952N46M1S",
  "X\t154366799\t+\tX\t154366561\t+\t-1\t0\t0\tD74RYQN1:328:C480EACXX:2:1103:3745:50806\t154366607\t1S33M92N67M\t154366562\t78M92N23M",
  "2\t38750054\t+\t2\t38748082\t+\t0\t0\t7\tD74RYQN1:328:C480EACXX:2:1103:16101:50753\t38749667\t39M308N40M22S\t38748083\t79S22M1423p101M",
  "21\t34375547\t-\t21\t34385552\t-\t0\t0\t5\tD74RYQN1:328:C480EACXX:2:1103:2995:53019\t34375548\t72S25M-19p14M3860N87M\t34385480\t72M29S",
  "9\t128613401\t-\t9\t128618048\t-\t0\t0\t4\tD74RYQN1:328:C480EACXX:2:1103:8540:51242\t128613402\t35S66M2259p101M\t128618013\t35M66S",
  "6\t75702644\t-\t6\t75703073\t-\t2\t0\t5\tD74RYQN1:328:C480EACXX:2:1103:5721:51585\t75702645\t73S28M\t75702753\t101M146p73M28S",
  "19\t14480670\t+\t19\t14480421\t+\t0\t0\t3\tD74RYQN1:328:C480EACXX:2:1103:19372:51552\t14480409\t77M107N24M-135p4M107N77M20S\t14480422\t81S20M"
)
chim <- readChimFile(file = paste0(l, collapse = "\n"))

test_that(desc = "Check readChimFile",{
  expect_equal(class(chim), c("tbl_df", "tbl", "data.frame"))
  expect_equal(nrow(chim), 10)
  expect_equal(ncol(chim), 14)
  expect_equivalent(
    unlist(lapply(chim, class)),
    c("character", "integer", "character", "character",
      "integer", "character", "integer", "integer", "integer", "character",
      "integer", "character", "integer", "character"
    )
  )
  expect_equal(colnames(chim), paste0("X", 1:14))
})

# test_that(desc = "Check generateCandidateFilters",{
#   expect_equal(generateCandidateFilters(chim, backSpliceDist = 1e5, seqAcrossSj = T), rep(c(F,T), each=5))
#   expect_equal(generateCandidateFilters(chim, backSpliceDist = 1e3, seqAcrossSj = T), rep(c(F,T), c(8,2)))
#   expect_equal(generateCandidateFilters(chim, backSpliceDist = 1e2, seqAcrossSj = T), rep(F, 10))
# })
#
# test_that(desc = "Check swapStrands",{
#   circ.chim <- chim[generateCandidateFilters(chim, backSpliceDist = 1e5, seqAcrossSj = T),]
#   expect_equal(swapStrands(circ.chim)$X1, circ.chim$X1)
#   expect_equal(swapStrands(circ.chim)$X2, circ.chim$X5)
#   expect_false(isTRUE(all.equal(swapStrands(circ.chim)$X3, circ.chim$X3)))
# })

## Do we need to test getCandidateBackspliceSites now that we have tested all "sub-functions"?????
# test_that(desc = "Check getCandidateBackspliceSites",{
# })
