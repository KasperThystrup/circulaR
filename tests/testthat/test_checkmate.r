# Import
context(desc = "Check helper functions")
library(circulaR)
library(GenomicRanges)
library(IRanges)

l <- c(
  "1\t180084017\t-\t1\t180098986\t-\t2\t0\t1\tD74RYQN1:328:C480EACXX:2:1103:1924:52354\t180084018\t22S79M\t180096085\t100M1S2779p22M79S",
  "1\t32187953\t-\t1\t32188125\t-\t2\t0\t2\tD74RYQN1:328:C480EACXX:2:1103:5915:54032\t32187954\t82S19M\t32187952\t8S93M-2p82M19S",
  "1\t19201816\t+\t1\t19197669\t+\t2\t0\t1\tD74RYQN1:328:C480EACXX:2:1103:4863:55474\t19201732\t84M17S\t19197670\t84S17M338p25M491N76M",
  "1\t40069801\t-\t1\t40070857\t-\t0\t0\t1\tD74RYQN1:328:C480EACXX:2:1103:10312:62718\t40069802\t62S39M\t40069788\t87M284N14M299p41M323N21M39S",
  "1\t51349606\t+\t1\t51349314\t+\t0\t0\t1\tD74RYQN1:328:C480EACXX:2:1103:16154:66240\t51349559\t47M54S\t51349315\t47S54M121p101M",
  "1\t235111222\t+\t1\t235110951\t+\t1\t1\t5\tD74RYQN1:328:C480EACXX:2:1213:10685:66708\t235111179\t43M58S\t235110952\t43S58M122p101M"
)
df <- readChimFile(file = paste0(l, collapse = "\n"))

## Does function 'checkPEReads' behave as expected
test_that("PEreadsOK", {
  expect_equal(class(checkPEReads(df)), "logical")
  expect_true(length(checkPEReads(df)) == nrow(df))
  #expect_true(all(checkPEReads(df, endTol = 0) == c(T, F, T, F, T, F)))
  #expect_true(all(checkPEReads(df, endTol = 10) == c(T, T, T, F, T, T)))
  expect_true(all(checkPEReads(df, endTol = 100)))
})

