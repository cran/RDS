

library(testthat)
library(RDS)

context("differential-activity.R")

test_that("differential-activity",{
  data(faux)
  est <- differential.activity.estimates(faux,"X",weight.type="RDS-II")
  # ratio value is unchanged by the 4.12 label fix
  expect_equal(as.numeric(est[1]), 1.03342283930911, tolerance = 1e-6)
})

# 4.12 Bug A: the printed label used max(outcome) *after* recoding to 0/1, so it
# always said "value 1" instead of the category name. The category name is now
# carried in attr(.,"level") and printed.
test_that("label names the category, not the recoded 1 (4.12 bug A)",{
  data(faux)
  est <- differential.activity.estimates(faux,"X",weight.type="RDS-II")
  expect_true(is.character(attr(est,"level")))
  expect_equal(attr(est,"level"), "red")
  out <- paste(capture.output(print(est)), collapse = " ")
  expect_true(grepl("value red", out, fixed = TRUE))
  expect_false(grepl("value 1", out, fixed = TRUE))
})

# 4.12 Bug B: for 3+ level factors the "0" group is everyone-not-the-top-group,
# not the second-most-frequent group the selection logic implies. The label now
# honestly says "all others"; the coded group is the alphabetical max of the two
# most-frequent levels.
test_that("3-level factor labels the coded level and 'all others' (4.12 bug B)",{
  data(faux)
  d <- faux
  n <- nrow(d)
  # deterministic counts: aa (50%) > bb (30%) > cc (20%)  -> top two are aa,bb
  na <- floor(n * 0.5); nb <- floor(n * 0.3)
  g <- rep("cc", n); g[seq_len(na)] <- "aa"; g[na + seq_len(nb)] <- "bb"
  d$G <- factor(g)
  est <- differential.activity.estimates(d,"G",weight.type="RDS-II")
  expect_equal(attr(est,"level"), "bb")   # max(c("aa","bb"))
  out <- paste(capture.output(print(est)), collapse = " ")
  expect_true(grepl("value bb", out, fixed = TRUE))
  expect_true(grepl("all others", out, fixed = TRUE))
})
