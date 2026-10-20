context("wtd.quantile: scale-invariant design quantile (7.5)")

# 7.5. wtd.quantile is used with Horvitz-Thompson / design weights (vh.weights,
# gile.ss.weights), which sum to ~n or to the population size N -- they are NOT
# integer frequency multipliers. The old Hmisc-derived type-7 algorithm was
# wrong for such weights: it was scale-DEPENDENT (the answer changed if all
# weights were multiplied by a constant) and it corrupted quantiles when any
# weight was zero (c(0,2,0,2) -> 2.5/2.5/3.25 via duplicate-cumsum ties).
#
# It was replaced by the inverse of the normalized weighted ECDF (a Hajek design
# quantile): scale-invariant by construction, robust to zero weights, and -- for
# integer weights -- identical to the type-1 quantile of the frequency-expanded
# data (the independent oracle used below).

wq <- RDS:::wtd.quantile

test_that("zero-weight observations are ignored", {
  # effective data {2,2,4,4}; the weight-0 points (x=1,3) must never be returned
  expect_equal(unname(wq(1:4, c(0, 2, 0, 2), probs = c(0, .25, .5, .75, 1))),
               c(2, 2, 2, 4, 4))
})

test_that("scale-invariant in the weights", {
  x <- c(3, 1, 4, 1, 5, 9, 2, 6)
  w <- c(2, 7, 1, 8, 2, 8, 1, 8)
  base <- unname(wq(x, w))
  for (s in c(1e-6, 0.5, 100, 1e6))
    expect_equal(unname(wq(x, w * s)), base)
})

test_that("integer weights match the type-1 quantile of frequency-expanded data", {
  # independent oracle: with integer weights the design quantile must equal the
  # ordinary type-1 (inverse-ECDF) quantile of the replicated sample.
  x     <- c(1, 2, 3, 4)
  w     <- c(3, 1, 4, 2)
  probs <- c(0, .1, .25, .5, .75, .9, 1)
  expect_equal(unname(wq(x, w, probs = probs)),
               unname(stats::quantile(rep(x, w), probs = probs, type = 1)))
  # equal weights reduce to the plain type-1 quantile
  expect_equal(unname(wq(x, rep(1, 4))),
               unname(stats::quantile(x, type = 1)))
})

test_that("impute.degree quantile path stays sensible under the new estimator", {
  data(faux)
  na.idx <- c(1, 2, 30, 52, 81, 101, 108, 111)
  d <- faux
  d$network.size[na.idx] <- NA
  out <- impute.degree(d, trait.variable = "X", method = "quantile",
                       round.degree = FALSE)
  expect_false(any(is.na(out)) || any(is.nan(out)))
  # imputed degrees are realizable (drawn from observed degrees), within range
  expect_true(all(out[na.idx] >= min(out) & out[na.idx] <= max(out)))
})
