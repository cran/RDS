context("impute.degree: stratified imputation alignment (2.1)")

# Finding 2.1. impute.degree stratified the imputation by a trait, but the
# per-stratum loop indexed the na.omit'd degree/weight vectors (length n - #NA)
# with `by.var == i`, a length-n mask. After the first missing degree everything
# shifted, so each stratum pooled degrees from the wrong trait groups -- the
# stratum means were silently wrong, and stratified imputation is the entire
# purpose of the function. The fix aligns degrees, weights, and grouping through
# a single `obs <- !is.na(degree)` mask (`by.var[obs] == i`).
#
# Oracle: the function's own documented example. The correct VH-weighted stratum
# means are blue = 14.5052, red = 15.0692; the misaligned code produced 15.53 /
# 14.63 (and a later weight-misalignment variant gave 16.64 / 17.17) -- all far
# from the correct values, so the check has strong teeth.

test_that("stratum means align degrees, weights, and trait groups", {
  data(faux)
  na.idx <- c(1, 2, 30, 52, 81, 101, 108, 111)
  d <- faux
  d$network.size[na.idx] <- NA

  out <- impute.degree(d, trait.variable = "X", round.degree = FALSE)

  expect_false(any(is.na(out)) || any(is.nan(out)))

  grp      <- as.numeric(as.factor(faux$X))[na.idx]     # blue = 1, red = 2
  expected <- c(14.50516, 15.06924)[grp]                # correct per-group means
  expect_equal(unname(out[na.idx]), expected, tolerance = 1e-4)

  # observed rows are returned untouched (only the recruitment lower bound)
  keep <- setdiff(seq_len(nrow(faux)), na.idx)
  base <- pmax(get.net.size(faux), get.number.of.recruits(faux) + 1)
  expect_equal(out[keep], as.numeric(base)[keep])
})

test_that("a fully-missing trait level falls back to the overall mean (no NaN)", {
  # If a stratum has no observed degrees its weighted mean is NaN; those rows
  # must not be left NaN -- they fall through to the overall weighted mean.
  data(faux)
  d <- faux
  blue <- as.character(faux$X) == "blue"
  d$network.size[blue] <- NA                            # wipe an entire stratum

  out <- impute.degree(d, trait.variable = "X", round.degree = FALSE)

  expect_false(any(is.na(out)) || any(is.nan(out)))
  # every wiped-stratum row gets the single overall fallback mean
  expect_equal(length(unique(out[blue])), 1L)
  expect_equal(unique(out[blue]), 15.0292, tolerance = 1e-3)
})
