context("Gile bootstrap: C and R paths")

# The Gile bootstrap has two implementations of the replicate generator:
# bsC.c (useC=TRUE, the default for categorical outcomes) and the R bsfn
# (useC=FALSE, and always used for continuous outcomes). Historically the
# two disagreed structurally: the R path's offby<0 fallback fired whenever
# ANY class had zero deficit (sum(dif==0) for sum(dif)==0), and both paths
# mishandled zero-count classes in the fallback (R spun in a rejection
# loop; C could decrement a class to -1). These tests pin the repaired
# behavior:
#  - the two paths must agree in distribution (parity guard), and
#  - heavy rebalancing (population size barely above sample size, so
#    round(n*props2) creates zero classes and offby churn) must complete
#    cleanly through both paths.
# Note: offby is only ever a rounding residue, so the historical bugs'
# footprint on tame data is small; the sharp regression test for the C
# clamp bug lives in test-gile-ss.R. This file guards path parity.

gile.bsests <- function(dat, var, N, useC, B = 60) {
  set.seed(42)
  iv <- suppressWarnings(RDS.bootstrap.intervals(
    dat, var, weight.type = "Gile's SS", uncertainty = "Gile",
    N = N, number.of.bootstrap.samples = B, useC = useC, verbose = FALSE))
  attr(iv, "bsresult")$bsests[, 1]
}

test_that("C and R bootstrap paths agree in distribution", {
  data(faux)   # 'X' is categorical, so useC is honored
  bC <- gile.bsests(faux, "X", N = 1000, useC = TRUE, B = 200)
  bR <- gile.bsests(faux, "X", N = 1000, useC = FALSE, B = 200)

  expect_length(bC, 200)
  expect_length(bR, 200)

  # Same distribution: KS test (deterministic under the pinned seeds;
  # measured p ~ 0.7) and means within 3 combined standard errors.
  ks <- suppressWarnings(stats::ks.test(bC, bR))
  expect_gt(ks$p.value, 0.001)
  se.diff <- sqrt(stats::var(bC) + stats::var(bR)) / sqrt(length(bC))
  expect_lt(abs(mean(bC) - mean(bR)), 4 * se.diff)
})

test_that("bootstrap is invariant to the seed placeholder's name", {
  # Regression: refs.to.trait was built with tapply(...)[-1], which assumed
  # the seed pseudo-recruiter id sorts first among recruiter ids. With
  # "seed"-style placeholders (all four shipped datasets) it sorts last, so
  # [-1] silently deleted the lexicographically-first *real* recruiter's
  # referrals from the tie matrix. The placeholder needs no positional
  # removal at all: match()/!is.na() filters it by identity.
  data(faux)
  mk <- function(placeholder) {
    df <- as.data.frame(faux)
    rid <- as.character(df$recruiter.id)
    rid[rid == "seed"] <- placeholder
    df$recruiter.id <- rid
    as.rds.data.frame(df, id = "id", recruiter.id = "recruiter.id",
                      network.size = "network.size", population.size = 1000)
  }
  bs <- function(dat) {
    set.seed(42)
    iv <- suppressWarnings(RDS.bootstrap.intervals(
      dat, "X", weight.type = "Gile's SS", uncertainty = "Gile",
      N = 1000, number.of.bootstrap.samples = 40, useC = TRUE,
      verbose = FALSE))
    attr(iv, "bsresult")$bsests
  }
  # "seed" sorts after numeric ids; "0" sorts before them. Same tree,
  # so results must be identical.
  expect_identical(bs(mk("seed")), bs(mk("0")))
})

# Finding 7.1. When a categorical trait has missing values, they are recoded to
# an "NA" factor level and refs.to.trait is built with one column per sorted
# level (including "NA"). The missing-value removal block then subset the matrix
# ROWS only, leaving the "NA" column in place. Downstream, 1:g column indexing
# (CC/tiemtx) mislabels a real group whenever "NA" does not sort last, and the
# row-sum coupon count (nrefs) always includes the stray column -- so results
# depended on the arbitrary sort position of "NA" relative to the group labels.
# Fixed by also column-subsetting: refs.to.trait[remvalues, match(new, old)].
test_that("Gile bootstrap is invariant to trait label under missing outcomes (7.1)", {
  data(fauxmadrona)
  base <- fauxmadrona
  set.seed(42)
  na.idx <- sample(nrow(base), 30)          # precomputed: data build consumes no RNG
  mk <- function(pos) {
    o <- ifelse(base$disease == 1, pos, "aa")
    o[na.idx] <- NA
    d <- base; d$trait <- o; d
  }
  boot_se <- function(pos) {
    d <- mk(pos)
    set.seed(7)                             # identical RNG for every labeling
    iv <- suppressWarnings(RDS.bootstrap.intervals(
      d, "trait", weight.type = "RDS-II", uncertainty = "Gile",
      N = 1000, number.of.bootstrap.samples = 40,
      number.ss.samples.per.iteration = 200, verbose = FALSE))
    unname(attr(iv, "bsresult")$se_estimate)
  }
  # Same partition; only the second group's label differs, changing where "NA"
  # sorts ("aa,bb,NA" -> last; "aa,NA,zz" -> middle). Results must match.
  se.bb <- boot_se("bb")
  se.zz <- boot_se("zz")
  expect_equal(se.bb, se.zz, tolerance = 1e-8)
})

test_that("fast=FALSE (per-replicate mapping re-estimation) works", {
  # Regression: fast=FALSE set data.mapping to NULL, but pis was
  # unconditionally interpolated from it before the useC split, so the
  # documented slow mode always errored ("need at least two non-NA values
  # to interpolate"). It now interpolates from the original-sample mapping
  # and forces useC=FALSE (the C sampler cannot re-estimate the mapping
  # per replicate).
  data(faux)
  set.seed(42)
  iv <- suppressWarnings(RDS.bootstrap.intervals(
    faux, "X", weight.type = "Gile's SS", uncertainty = "Gile",
    N = 1000, number.of.bootstrap.samples = 10, fast = FALSE,
    useC = TRUE,   # must be silently rerouted, not crash
    verbose = FALSE))
  se <- matrix(iv$interval, ncol = 6)[, 5]
  expect_true(all(is.finite(se)))
  expect_true(all(se > 0))
})

test_that("heavy rebalancing (N barely above n) completes through both paths", {
  data(fauxmadrona)
  fm <- fauxmadrona
  fm$disease <- as.factor(fm$disease)  # categorical so useC routes

  # N = 520 vs n = 500: round(n*props2) zeroes rare (degree x disease)
  # classes and offby rebalancing runs hot, exercising the <=1 guards and
  # the last-resort branches. Pre-fix, the R path could spin in its
  # rejection loop and the C path could drive a class count negative.
  for (useC in c(TRUE, FALSE)) {
    b <- gile.bsests(fm, "disease", N = 520, useC = useC, B = 25)
    expect_length(b, 25)
    expect_true(all(is.finite(b)))
    expect_true(all(b >= 0 & b <= 1))
  }
})
