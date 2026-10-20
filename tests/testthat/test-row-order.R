context("Estimator invariance to row order")

# Regression tests for the hcg.estimate seed-flag bug: the seed indicator was
# computed on the input-order vectors while everything it was combined with
# came from the frame sorted by recruitment time, so any data set whose rows
# were not already time-sorted attached seed flags to the wrong subjects.
# hcg.estimate now also breaks time ties by subject id, so results are fully
# deterministic under row permutation even when interview times are tied
# (fauxtime has 511 rows but only 27 distinct times).

permuted.fauxtime <- function(perm) {
  data(fauxtime, envir = environment())
  df <- as.data.frame(fauxtime)[perm, ]
  suppressWarnings(
    as.rds.data.frame(df,
                      id = "SER",
                      recruiter.id = "recruiter.id",
                      network.size = "NETWORK",
                      time = "DATEINTERVIEW",
                      population.size = get.population.size(fauxtime))
  )
}

test_that("HCG (small.fraction) is invariant to row order", {
  data(fauxtime)
  set.seed(11)
  perm <- sample(nrow(fauxtime))
  fx <- permuted.fauxtime(perm)

  # This path is closed-form (no RNG), so equality must be exact.
  e0 <- RDS.HCG.estimates(fauxtime, "var1", small.fraction = TRUE)@estimate
  e1 <- RDS.HCG.estimates(fx, "var1", small.fraction = TRUE)@estimate
  expect_equal(e0, e1)
})

test_that("wave-fallback recruitment times respect the recruitment tree", {
  # Regression test for get.recruitment.time(wave.fallback=TRUE), which
  # returned order(wave) -- the sorting permutation -- instead of per-subject
  # ranks. The two coincide only when rows happen to be pre-sorted by wave,
  # so any other row order silently corrupted the fallback times feeding
  # hcg.weights. These assertions hold for any monotone-in-wave scheme
  # (rank with any ties.method, or raw waves), so they don't pin down the
  # tie-breaking choice.
  data(faux)
  set.seed(5)
  perm <- sample(nrow(faux))
  df <- as.data.frame(faux)[perm, ]
  fx <- as.rds.data.frame(df, id = "id", recruiter.id = "recruiter.id",
                          network.size = "network.size")

  for (dat in list(faux, fx)) {
    tt <- get.recruitment.time(dat, wave.fallback = TRUE)
    wv <- get.wave(dat)

    # No subject may have a fallback time earlier than their recruiter's.
    ri <- match(get.rid(dat), get.id(dat))
    ok <- !is.na(ri)
    expect_equal(sum(tt[ok] < tt[ri[ok]]), 0)

    # Times must be ordered by wave: the bands for successive waves
    # must not overlap.
    hi <- tapply(tt, wv, max)
    lo <- tapply(tt, wv, min)
    ord <- order(as.numeric(names(hi)))
    expect_true(all(hi[ord][-length(ord)] < lo[ord][-1]))
  }
})

test_that("HCG (small.fraction) on wave-fallback data is invariant to row order", {
  data(faux)
  set.seed(5)
  perm <- sample(nrow(faux))
  df <- as.data.frame(faux)[perm, ]
  fx <- as.rds.data.frame(df, id = "id", recruiter.id = "recruiter.id",
                          network.size = "network.size")

  e0 <- RDS.HCG.estimates(faux, "X", small.fraction = TRUE)@estimate
  e1 <- RDS.HCG.estimates(fx, "X", small.fraction = TRUE)@estimate
  expect_equal(e0, e1[names(e0)])
})

test_that("HCG (full likelihood) estimate and weights are invariant to row order", {
  data(fauxtime)
  set.seed(11)
  perm <- sample(nrow(fauxtime))
  fx <- permuted.fauxtime(perm)

  # The full path consumes RNG (gile.ss.weights simulates inclusion
  # probabilities), so pin the seed identically for each call; given that,
  # results must be exactly invariant to the input row order.
  set.seed(42)
  e0 <- RDS.HCG.estimates(fauxtime, "var1", N = 3000)@estimate
  set.seed(42)
  e1 <- RDS.HCG.estimates(fx, "var1", N = 3000)@estimate
  expect_equal(e0, e1)

  set.seed(42)
  w0 <- hcg.weights(fauxtime, "var1", N = 3000)
  set.seed(42)
  w1 <- hcg.weights(fx, "var1", N = 3000)
  # w0 is in original fauxtime row order; w1 in permuted order.
  expect_equal(w0[perm], w1)
})
