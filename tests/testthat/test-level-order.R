context("Estimator invariance to factor level order")

# Regression tests for the rds.I.weights level-order bug: the per-group count
# vector was indexed in the outcome factor's own level order while the
# stationary-distribution estimates are in sorted level order, so any
# non-alphabetical level ordering silently divided each weight by the wrong
# group's count.

test_that("RDS-I estimates are invariant to 2-level factor order", {
  data(faux)

  e.alpha <- RDS.I.estimates(faux, "X")@estimate

  faux$X1 <- factor(faux$X, levels = c("red", "blue"))
  e.rev <- RDS.I.estimates(faux, "X1")@estimate

  expect_equal(e.rev[names(e.alpha)], e.alpha)
})

test_that("RDS-I estimates are invariant to 3-level factor permutations", {
  data(faux)

  # Build a 3-level outcome deterministically from the data.
  set.seed(7)
  faux$Z <- ifelse(faux$X == "blue", "apple",
                   sample(c("banana", "cherry"), nrow(faux), replace = TRUE))

  orders <- list(
    c("apple", "banana", "cherry"),   # alphabetical (reference)
    c("banana", "cherry", "apple"),   # 3-cycle: catches inverse-permutation bugs
    c("cherry", "banana", "apple")    # reversal
  )

  ests <- lapply(orders, function(o) {
    dat <- faux
    dat$Z <- factor(dat$Z, levels = o)
    RDS.I.estimates(dat, "Z")@estimate
  })

  ref <- ests[[1]]
  expect_equal(ests[[2]][names(ref)], ref)
  expect_equal(ests[[3]][names(ref)], ref)
})

test_that("RDS-I (data smoothed) estimates are invariant to level order", {
  data(faux)

  e.alpha <- RDS.I.estimates(faux, "X", smoothed = TRUE)@estimate

  faux$X1 <- factor(faux$X, levels = c("red", "blue"))
  e.rev <- RDS.I.estimates(faux, "X1", smoothed = TRUE)@estimate

  expect_equal(e.rev[names(e.alpha)], e.alpha)
})

test_that("2-level HCG (small.fraction) matches RDS-I and is level-order invariant", {
  data(faux)

  e.rdsi <- RDS.I.estimates(faux, "X")@estimate
  e.hcg <- RDS.HCG.estimates(faux, "X", small.fraction = TRUE)@estimate

  # For a 2-level outcome both estimators compute the same closed-form
  # expression, but with one definitional difference: RDS-I's per-group
  # harmonic mean degree (get.h.hat) includes seeds, while HCG's dbar
  # excludes them (hcg-util.R, dg[!seed]). faux has a single seed, giving an
  # observed difference of ~5e-5; the gap scales with the number and
  # atypicality of seeds, so agreement is exact only in the no-seed limit.
  expect_equal(e.hcg[names(e.rdsi)], e.rdsi, tolerance = 1e-3)

  # HCG must also be invariant to the factor's level order.
  faux$X1 <- factor(faux$X, levels = c("red", "blue"))
  e.hcg.rev <- RDS.HCG.estimates(faux, "X1", small.fraction = TRUE)@estimate
  expect_equal(e.hcg.rev[names(e.hcg)], e.hcg)
})
