context("measurement-error: optimism=FALSE means no optimism (1.3)")

# Finding 1.3. The visibility measurement-error log-likelihood (llcmpme/llnbme)
# appended the literal 1 to the parameter vector when optimism=FALSE, but the C
# likelihood reads that slot on the LOG scale (opt = exp(v)), so "optimism off"
# actually fit opt = exp(1) = e ~ 2.718 instead of the intended opt = 1. The
# imputation pdf (dmepdf) reads the same slot RAW and used opt = 1, so the
# fitted model and the model used to impute disagreed. optimism=FALSE is the
# default of impute.visibility_mle(), so this is the default path.
#
# Invariant: turning optimism off must equal fitting with the optimism
# parameter pinned to "no optimism" (a natural multiplier of 1). In the
# optimism=TRUE parametrization the C code exponentiates the slot, so the
# no-optimism value is the slot value 0 (exp(0) = 1) -- NOT 1 (exp(1) = e).
# Before the fix off == (slot 1); after the fix off == (slot 0).

test_that("optimism=FALSE equals the optimism parameter pinned to no-optimism", {
  set.seed(1)
  n <- 50
  network.size  <- sample(4:25, n, replace = TRUE)
  num.recruits  <- sample(0:3, n, replace = TRUE)
  recruit.times <- sort(runif(n))
  gmean <- length(network.size) / sum(1 / network.size)

  # Parameter layout of the `v` argument (before llcmpme/llnbme prepend the mean
  # and, with unit.scale < 0, estimate the scale):
  #   beta0, beta1(recruit-time coef), log(sd), log(scale)  [+ optimism slot]
  v4 <- c(0.1, 0.05, log(0.6), log(1.3))

  ll <- function(fn, v, optimism)
    fn(v, n = n, network.size = network.size, num.recruits = num.recruits,
       recruit.time = TRUE, recruit.times = recruit.times, max.coupons = 3,
       K = 200, unit.scale = -1, optimism = optimism, gmean = gmean)

  for (fn in list(RDS:::llcmpme, RDS:::llnbme)) {
    off <- ll(fn, v4,        optimism = FALSE)  # "no optimism"
    on0 <- ll(fn, c(v4, 0),  optimism = TRUE)   # opt slot 0 -> natural opt = exp(0) = 1
    on1 <- ll(fn, c(v4, 1),  optimism = TRUE)   # opt slot 1 -> natural opt = exp(1) = e

    # sanity: a finite, non-degenerate likelihood and the slot genuinely matters
    expect_true(is.finite(off) && off > -1e9)
    expect_false(isTRUE(all.equal(on0, on1)))

    # the fix: optimism off == optimism pinned to no-optimism (opt = 1)
    expect_equal(off, on0, tolerance = 1e-8)
  }
})
