context("measurement-error: CMP fit/impute scale consistency (1.4)")

# Finding 1.4. vtrans/ltrans apply a +1 offset to the scale slot and do NOT
# branch on unit.model: the natural scale reported in fit$coef is exp(l)+1 while
# the fitted C likelihood uses scale = exp(v[4]) = exp(l) = natural - 1. The
# nbinom imputation pdf (gnbmepdfC) undoes the offset with shape = dmean/(v[4]-1)
# and stays consistent with its fit, but the CMP pdf (gcmpmepdfC) read the slot
# raw (scale = v[4] = natural), so it imputed at a CMP scale one larger than the
# model was fit at. CMP is the default unit.model of impute.visibility_mle().
#
# Invariant: at a FIXED scale the fit and the pdf both read *cmpscale raw, so the
# fixed-scale pdf is an unambiguous reference for "internal CMP scale s". The
# estimated-scale pdf reports the natural scale (internal + 1), so evaluating it
# at natural scale s must equal the fixed-scale pdf at s - 1. This holds only
# when the estimated branch subtracts 1 (the fix); before the fix it used the
# raw natural scale and the two disagreed.

test_that("CMP estimated-scale pdf at natural s equals fixed-scale pdf at s-1", {
  set.seed(3)
  n <- 40
  network.size  <- sample(4:25, n, replace = TRUE)
  num.recruits  <- sample(0:3, n, replace = TRUE)
  recruit.times <- sort(runif(n))
  gmean <- length(network.size) / sum(1 / network.size)

  beta0 <- 0.1; beta1 <- 0.05; sd <- 0.6; opt <- 1.2
  s <- 1.7                     # natural CMP scale (reported in fit$coef)

  pdf <- function(coef, unit.scale)
    RDS:::dmepdf(coef, network.size = network.size, num.recruits = num.recruits,
                 recruit.time = TRUE, recruit.times = recruit.times,
                 max.coupons = 3, K = 200, unit.scale = unit.scale,
                 unit.model = "cmp", optimism = TRUE, gmean = gmean)

  # estimated-scale coef layout: beta0, beta1, sd, scale(natural), opt
  p_est   <- pdf(c(beta0, beta1, sd, s, opt), unit.scale = NULL)
  # fixed-scale coef layout drops the scale slot; feed s-1 as the fixed scale
  p_fixed <- pdf(c(beta0, beta1, sd, opt),    unit.scale = s - 1)

  expect_true(all(is.finite(p_est)) && sum(p_est) > 0)
  expect_equal(p_est, p_fixed, tolerance = 1e-10)

  # guard: the scale genuinely matters -- the buggy behaviour (no -1, i.e. the
  # estimated pdf at natural s) must NOT match the fixed pdf at s-1.
  p_fixed_raw <- pdf(c(beta0, beta1, sd, opt), unit.scale = s)   # internal scale s
  expect_false(isTRUE(all.equal(p_fixed, p_fixed_raw)))
})

test_that("nbinom fit/impute scale handling is unchanged by the CMP fix", {
  # gnbmepdfC already subtracts 1; the CMP edit must not have touched it. Pin a
  # snapshot of the nbinom imputation pdf so an accidental nbinom change trips.
  set.seed(3)
  n <- 40
  network.size  <- sample(4:25, n, replace = TRUE)
  num.recruits  <- sample(0:3, n, replace = TRUE)
  recruit.times <- sort(runif(n))
  gmean <- length(network.size) / sum(1 / network.size)

  p_nb <- RDS:::dmepdf(c(0.1, 0.05, 0.6, 1.7, 1.2),
                       network.size = network.size, num.recruits = num.recruits,
                       recruit.time = TRUE, recruit.times = recruit.times,
                       max.coupons = 3, K = 200, unit.scale = NULL,
                       unit.model = "nbinom", optimism = TRUE, gmean = gmean)
  expect_true(all(is.finite(p_nb)) && sum(p_nb) > 0)
})
