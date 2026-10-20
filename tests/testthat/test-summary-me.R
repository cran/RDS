context("summary.me: AIC/BIC use parameter count, not n (2.3)")

# Finding 2.3. summary.me scaled the AIC/BIC penalties by `dyads` (= object$df =
# length(network.size) = n, the number of respondents) instead of by the
# parameter count p = length(coef). AIC was -2*ll + 2*n (should be 2*p) and BIC
# was -2*ll + log(p)*n (should be log(n)*p) -- both the multiplier and, for BIC,
# the log's argument were the wrong quantity. With n >> p the reported values
# were wildly inflated and useless for model comparison.
#
# Fit-free regression: a minimal fake "me" object exercises the formulas
# directly. It calls the `summary` generic (not RDS:::summary.me) to also cover
# the S3 registration from 2.6 item 7 -- summary.me/print.summary.me are now
# registered via @exportS3Method, so a user-level summary() dispatches here
# instead of falling back to summary.default.

test_that("AIC and BIC penalise the parameter count, with n = respondent count", {
  n <- 500L
  p <- 4L
  ll   <- -1234.5
  llnull <- -1500
  obj <- list(
    coef        = setNames(c(0.1, 0.2, 0.3, 0.4),
                           c("Recruitment Odds", "Error log-s.d.", "CMP skew", "Optimism")),
    covar       = diag(0.01, p),
    df          = n,             # summary.me reads this as `dyads` (the sample size)
    loglik      = ll,
    loglik.null = llnull,
    formula     = NULL, iterations = 10L, control = NULL
  )
  class(obj) <- "me"

  s <- summary(obj)                                # generic dispatch (2.6 item 7)
  expect_s3_class(s, "summary.me")                 # not summary.default

  expect_equal(s$aic, -2 * ll + 2 * p)             # 2*p, not 2*n
  expect_equal(s$bic, -2 * ll + log(n) * p)        # log(n)*p, not log(p)*n

  # teeth: the old (buggy) values were multiplied by n
  expect_false(isTRUE(all.equal(s$aic, -2 * ll + 2 * n)))
  expect_false(isTRUE(all.equal(s$bic, -2 * ll + log(p) * n)))
})
