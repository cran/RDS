context("impute.visibility_mle: all-NA recruitment times error, not loop (2.5)")

# Finding 2.5. When every recruitment time is NA, the moving-median fill loop
# never terminates: median() of all-NA neighbours is NA, so
# `while(any(is.na(recruit.times)))` can never clear and the call hangs. The fix
# stops early ("No valid recruitment times.") before the loop when no finite
# time exists. This test would HANG on the pre-fix code, so it also documents
# that the guard fires fast (before any model fit).

test_that("all-NA recruitment times raise an error instead of looping forever", {
  data(faux)
  d <- faux
  d$badtime <- NA_real_
  expect_error(
    suppressWarnings(impute.visibility_mle(d, recruit.time = "badtime",
                                           verbose = FALSE)),
    "No valid recruitment times")
})
