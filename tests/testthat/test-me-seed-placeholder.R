context("impute.visibility_mle: seed identification is placeholder-invariant (2.4)")

# Finding 2.4. `is.seed` hardcoded the literal "seed" recruiter marker, but the
# package allows any seed placeholder (get.seed.rid, e.g. "0"). With
# include.tree=TRUE a non-"seed" placeholder made `is.seed` all-FALSE, so seeds
# were treated as non-seeds and their augmented visibility floor
# (nr + !is.seed) was inflated by one. The fix identifies seeds via
# get.seed.rid.
#
# Invariance test: relabelling the seed placeholder must not change the
# imputation. The built-in datasets never activate the include.tree floor
# (network.size >= nr+1 everywhere), so we lower the seeds' reported sizes to
# make the floor active -- otherwise the bug is unobservable. With the floor
# active the placeholder affects ~100+ imputed rows, so reverting the fix
# (== "seed") breaks this exact-equality check. Invariance holds for any
# K/maxit because, once seeds are identified correctly, the two relabelled
# datasets produce byte-identical inputs and hence identical fits.

test_that("imputation does not depend on the seed placeholder label", {
  data(faux)
  nsv <- attr(faux, "network.size.variable")

  d_seed <- faux
  d_seed[[nsv]][get.wave(faux) == 0] <- 1L        # activate the seed floor (ns < nr)

  # identical data, seed placeholder relabelled "seed" -> "0"
  df <- as.data.frame(d_seed)
  df$recruiter.id <- as.character(df$recruiter.id)
  df$recruiter.id[df$recruiter.id == as.character(get.seed.rid(d_seed))] <- "0"
  d_zero <- as.rds.data.frame(df, id = "id", recruiter.id = "recruiter.id",
                              network.size = "network.size")
  expect_equal(get.seed.rid(d_zero), "0")          # sanity: placeholder really changed

  imp <- function(x) {
    set.seed(1)
    suppressWarnings(impute.visibility_mle(
      x, include.tree = TRUE, type.impute = "mode",
      K = 80, maxit = 25, verbose = FALSE))
  }
  expect_equal(imp(d_seed), imp(d_zero))
})
