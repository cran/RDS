context("Bootstrap agreement across schemes (small sample fraction)")

# Cross-implementation sanity check across the three model-based resampling
# schemes -- Salganik, HCG, and Gile's SS. At a small sample fraction they
# share the same first-order recruiter->recruit transition model and each
# draws recruit degrees proportional to degree (size-biased population), so
# their bootstrap distributions should agree in center and be within a modest
# variance band. To isolate the resampling machinery from the choice of point
# estimator, all three evaluate the SAME estimator (RDS-II) on their replicates.
#
# They are not identical bootstraps (different seed composition, tree vs. chain
# structure, PPSWOR depletion that vanishes only in the limit), so this asserts
# approximate agreement, not equality. It guards against a regression that
# GROSSLY breaks the agreement of a scheme -- e.g. a group ordering/labeling
# bug, or applying the depletion regime at a small fraction. It is not tuned to
# catch subtle degree-resampling tilts (a ~12% SD shift stays inside these
# bands, since one scheme is typically the high-variance one regardless);
# calibration: the assertions PASS at n/N ~= 0.04 (sd max/min <= 1.15) and FAIL
# hard at n/N ~= 0.86 (mean gaps to 0.21, sd max/min to 3.7). Deliberately a
# "less favorable" case: 3 groups with strongly unequal degrees.

test_that("Salganik, HCG and Gile bootstraps agree at small sample fraction", {
  # The agreement bands are computed from three bootstrap schemes, one of which
  # (HCG) runs an optimizer and another (Gile) a C sampler, so the bootstrap
  # SDs/means depend on BLAS/LAPACK, compiler, and sample.kind. The bands have
  # margin but are not guaranteed reproducible across CRAN's platforms; kept as
  # a local/CI guard.
  skip_on_cran()
  data(faux)

  set.seed(202)
  n <- nrow(faux)
  grp <- ifelse(faux$X == "blue", 1L, sample(2:3, n, replace = TRUE))
  deg <- ceiling(c(26, 13, 6)[grp] * rexp(n, 4))#exp(rnorm(n, 0, 0.30)))
  deg[deg < 1] <- 1L
  #deg[grp == 1][5:30] <- 400
  N <- 10000L          # n/N ~= 0.04  (small sample fraction)
  B <- 300L

  bs.dist <- function(uncertainty) {
    set.seed(1)
    d <- faux
    d[[attr(d, "network.size.variable")]] <- deg
    d$grp <- factor(c("a", "b", "c")[grp], levels = c("a", "b", "c"))
    if(uncertainty == "Salganik"){
      #d[[attr(d, "network.size.variable")]] <- deg^2
    }
    iv <- suppressWarnings(RDS.bootstrap.intervals(
      d, "grp", weight.type = "RDS-II", uncertainty = uncertainty,
      N = N, number.of.bootstrap.samples = B, verbose = FALSE))
    bs <- attr(iv, "bsresult")
    if (is.list(bs) && !is.null(bs$bsests)) bs$bsests else bs   # reps x groups
  }

  S <- bs.dist("Salganik")
  H <- bs.dist("HCG")
  G <- bs.dist("Gile")

  expect_equal(ncol(S), 3L)
  expect_equal(ncol(H), 3L)
  expect_equal(ncol(G), 3L)

  for (g in 1:3) {
    mns <- c(mean(S[, g]), mean(H[, g]), mean(G[, g]))
    sds <- c(sd(S[, g]),   sd(H[, g]),   sd(G[, g]))
    # centers agree (all three consistent for RDS-II in the small-fraction limit)
    expect_lt(max(mns) - min(mns), 0.03)
    # spreads agree within a modest band (distinct schemes + MC noise)
    expect_lt(max(sds) / min(sds), 1.7)
  }
})
