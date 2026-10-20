context("bootstrap-intervals: CI assembly (5.3 refactor guard)")

# Guard tests for the planned refactor of the copy-pasted `switch(ci.type,...)`
# CI-assembly block in RDS.bootstrap.intervals (finding 5.3). These pin the
# CURRENT, correct behaviour so the refactor cannot silently change it:
#   * exact snapshots of the assembled $interval for every uncertainty x
#     ci.type (behaviour-preserving guard),
#   * relabel-invariance of the categorical CIs (the core correctness the
#     inverse-`match` reorder must preserve),
#   * structural sanity (point column == the observed estimate).
# All four uncertainty branches are exercised with a common RDS-II estimator:
# the ci.type assembly branch is selected by `uncertainty`, independent of
# weight.type, so this covers the refactored code cheaply. The snapshots assume
# the refactor touches only the post-sampling assembly (it must not change the
# bootstrap RNG stream).

make_data <- function() {
  data(faux, envir = environment())
  set.seed(101)
  part <- sample(1:3, nrow(faux), replace = TRUE)
  deg  <- c(35L, 18L, 6L)[part] + rpois(nrow(faux), 2)
  catd <- faux
  catd[[attr(catd, "network.size.variable")]] <- deg
  catd$grp <- factor(c("a", "b", "c")[part], levels = c("a", "b", "c"))
  ctsd <- faux
  ctsd[[attr(ctsd, "network.size.variable")]] <- deg
  ctsd$y <- round(deg + rnorm(nrow(faux), 0, 3), 3)   # keep RNG order stable
  list(cat = catd, cts = ctsd, part = part, deg = deg)
}

bsi <- function(d, var, unc, ci) {
  set.seed(1)
  suppressWarnings(RDS.bootstrap.intervals(
    d, var, weight.type = "RDS-II", uncertainty = unc,
    N = 1000, number.of.bootstrap.samples = 100, ci.type = ci, verbose = FALSE))
}

# Snapshots recorded from the pre-refactor code (RDS-II estimator, seed 1,
# B = 100). Treeboot + continuous is a separate crash (finding 5.12) and is
# intentionally omitted here.
baseline <- list(
  "Salganik.t"         = c(0.14411, 0.24299, 0.6129, 0.11339, 0.19763, 0.55528, 0.17482, 0.28835, 0.67052, 1.2665, 1.85214, 2.31675, 0.01567, 0.02314, 0.0294, 140, 125, 124),
  "HCG.t"              = c(0.14411, 0.24299, 0.6129, 0.12685, 0.20369, 0.5667, 0.16136, 0.28229, 0.65911, 0.39964, 1.39033, 1.48982, 0.0088, 0.02005, 0.02357, 140, 125, 124),
  "Treeboot.t"         = c(0.14411, 0.24299, 0.6129, 0.09406, 0.17069, 0.51946, 0.19415, 0.3153, 0.70634, 3.36215, 4.70569, 6.09294, 0.02553, 0.03689, 0.04767, 140, 125, 124),
  "Gile.t"             = c(0.14411, 0.24299, 0.6129, 0.17584, 0.22678, 0.50453, 0.21708, 0.29093, 0.58484, 0.57078, 0.92604, 1.12531, 0.01052, 0.01637, 0.02049, 140, 125, 124),
  "Salganik.quantile"  = c(0.14411, 0.24299, 0.6129, 0.112, 0.20261, 0.56068, 0.17683, 0.28333, 0.66755, 1.41053, 1.4665, 1.99258, 0.01654, 0.02059, 0.02726, 140, 125, 124),
  "HCG.quantile"       = c(0.14411, 0.24299, 0.6129, 0.13042, 0.21093, 0.56417, 0.1609, 0.28789, 0.65115, 0.31178, 1.33294, 1.31975, 0.00778, 0.01963, 0.02219, 140, 125, 124),
  "Treeboot.quantile"  = c(0.14411, 0.24299, 0.6129, 0.10018, 0.18008, 0.51158, 0.19565, 0.31937, 0.69622, 3.05878, 4.36565, 5.9482, 0.02436, 0.03553, 0.0471, 140, 125, 124),
  "Gile.quantile"      = c(0.14411, 0.24299, 0.6129, 0.15928, 0.21931, 0.52999, 0.20032, 0.28253, 0.61094, 0.56537, 0.89937, 1.14331, 0.01047, 0.01613, 0.02065, 140, 125, 124),
  "Salganik.proportion"= c(0.14411, 0.24299, 0.6129, 0.112, 0.20261, 0.56068, 0.17683, 0.28333, 0.66755, 1.41053, 1.4665, 1.99258, 0.01654, 0.02059, 0.02726, 140, 125, 124),
  "HCG.proportion"     = c(0.14411, 0.24299, 0.6129, 0.13042, 0.21093, 0.56417, 0.1609, 0.28789, 0.65115, 0.31178, 1.33294, 1.31975, 0.00778, 0.01963, 0.02219, 140, 125, 124),
  "Treeboot.proportion"= c(0.14411, 0.24299, 0.6129, 0.10018, 0.18008, 0.51158, 0.19565, 0.31937, 0.69622, 3.05878, 4.36565, 5.9482, 0.02436, 0.03553, 0.0471, 140, 125, 124),
  "Gile.proportion"    = c(0.14411, 0.24299, 0.6129, 0.15928, 0.21931, 0.52999, 0.20032, 0.28253, 0.61094, 0.56537, 0.89937, 1.14331, 0.01047, 0.01613, 0.02065, 140, 125, 124),
  "Salganik.cts.t"     = c(14.56013, 13.46244, 15.65783, 1.68739, 0.56006, 389),
  # HCG/Gile continuous CIs re-recorded after the 7.2 deg.out.map fix (`<=`->`<`):
  # the imputation pools were previously depleted, so the bootstrap SD/CI moved.
  # Point estimate (col 1) and n (col 6) are unchanged; Salganik.cts.t above is
  # unaffected (Salganik does not use deg.out.map).
  "HCG.cts.t"          = c(14.56013, 13.62911, 15.49116, 1.21387, 0.47502, 389),
  "Gile.cts.t"         = c(16.36395, 15.40022, 17.32768, 1.30065, 0.49171, 389),
  "Gile.cts.quantile"  = c(16.36395, 14.81641, 16.57976, 1.08859, 0.44984, 389)
)

test_that("CI assembly is unchanged (snapshot guard)", {
  # Exact bootstrap-interval snapshots to 1e-3. The HCG branch runs an optimizer
  # and the Gile branch a C sampler, so these values depend on BLAS/LAPACK,
  # compiler, and sample.kind and are not reproducible across platforms. Kept as
  # a local/CI refactor guard; skipped on CRAN's diverse test farm.
  skip_on_cran()
  dat <- make_data()
  specs <- list(
    Salganik.t = c("cat","grp","Salganik","t"),        HCG.t = c("cat","grp","HCG","t"),
    Treeboot.t = c("cat","grp","Treeboot","t"),        Gile.t = c("cat","grp","Gile","t"),
    Salganik.quantile = c("cat","grp","Salganik","quantile"), HCG.quantile = c("cat","grp","HCG","quantile"),
    Treeboot.quantile = c("cat","grp","Treeboot","quantile"), Gile.quantile = c("cat","grp","Gile","quantile"),
    Salganik.proportion = c("cat","grp","Salganik","proportion"), HCG.proportion = c("cat","grp","HCG","proportion"),
    Treeboot.proportion = c("cat","grp","Treeboot","proportion"), Gile.proportion = c("cat","grp","Gile","proportion"),
    Salganik.cts.t = c("cts","y","Salganik","t"),      HCG.cts.t = c("cts","y","HCG","t"),
    Gile.cts.t = c("cts","y","Gile","t"),              Gile.cts.quantile = c("cts","y","Gile","quantile")
  )
  for (tag in names(specs)) {
    s  <- specs[[tag]]
    iv <- bsi(dat[[s[1]]], s[2], s[3], s[4])
    expect_lt(max(abs(as.numeric(iv$interval) - baseline[[tag]])), 1e-3)  # info: tag
  }
})

test_that("categorical bootstrap CIs are invariant to group relabeling", {
  # Restricted to Salganik and Treeboot: these resample in a way whose RNG
  # stream is invariant to group order (Treeboot resamples the tree by row;
  # Salganik draws groups via R's prob-weighted sample.int, which sorts
  # internally), so with a fixed seed the whole CI is exactly relabel-invariant
  # -- a clean check that the ci.type assembly attaches each CI to the right
  # group. HCG/Gile use order-dependent PPSWOR depletion, so a fixed seed gives
  # MC-different (still correctly labeled) draws; their labeling is guarded by
  # the exact snapshots above instead.
  data(faux)
  set.seed(101)
  part <- sample(1:3, nrow(faux), replace = TRUE)
  deg  <- c(35L, 18L, 6L)[part] + rpois(nrow(faux), 2)
  build <- function(labels) {
    d <- faux
    d[[attr(d, "network.size.variable")]] <- deg
    d$grp <- factor(labels[part], levels = sort(unique(labels)))
    d
  }
  ci_matrix <- function(d, unc, ci) {
    set.seed(1)
    iv <- suppressWarnings(RDS.bootstrap.intervals(
      d, "grp", weight.type = "RDS-II", uncertainty = unc,
      N = 1000, number.of.bootstrap.samples = 100, ci.type = ci, verbose = FALSE))
    m <- matrix(iv$interval, ncol = 6)[, 1:3]         # point, lower, upper
    rownames(m) <- names(iv$estimate)
    m[order(rownames(m)), ]
  }
  remap <- c(m = "a", p = "b", c = "c")               # partition 1->m,2->p,3->c
  for (unc in c("Salganik", "Treeboot")) {
    A <- ci_matrix(build(c("a", "b", "c")), unc, "quantile")
    B <- ci_matrix(build(c("m", "p", "c")), unc, "quantile")
    rownames(B) <- remap[rownames(B)]
    B <- B[c("a", "b", "c"), ]
    expect_equal(unname(B), unname(A), tolerance = 1e-4)   # info: unc
  }
})

test_that("point-estimate column equals the observed RDS-II estimate", {
  dat <- make_data()
  obs <- RDS.II.estimates(dat$cat, "grp")@estimate
  for (unc in c("Salganik", "HCG", "Treeboot", "Gile")) {
    iv <- bsi(dat$cat, "grp", unc, "t")
    expect_equal(unname(iv$estimate), unname(obs), tolerance = 1e-8)
    pt <- matrix(iv$interval, ncol = 6)[, 1]
    expect_equal(unname(pt), unname(obs), tolerance = 1e-8)
  }
})

test_that("continuous outcome works with quantile/pivotal ci.type (5.3)", {
  # Before the ci.type refactor these crashed ("object 'outclasses' not found"
  # for Salganik; 0-column bsests for HCG). Treeboot + continuous is excluded:
  # it remains broken upstream in treeboot.bootstrap.estimates (finding 5.12).
  dat <- make_data()
  for (unc in c("Salganik", "HCG", "Gile")) for (ci in c("quantile", "pivotal")) {
    iv <- bsi(dat$cts, "y", unc, ci)
    expect_true(all(is.finite(iv$interval)))
    expect_equal(nrow(matrix(iv$interval, ncol = 6)), 1L)
  }
})
