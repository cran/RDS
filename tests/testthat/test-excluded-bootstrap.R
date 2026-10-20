# NA/subset boundary handling in the bootstrap uncertainty methods
# (NA_SUBSET_HANDLING.md): factor outcomes carry excluded subjects through
# the engines as an internal group that is simulated like any other and
# removed (with renormalization) from every replicate estimate; continuous
# outcomes fall back per method (HCG masks natively; Salganik/Gile drop the
# rows with an explicit warning).

est.of <- function(x) tryCatch(x@estimate, error = function(e) x$estimate)

test_that("factor outcomes with NAs run through all four engines (T5)", {
  skip_on_cran()
  data(fauxmadrona, package="RDS")
  fm <- fauxmadrona
  fm$disease <- factor(ifelse(fm$disease == 1, "pos", "neg"))
  set.seed(3); fm$disease[sample(nrow(fm), 25)] <- NA
  cfgs <- list(list(unc="Salganik", wt="RDS-II"),
               list(unc="Treeboot", wt="RDS-II"),
               list(unc="Gile",     wt="Gile's SS"),
               list(unc="HCG",      wt="HCG"))
  for(cf in cfgs){
    set.seed(11)
    x <- suppressWarnings(RDS.bootstrap.intervals(fm, "disease",
          weight.type=cf$wt, uncertainty=cf$unc, N=1000,
          number.of.bootstrap.samples=10, verbose=FALSE))
    est <- est.of(x)
    expect_false(RDS:::.rds.excluded %in% names(est))
    expect_equal(names(est)[1:2], c("neg","pos"), info=cf$unc)
    expect_equal(sum(est[1:2]), 1, tolerance=1e-8, info=cf$unc)
    expect_true(all(is.finite(est)), info=cf$unc)
    if(cf$unc == "Gile"){
      bsr <- attr(x, "bsresult")
      expect_false(RDS:::.rds.excluded %in% colnames(bsr$bsests))
      expect_true(all(abs(rowSums(bsr$bsests) - 1) < 1e-8))
    }
  }
})

test_that("continuous outcome with NAs runs in every engine without dropping rows (T4)", {
  skip_on_cran()
  data(fauxmadrona, package="RDS")
  fm <- fauxmadrona
  set.seed(2); fm$Y <- rnorm(nrow(fm), 10, 2)
  fm$Y[sample(nrow(fm), 25)] <- NA
  ivl <- function(x) x$interval[1, 1:3]  # point, lower, upper
  no.drop.warning <- function(expr){
    wmsgs <- character()
    r <- withCallingHandlers(expr, warning = function(w){
      wmsgs <<- c(wmsgs, conditionMessage(w)); invokeRestart("muffleWarning")})
    expect_false(any(grepl("removed", wmsgs)))
    r
  }
  # Gile: this configuration crashed pre-rework (sample.int(ncrefs):
  # refs.to.trait corrupted by the legacy raw-value reset + removal)
  set.seed(4)
  xg <- no.drop.warning(RDS.bootstrap.intervals(fm, "Y", weight.type="Gile's SS",
        uncertainty="Gile", N=1000, number.of.bootstrap.samples=8, verbose=FALSE))
  expect_true(all(is.finite(ivl(xg))))
  # Salganik: excluded subjects ride as the internal cut class; the bootstrap
  # point equals the point estimator exactly (same boundary math)
  set.seed(4)
  xs <- no.drop.warning(RDS.bootstrap.intervals(fm, "Y", weight.type="RDS-II",
        uncertainty="Salganik", N=1000, number.of.bootstrap.samples=8, verbose=FALSE))
  expect_true(all(is.finite(ivl(xs))))
  pe <- suppressWarnings(RDS.II.estimates(fm, "Y"))
  expect_equal(unname(ivl(xs)[1]), unname(pe@estimate), tolerance = 1e-8)
  # HCG masks continuous NAs natively
  set.seed(4)
  xh <- suppressWarnings(RDS.bootstrap.intervals(fm, "Y", weight.type="HCG",
        uncertainty="HCG", N=1000, number.of.bootstrap.samples=6, verbose=FALSE))
  expect_true(all(is.finite(ivl(xh))))
})

test_that("subset= flows through a bootstrap without orphaning recruits", {
  skip_on_cran()
  data(fauxmadrona, package="RDS")
  fm <- fauxmadrona
  set.seed(7)
  fm$grp <- factor(sample(c("a","b","c"), nrow(fm), TRUE, prob=c(.5,.3,.2)))
  set.seed(21)
  x <- suppressWarnings(RDS.bootstrap.intervals(fm, "grp", weight.type="RDS-II",
        uncertainty="Salganik", N=1000, number.of.bootstrap.samples=8,
        verbose=FALSE, subset = disease==0))
  est <- est.of(x)
  expect_equal(names(est)[1:3], c("a","b","c"))
  expect_equal(sum(est[1:3]), 1, tolerance=1e-8)
  expect_false(RDS:::.rds.excluded %in% names(est))
})

test_that("HCG.bootstrap replicates expose NA, never the internal label", {
  skip_on_cran()
  data(fauxtime, package="RDS")
  ft <- fauxtime
  ft$var1 <- factor(ifelse(ft$var1 == 1, "y", "n"))
  set.seed(5); ft$var1[sample(nrow(ft), 40)] <- NA
  set.seed(2)
  reps <- suppressWarnings(RDS:::HCG.bootstrap(ft, "var1",
          number.of.bootstrap.samples=3, N=3000, small.fraction=TRUE,
          verbose=FALSE, fun=function(x) x$var1))
  expect_false(RDS:::.rds.excluded %in% unlist(lapply(reps, levels)))
  expect_false(RDS:::.rds.excluded %in% unlist(lapply(reps, as.character)))
  expect_true(any(sapply(reps, anyNA)))
})

test_that("bootstrap.contingency.test runs with NAs; empty-cell case intact (T6)", {
  skip_on_cran()
  data(fauxtime, package="RDS")
  f4 <- fauxtime
  set.seed(4)
  f4$v2 <- ifelse(runif(nrow(f4)) < .5, "x", "y")
  f4$var1[sample(nrow(f4), 30)] <- NA
  set.seed(9)
  ct <- suppressWarnings(bootstrap.contingency.test(f4, "var1", "v2",
        number.of.bootstrap.samples=12, verbose=FALSE))
  expect_true(is.finite(ct$p.value))
  # the structurally-empty-cell regression (5.6) still passes
  f5 <- fauxtime
  f5$v2 <- ifelse(seq_len(nrow(f5)) %% 2 == 0, "x", "y")
  f5$v2[f5$var1 == 1] <- "x"
  set.seed(1)
  ct2 <- suppressWarnings(bootstrap.contingency.test(f5, "var1", "v2",
         number.of.bootstrap.samples=12, verbose=FALSE))
  expect_true(is.finite(ct2$p.value))
})

test_that("bootstrap.incidence runs with missing HIV/recency (T10)", {
  skip_on_cran()
  data(faux, package="RDS")
  f <- faux
  set.seed(6)
  f$recent <- rbinom(nrow(f), 1, .3)
  f$hiv <- rbinom(nrow(f), 1, .4)
  f$recent[f$hiv == 0] <- NA
  f$hiv[c(1, 6, 10, 21)] <- NA
  set.seed(10)
  r <- suppressWarnings(bootstrap.incidence(f, "recent", "hiv",
        weight.type="RDS-II", number.of.bootstrap.samples=12, verbose=FALSE))
  expect_true(is.finite(r$Incidence))
  expect_true(is.finite(r$SE))
})

test_that("hcg.replicate.weights zeroes excluded rows in all replicates (T10)", {
  skip_on_cran()
  data(fauxtime, package="RDS")
  ft <- fauxtime
  ft$var1 <- factor(ifelse(ft$var1 == 1, "y", "n"))
  set.seed(5); ft$var1[sample(nrow(ft), 30)] <- NA
  set.seed(13)
  rw <- suppressWarnings(hcg.replicate.weights(ft, "var1", 5,
        include.sample.weights=TRUE, N=3000, small.fraction=TRUE))
  expect_equal(ncol(rw), 6)
  ex <- is.na(ft$var1)
  expect_true(all(rw[ex, 1] == 0))
  expect_true(all(as.matrix(rw[ex, -1]) == 0))
  # non-seed observed subjects keep positive sample weights
  seeds <- get.rid(ft) == get.seed.rid(ft)
  expect_true(all(rw[!ex & !seeds, 1] > 0))
})
