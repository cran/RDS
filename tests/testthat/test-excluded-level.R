# NA/subset boundary handling (NA_SUBSET_HANDLING.md): missing outcomes and
# subset-excluded subjects are retained internally as their own factor level
# so the recruitment tree stays intact through the weighting machinery; the
# internal group receives weight 0 and is removed at the output boundary.
# Weighted estimates built from compute.weights() must equal those from the
# RDS.*.estimates() entry points for the same data/subset.

ht.prop <- function(w, out, keep){ a <- tapply(w[keep], out[keep], sum); a / sum(a) }
ht.mean <- function(w, out, keep) sum(w[keep] * out[keep]) / sum(w[keep])

test_that("excluded-level helper primitives", {
  ex <- RDS:::.rds.excluded
  x <- factor(c("a","b",NA,"a"))
  m <- RDS:::.excluded.mask(x)
  expect_equal(m, c(FALSE,FALSE,TRUE,FALSE))
  # minting is a no-op when nothing is excluded (bit-compatibility)
  y <- factor(c("a","b"))
  expect_identical(RDS:::.mint.excluded(y, c(FALSE,FALSE)), y)
  z <- RDS:::.mint.excluded(x, m)
  expect_equal(levels(z), c("a","b",ex))
  expect_equal(as.character(z)[3], ex)
  # existing level order is preserved (numeric-code factors stay numeric-ordered)
  num <- factor(c(1,2,10,NA))
  zn <- RDS:::.mint.excluded(num, is.na(num))
  expect_equal(levels(zn), c("1","2","10",ex))
  # .sort.levels: observed alphabetized, internal placed last
  f <- factor(c("b","a",ex), levels=c("b",ex,"a"))
  expect_equal(levels(RDS:::.sort.levels(f)), c("a","b",ex))
  # .drop.excluded matches by name and renormalizes
  est <- c(a=.2, b=.3); est[ex] <- .5
  expect_equal(RDS:::.drop.excluded(est), c(a=.4, b=.6))
  # .cut2.excl computes breaks from retained values only
  v <- c(1:10, 1000)
  cf <- RDS:::.cut2.excl(v, 2, c(rep(FALSE,10), TRUE))
  expect_equal(as.character(cf)[11], ex)
  expect_equal(nlevels(cf), 3)
  expect_equal(RDS:::.cut2.excl(1:10, 2, rep(FALSE,10)), RDS:::.cut2(1:10, 2))
})

test_that("count.transitions pools pairs through the excluded group", {
  d <- data.frame(id=as.character(1:6), recruiter.id=c("s","1","2","3","4","5"),
                  net.size=c(5L,4L,6L,3L,5L,4L),
                  grp=factor(c("a","a","b","a","a","a"), levels=c("a","b")))
  rds <- as.rds.data.frame(d, id="id", recruiter.id="recruiter.id",
                           network.size="net.size", population.size=1000)
  rds$grp[3] <- NA
  rds$grp2 <- RDS:::.mint.excluded(rds$grp, RDS:::.excluded.mask(rds$grp))
  t2 <- count.transitions(rds, "grp2")
  ex <- RDS:::.rds.excluded
  # the pairs through the NA subject survive via the internal group
  # (raw table() would delete both and fabricate a 0.1 row for group b)
  expect_equal(unname(t2["a", ex]), 1)
  expect_equal(unname(t2[ex, "a"]), 1)
  expect_equal(unname(t2["a", "a"]), 3)
})

test_that("weights and estimates paths agree, with NAs and subset (T1)", {
  data(fauxmadrona, package="RDS")
  fm <- fauxmadrona
  set.seed(7)
  fm$grp <- factor(sample(c("a","b","c"), nrow(fm), TRUE, prob=c(.5,.3,.2)))
  fm$grp[sample(nrow(fm), 20)] <- NA
  sub <- fm$disease == 0
  keep <- !is.na(fm$grp)[sub]
  for(wt in c("RDS-II","RDS-I","RDS-I (DS)","Gile's SS")){
    N <- if(wt == "Gile's SS") 1000 else NULL
    set.seed(42)
    ee <- suppressWarnings(switch(wt,
      "RDS-II"     = RDS.II.estimates(fm, "grp", subset = disease==0),
      "RDS-I"      = RDS.I.estimates(fm, "grp", subset = disease==0),
      "RDS-I (DS)" = RDS.I.estimates(fm, "grp", subset = disease==0, smoothed=TRUE),
      "Gile's SS"  = RDS.SS.estimates(fm, "grp", N=1000, subset = disease==0)))
    set.seed(42)
    ww <- suppressWarnings(compute.weights(fm, weight.type=wt, outcome.variable="grp",
                                           N=N, subset = disease==0))
    expect_equal(as.numeric(ee@estimate),
                 as.numeric(ht.prop(ww, fm$grp[sub], keep)),
                 tolerance = 1e-8, info = wt)
  }
  # subjects excluded by NA get weight 0 without subsetting
  w.na <- suppressWarnings(compute.weights(fm, weight.type="RDS-I", outcome.variable="grp"))
  expect_true(all(w.na[is.na(fm$grp)] == 0))
  expect_true(all(w.na[!is.na(fm$grp)] > 0))
})

test_that("HCG paths agree and the fit is de-truncated (T1/T8)", {
  data(fauxtime, package="RDS")
  ft <- fauxtime
  ft$var1f <- factor(ifelse(ft$var1 == 1, "y", "n"))
  set.seed(5); ft$var1f[sample(nrow(ft), 40)] <- NA
  keep <- !is.na(ft$var1f)
  set.seed(42)
  ee <- suppressWarnings(RDS.HCG.estimates(ft, "var1f", N=3000))
  set.seed(42)
  ww <- suppressWarnings(hcg.weights(ft, "var1f", N=3000))
  expect_true(all(ww[!keep] == 0))
  expect_equal(as.numeric(ee@estimate),
               as.numeric(ht.prop(ww, ft$var1f, keep)), tolerance = 1e-8)
  # linkage preservation: the internal-group fit differs from refitting on
  # data with the NA rows dropped (which orphans their recruits)
  set.seed(42)
  et <- suppressWarnings(RDS.HCG.estimates(ft[keep, , warn=FALSE], "var1f", N=3000))
  expect_false(isTRUE(all.equal(as.numeric(ee@estimate),
                                as.numeric(et@estimate), tolerance = 1e-6)))
})

test_that("continuous outcomes agree between the two paths (T1 cts)", {
  data(fauxmadrona, package="RDS")
  fm <- fauxmadrona
  set.seed(3); fm$disease[sample(nrow(fm), 25)] <- NA
  keep <- !is.na(fm$disease)
  ee <- suppressWarnings(RDS.I.estimates(fm, "disease"))
  ww <- suppressWarnings(rds.I.weights(fm, "disease"))
  expect_equal(as.numeric(ee@estimate), ht.mean(ww, fm$disease, keep), tolerance = 1e-10)
})

test_that("HCG requires recruitment time on the weights path; bootstraps stay permissive", {
  data(fauxmadrona, package="RDS")  # no recruitment time
  expect_error(RDS.HCG.estimates(fauxmadrona, "disease", N=1000),
               "recruitment time")
  expect_error(suppressWarnings(hcg.weights(fauxmadrona, "disease", N=1000)),
               "recruitment time")
  # small.fraction path does not require time
  expect_silent(suppressWarnings(
    hcg.weights(fauxmadrona, "disease", N=1000, small.fraction=TRUE)))
})

test_that("reserved level name in user data is treated as excluded", {
  data(faux, package="RDS")
  f2 <- faux
  f2$X <- as.character(f2$X)
  f2$X[c(3,10)] <- RDS:::.rds.excluded
  e <- suppressWarnings(RDS.I.estimates(f2, "X"))
  expect_false(RDS:::.rds.excluded %in% names(e@estimate))
  expect_equal(sum(e@estimate), 1, tolerance = 1e-10)
})

test_that("engines fail fast on missing outcomes (T9)", {
  data(fauxmadrona, package="RDS")
  fmf <- fauxmadrona
  fmf$disease <- factor(fmf$disease)
  fmf$disease[5] <- NA
  expect_error(RDS:::salganik.bootstrap.estimates(rds.data=fmf, group.variable="disease",
      number.of.bootstrap.samples=2, estimator.name="RDS-II", N=1000),
      "estimation engine")
  expect_error(RDS:::treeboot.bootstrap.estimates(rds.data=fmf, group.variable="disease",
      number.of.bootstrap.samples=2, estimator.name="RDS-II", N=1000),
      "estimation engine")
})

test_that("impute.degree imputes for subjects with missing trait (P13)", {
  data(faux, package="RDS")
  f <- faux
  f$X[5] <- NA
  f$network.size[5] <- NA
  f$network.size[10] <- NA
  imp <- suppressWarnings(impute.degree(f, trait.variable="X", N=1000))
  expect_true(is.finite(imp[10]))  # observed trait, missing degree
  expect_true(is.finite(imp[5]))   # missing trait AND degree: was NA pre-fix
})
