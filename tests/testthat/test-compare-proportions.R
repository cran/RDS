context("RDS.compare.proportions: reversed match + swapped indices (4.6)")

# Minimal mock rds.interval.estimate: RDS.compare.proportions reads
# attr(x,"bsresult") (a bootstrap matrix, colnames = level names) and
# x$estimate (named observed proportions). RDS.compare.two.proportions also
# reads x$interval and the list attr(data,"variables").
mk_int <- function(est, bs) {
  o <- list(estimate = est, interval = rep(0, 6 * length(est)))
  attr(o, "bsresult") <- bs
  class(o) <- "rds.interval.estimate"
  o
}

# Bug 1 (lines 30/39/83/93): obs <- tmp[match(names(tmp), levs)] applies the
# INVERSE permutation -- it should be tmp[match(levs, names(tmp))]. So obs is
# mis-assigned to the bootstrap columns whenever the estimate's name order
# differs from the bootstrap column order (>= 3 levels). Test: permuting only
# the bootstrap columns (estimate fixed) must not change the name-matched
# p-values.
test_that("p-values are invariant to bootstrap column order (bug 1)", {
  set.seed(1); n <- 600
  d_a <- rnorm(n, .2, .03); d_b <- rnorm(n, .3, .03); d_c <- rnorm(n, .5, .03)
  est1 <- c(a = .2, b = .3, c = .5)                       # fixed name order a,b,c
  bs_abc <- cbind(a = d_a, b = d_b, c = d_c)              # cols a,b,c
  bs_cab <- cbind(c = d_c, a = d_a, b = d_b)              # same data, cols c,a,b
  i2 <- mk_int(c(x = .4, y = .6),
               cbind(x = rnorm(n, .4, .03), y = rnorm(n, .6, .03)))

  set.seed(9); r1 <- as.matrix(RDS.compare.proportions(mk_int(est1, bs_abc), i2, M = 4000))
  set.seed(9); r2 <- as.matrix(RDS.compare.proportions(mk_int(est1, bs_cab), i2, M = 4000))
  lv <- c("a", "b", "c")
  expect_equal(r1[lv, ], r2[lv, ], tolerance = 1e-6)
})

# Bug 2 (line 101, RDS.compare.two.proportions only): bs <- bs2[,i] - bs1[,j]
# with i over nlev1 and j over nlev2 -- but bs2 has nlev2 columns and bs1 has
# nlev1. When the level counts differ it is a subscript error. Should be
# bs1[,i] - bs2[,j] (as in the first function, line 47).
test_that("differing level counts do not crash (bug 2)", {
  set.seed(1); n <- 400
  i1 <- mk_int(c(a = .2, b = .3, c = .5),
               matrix(rnorm(n * 3, rep(c(.2, .3, .5), each = n), .03),
                      ncol = 3, dimnames = list(NULL, c("a", "b", "c"))))   # 3 levels
  i2 <- mk_int(c(x = .4, y = .6),
               matrix(rnorm(n * 2, rep(c(.4, .6), each = n), .03),
                      ncol = 2, dimnames = list(NULL, c("x", "y"))))        # 2 levels
  dat <- data.frame(v1 = 1, v2 = 1)
  attr(dat, "variables") <- list(v1 = i1, v2 = i2)

  grDevices::pdf(NULL); oldpar <- graphics::par(no.readonly = TRUE)
  on.exit({ graphics::par(oldpar); grDevices::dev.off() })
  expect_error(
    RDS.compare.two.proportions(dat, c("v1", "v2"),
                                number.of.bootstrap.samples = 1000, plot = FALSE),
    NA)
})
