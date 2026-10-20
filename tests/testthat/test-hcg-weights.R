context("hcg.weights degenerate inputs")

# Regression tests for the single-outcome-level early exit in hcg.weights.
# It previously returned N * rep(1/n, n), which with N=NULL silently produced
# numeric(0) (NULL arithmetic). It now falls back to RDS-II (vh.weights) when
# N is unknown and Gile's SS weights when N is supplied; both handle zero and
# missing degrees.

single.level.faux <- function(degree.fix = NULL) {
  data(faux, envir = environment())
  df <- as.data.frame(faux)
  df$X <- "blue"
  # strip the population size so N=NULL survives end to end
  for (a in c("population.size.low", "population.size.mid", "population.size.high"))
    attr(df, a) <- NULL
  if (!is.null(degree.fix))
    df$network.size[degree.fix$i] <- degree.fix$v
  suppressWarnings(
    as.rds.data.frame(df, id = "id", recruiter.id = "recruiter.id",
                      network.size = "network.size")
  )
}

test_that("single-level outcome with N=NULL gives RDS-II style weights", {
  f0 <- single.level.faux()
  w <- suppressWarnings(hcg.weights(f0, "X", N = NULL, small.fraction = TRUE))

  n <- nrow(f0)
  expect_length(w, n)                       # was numeric(0) pre-fix
  expect_true(all(is.finite(w)))
  expect_true(all(w > 0))
  # RDS-II fallback: weights proportional to 1/degree
  expect_equal(w / sum(w),
               (1 / f0$network.size) / sum(1 / f0$network.size))
})

test_that("single-level outcome with N supplied gives Gile's SS weights", {
  f0 <- single.level.faux()
  set.seed(1)
  w <- suppressWarnings(hcg.weights(f0, "X", N = 3000, small.fraction = TRUE))

  expect_length(w, nrow(f0))
  expect_true(all(is.finite(w)))
  expect_true(all(w > 0))
  expect_gt(sum(w), 0)
})

test_that("single-level fallback tolerates zero and missing degrees", {
  # One NA degree: previously the N=NULL branch of an intermediate fix
  # poisoned every weight (all NA); vh.weights assigns that subject weight 0.
  f.na <- single.level.faux(list(i = 5, v = NA))
  w <- suppressWarnings(hcg.weights(f.na, "X", N = NULL, small.fraction = TRUE))
  expect_length(w, nrow(f.na))
  expect_true(all(is.finite(w)))
  set.seed(1)
  wg <- suppressWarnings(hcg.weights(f.na, "X", N = 3000, small.fraction = TRUE))
  expect_true(all(is.finite(wg)))

  # One zero degree: raw 1/degree would give Inf/NaN; both branches must
  # stay finite.
  f.z <- single.level.faux(list(i = 7, v = 0))
  w <- suppressWarnings(hcg.weights(f.z, "X", N = NULL, small.fraction = TRUE))
  expect_true(all(is.finite(w)))
  set.seed(1)
  wg <- suppressWarnings(hcg.weights(f.z, "X", N = 3000, small.fraction = TRUE))
  expect_true(all(is.finite(wg)))
})

test_that("RDS.HCG.estimates works end to end on a single-level outcome", {
  f0 <- single.level.faux()
  est <- suppressWarnings(RDS.HCG.estimates(f0, "X", small.fraction = TRUE))
  expect_equal(as.numeric(est@estimate), 1)
  expect_length(est@weights, nrow(f0))
})
