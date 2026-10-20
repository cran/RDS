context("Gile's SS C samplers")

# Regression test for the clamp-ordering bug in getinclCstacked.c (and its
# twin in bsC.c): in the offby > 0 rebalancing loop, cdif[0] was copied from
# dif[0] BEFORE dif[0] was clamped at zero, so a negative first-class deficit
# (smallest degree class over-represented) shifted the whole cumulative
# probability vector down and starved early classes of selection mass.

test_that("offby>0 rebalancing allocates by deficit, not by position", {
  # Crafted case: class 1 is over-represented (props2[1]=0, nbyclass 10/20),
  # classes 2 and 3 are equally deficient. The one added member must go to
  # class 2 or class 3 with equal probability. Class identity is recovered
  # from the returned (sorted) size vector: adding to class 2 produces the
  # normalized product 12/42; adding to class 3 produces 24/44.
  set.seed(1)
  picks <- replicate(2000, {
    ret <- .C("getinclCstacked",
              nbyclass = as.integer(c(10, 5, 5)),
              size = as.double(c(1, 2, 4)),
              props2 = as.double(c(0, 0.5, 0.5)),
              offby = as.integer(1),
              N = as.integer(20),
              K = as.integer(3),
              n = as.integer(1),
              samplesize = as.integer(1),
              Nk = integer(3), PACKAGE = "RDS")
    s <- sort(ret$size, decreasing = TRUE)
    if (any(abs(s - 12/42) < 1e-9)) 2L else 3L
  })
  n2 <- sum(picks == 2L)

  # Binomial(2000, 0.5): sd ~ 22. Bounds are ~6.7 sd wide; the pre-fix code
  # deterministically gave n2 == 0.
  expect_gt(n2, 850)
  expect_lt(n2, 1150)
})

test_that("gile.ss.weights returns sane weights end to end", {
  data(faux)
  deg <- as.numeric(faux[[attr(faux, "network.size.variable")]])
  set.seed(1)
  w <- gile.ss.weights(deg, N = 1000)
  expect_length(w, length(deg))
  expect_true(all(is.finite(w)))
  expect_true(all(w > 0))
  # higher degree => lower weight overall (Monte Carlo noise in the simulated
  # inclusion probabilities permits tiny local wiggles, so test the trend,
  # not strict monotonicity)
  expect_lt(stats::cor(deg, w, method = "spearman"), -0.9)
})

# Finding 7.3. In the offby<0 (removal) loop, the temp==0.0 fallback protected
# only count-==1 classes, so a count-0 class received a positive selection
# weight and could be decremented to -1, -2, ... Negative class counts flowed
# back to R (1/probs weights). The fix ports bsC.c's treatment: clamp <=1 in the
# fallback plus a last-resort branch that removes only from non-empty classes.
# Reachable via gile.ss.weights(N < n), which is now also rejected up front.
test_that("gile.ss.weights rejects a population smaller than the sample (7.3)", {
  expect_error(gile.ss.weights(1:8, N = 5),
               "cannot be smaller than the sample")
  # census / normal cases are unaffected
  set.seed(1); expect_true(all(gile.ss.weights(1:8, N = 8)  >= 0))
  set.seed(1); expect_true(all(gile.ss.weights(1:8, N = 20) >  0))
})

test_that("getinclCstacked never drives a class count negative (7.3)", {
  # Directly stress the C allocator under N < #classes (the R guard's blind
  # spot), inspecting the in-place nbyclass. Pre-fix this reached -4.
  worst <- 0L
  for (K in 5:8) for (Nt in 2:(K - 1)) for (s in 1:20) {
    set.seed(s + 100 * K + 10 * Nt)
    S  <- Nt + sample(1:5, 1)
    nb <- as.integer(tabulate(sample(1:K, S, replace = TRUE), nbins = K))
    g  <- rgamma(K, 0.3); p2 <- g / sum(g)
    ret <- .C("getinclCstacked",
      nbyclass = nb, size = as.double(1:K), props2 = as.double(p2),
      offby = as.integer(-(S - Nt)), N = as.integer(Nt), K = as.integer(K),
      n = as.integer(Nt), samplesize = as.integer(1), Nk = integer(K),
      PACKAGE = "RDS")
    worst <- min(worst, min(ret$nbyclass))
  }
  expect_gte(worst, 0L)
})
