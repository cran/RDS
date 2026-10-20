context("HCG bootstrap replicate generation")

# Regression tests for the HCG.bootstrap fixes:
# - 5.5: as.factor(out.boot); levels(out.boot) <- lev relabeled groups
#   positionally among *observed* codes, so a replicate that missed a rare
#   group shifted every higher group's label down by one.
# - 5.6: the degree-draw fallback compared a factor to an integer index,
#   used the sample(d, 1, prob=d) scalar gotcha, and crashed on groups with
#   no observed degrees (e.g. empty contingency-table cells under the
#   independence null).

three.level.fauxtime <- function(p.rare = 0.006, seed = 7) {
  data(fauxtime, envir = environment())
  set.seed(seed)
  r <- runif(nrow(fauxtime))
  fauxtime$Z <- factor(ifelse(r < p.rare, "b",
                              ifelse(fauxtime$var1 == 1, "a", "c")),
                       levels = c("a", "b", "c"))
  fauxtime
}

test_that("replicates that miss a group do not relabel the remaining groups", {
  fx <- three.level.fauxtime()   # group b has a single subject
  set.seed(99)
  tabs <- suppressWarnings(RDS:::HCG.bootstrap(
    fx, "Z", number.of.bootstrap.samples = 200, N = 3000,
    small.fraction = TRUE, verbose = FALSE,
    fun = function(x) table(factor(x[["Z"]], levels = c("a", "b", "c")))))
  m <- do.call(rbind, tabs)

  # The rare group must actually be missed in some replicates, otherwise
  # this test exercises nothing.
  expect_gt(sum(m[, "b"] == 0), 0)
  # Pre-fix, exactly those replicates reported the *largest* group as absent
  # (all its members relabeled). The large groups must never vanish.
  expect_equal(sum(m[, "a"] == 0), 0)
  expect_equal(sum(m[, "c"] == 0), 0)
})

test_that("degree draws from a single-member group return its degree", {
  fx <- three.level.fauxtime()
  b.degree <- fx[[attr(fx, "network.size.variable")]][fx$Z == "b"]
  expect_length(b.degree, 1)

  set.seed(99)
  degs <- suppressWarnings(RDS:::HCG.bootstrap(
    fx, "Z", number.of.bootstrap.samples = 100, N = 3000,
    small.fraction = TRUE, verbose = FALSE,
    fun = function(x) {
      i <- x[["Z"]] == "b" & !is.na(x[["Z"]])
      unique(x[[attr(x, "network.size.variable")]][i])
    }))
  drawn <- unique(unlist(degs))
  # The sample(d, 1, prob=d) scalar gotcha drew uniformly from 1:d instead.
  expect_equal(drawn, b.degree)
})

test_that("a factor level with zero observations does not crash the bootstrap", {
  data(fauxtime)
  fe <- fauxtime
  fe$Z <- factor(ifelse(fe$var1 == 1, "a", "c"), levels = c("a", "b", "c"))

  set.seed(1)
  tabs <- suppressWarnings(RDS:::HCG.bootstrap(
    fe, "Z", number.of.bootstrap.samples = 25, N = 3000,
    small.fraction = TRUE, verbose = FALSE,
    fun = function(x) table(factor(x[["Z"]], levels = c("a", "b", "c")))))
  m <- do.call(rbind, tabs)
  expect_equal(nrow(m), 25)
  # observed groups keep sane counts
  expect_true(all(m[, "a"] + m[, "b"] + m[, "c"] == sum(!is.na(fe$Z))))
})

test_that("bootstrap.contingency.test handles structurally empty cells", {
  data(fauxtime)
  f4 <- fauxtime
  f4$v2 <- ifelse(seq_len(nrow(f4)) %% 2 == 0, "x", "y")
  f4$v2[f4$var1 == 1] <- "x"   # cell (var1=1, v2=y) is empty by construction

  set.seed(1)
  ct <- suppressWarnings(bootstrap.contingency.test(
    f4, "var1", "v2", number.of.bootstrap.samples = 50, verbose = FALSE))
  expect_true(is.finite(ct$p.value))
  # the empty cell here is strong evidence of dependence
  expect_lt(ct$p.value, 0.1)
})
