context("Salganik bootstrap: group-order invariance")

# Regression test for the group-index ordering bug (5.1). The bootstrap built
# `group.names` in first-appearance order but indexed the transition matrix
# `tij` (from count.transitions, which is alphabetical by label) and the
# per-group degree/outcome pools inconsistently, so the recruitment chain used
# the wrong transition row and drew degrees from the wrong group whenever
# appearance order != alphabetical order.
#
# The invariant: the bootstrap distribution of "proportion in a given
# partition group" cannot depend on what we *call* the groups. Relabeling the
# same partition with labels that sort differently must leave the per-group
# bootstrap SEs unchanged (up to the ~1e-8 floating-point noise that a
# different factor-label string set introduces). Relabeling touches only the
# group labels -- not the tree, wave order, ids, degrees or partition -- so
# the wave-sorted internal processing is unaffected.

salganik_se <- function(rds.data, group.variable, seed = 42, B = 120) {
  set.seed(seed)
  iv <- suppressWarnings(RDS.bootstrap.intervals(
    rds.data, group.variable, weight.type = "RDS-II",
    uncertainty = "Salganik", N = 1000,
    number.of.bootstrap.samples = B, verbose = FALSE))
  bs <- attr(iv, "bsresult")            # replicate x group matrix
  se <- apply(bs, 2, stats::sd)
  se[order(names(se))]
}

test_that("Salganik bootstrap SEs are invariant to group relabeling", {
  data(faux)

  # A 3-group partition on the faux tree, with strongly group-separated
  # degrees so that drawing from the wrong group's degree pool (the bug)
  # produces a large, reliable effect.
  set.seed(101)
  part <- sample(1:3, nrow(faux), replace = TRUE)
  deg  <- c(35L, 18L, 6L)[part] + rpois(nrow(faux), 2)

  build <- function(labels) {
    d <- faux
    d[[attr(d, "network.size.variable")]] <- deg
    # alphabetical levels so as.factor/count.transitions and our labels agree
    d$grp <- factor(labels[part], levels = sort(unique(labels)))
    d
  }

  # Labeling A: a,b,c  -> sorted a<b<c, partition k -> group k (identity).
  # Labeling B: m,p,c  -> sorted c<m<p, so partition 1->m(idx2), 2->p(idx3),
  #                       3->c(idx1): appearance order != alphabetical order,
  #                       which is exactly what triggered the bug.
  A <- salganik_se(build(c("a", "b", "c")), "grp")
  B <- salganik_se(build(c("m", "p", "c")), "grp")

  # Map B's labels back onto A's labels via the shared partition.
  remap <- c(m = "a", p = "b", c = "c")
  B_on_A <- setNames(as.numeric(B), remap[names(B)])[c("a", "b", "c")]
  A_ord  <- A[c("a", "b", "c")]

  expect_equal(unname(B_on_A), unname(A_ord), tolerance = 1e-4)
})

test_that("Treeboot bootstrap SEs are invariant to group relabeling", {
  # Regression for the treeboot half of 5.2: value columns were labeled in
  # first-appearance order (`group.names`) while filled positionally from
  # `bs.results` rows, which are in RDS.estimates.local's sorted-label order.
  # For a non-alphabetical appearance order the per-group SEs were attached to
  # the wrong group labels (a cyclic scramble). Treeboot resamples the tree
  # structure, so relabeling must leave the per-group SEs unchanged.
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
  tb_se <- function(d) {
    set.seed(42)
    iv <- suppressWarnings(RDS.bootstrap.intervals(
      d, "grp", weight.type = "RDS-II", uncertainty = "Treeboot",
      N = 1000, number.of.bootstrap.samples = 120, verbose = FALSE))
    se <- apply(attr(iv, "bsresult"), 2, stats::sd)
    se[order(names(se))]
  }

  A <- tb_se(build(c("a", "b", "c")))   # sorted a<b<c = partition order
  B <- tb_se(build(c("m", "p", "c")))   # sorted c<m<p : appearance != sorted
  remap  <- c(m = "a", p = "b", c = "c")
  B_on_A <- setNames(as.numeric(B), remap[names(B)])[c("a", "b", "c")]

  expect_equal(unname(B_on_A), unname(A[c("a", "b", "c")]), tolerance = 1e-4)
})
