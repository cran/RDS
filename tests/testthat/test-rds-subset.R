context("rds.data.frame subsetting dispatches to the registered [ method (3.1)")

# Finding 3.1. `[.rds.data.frame` / `[<-.rds.data.frame` were defined with
# @method but no @export, so no S3method line was emitted. A user-level
# `x[i, ]` then dispatched to base `[.data.frame`, which preserves the class
# (masking the problem) but skips the custom method's tree pruning: a subset
# that orphans a recruiter was left with a dangling recruiter.id and no warning,
# silently producing a structurally-invalid rds.data.frame. Fix: add @export so
# the methods register and user-level subsetting runs the RDS methods.

test_that("[ and [<- are registered S3 methods", {
  expect_false(is.null(getS3method("[",   "rds.data.frame", optional = TRUE)))
  expect_false(is.null(getS3method("[<-", "rds.data.frame", optional = TRUE)))
})

test_that("subsetting that orphans a recruiter reseeds it (registered dispatch)", {
  data(faux)
  rid <- attr(faux, "recruiter.id")
  sid <- get.seed.rid(faux)

  # rows 1 & 5: row 5's recruiter (id 134) is not in the subset -> orphaned.
  # Pre-fix (base [.data.frame) left recruiter.id = "134" with no warning;
  # the registered method reseeds it to the seed placeholder and warns.
  expect_warning(faux[c(1, 5), ], "Recruiter removed but child remains")

  g <- suppressWarnings(faux[c(1, 5), ])
  expect_s3_class(g, "rds.data.frame")
  expect_true(all(g[[rid]] == sid))          # no dangling recruiter references

  # warn=FALSE still reseeds, just silently (documented behaviour)
  expect_warning(faux[c(1, 5), warn = FALSE], NA)
  expect_true(all(suppressWarnings(faux[c(1, 5), warn = FALSE])[[rid]] == sid))
})
