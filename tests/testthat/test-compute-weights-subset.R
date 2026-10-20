context("compute.weights: subset= expression support (4.5)")

# Finding 4.5. compute.weights only eval'd `subset` when it was a bare name, and
# even then eval'd the *forced promise* in the caller rather than the captured
# expression in rds.data. So any column-based subset -- an expression
# (X=="blue", the documented seed>0) or even a bare column name -- errored with
# "object not found". Fix: eval the substituted expression `se` inside rds.data
# with the caller as fallback (the pattern already used in RDS.estimates.local /
# differential.activity): subset <- eval(se, rds.data, parent.frame()).

test_that("subset= accepts a column expression and a bare column name", {
  data(faux)
  full <- compute.weights(faux, "RDS-II")
  blue <- faux$X == "blue"

  # expression form (previously: object 'X' not found)
  we <- compute.weights(faux, "RDS-II", subset = X == "blue")
  expect_equal(length(we), sum(blue))
  expect_equal(unname(we), unname(full[blue]))

  # bare column-name form (previously: object 'keep' not found)
  faux$keep <- blue
  wn <- compute.weights(faux, "RDS-II", subset = keep)
  expect_equal(unname(wn), unname(full[blue]))
})
