context("[.rds.data.frame + recruitment-time repair (3.7 parts 1 & 3)")

# 3.7 part 1: `rid <- id <- attr(x,"recruiter.id")` clobbered `id` with the
# recruiter column, so the seed-repair guard checked recruiter.id twice and
# subsetting to columns without the id column errored obscurely instead of
# returning a plain data.frame.
test_that("column subset without the id column returns a plain data.frame", {
  data(faux)
  ridcol <- attr(faux, "recruiter.id")
  nscol  <- attr(faux, "network.size.variable")
  z <- faux[, c(ridcol, nscol)]
  expect_false(is.rds.data.frame(z))
  expect_equal(ncol(z), 2L)
})

# 3.7 part 3: make.valid.recruitment.time terminated on stagnant invalid COUNT
# rather than the stagnant invalid SET, so a chain that repairs one subject per
# pass (constant count, shifting set) was clobbered with min(time) and errored
# ("Some subjects recruited before their recruiter"). The linear chain
# seed(t=5) -> a(1) -> b(2) -> c(3) exercises exactly that.
test_that("cascading out-of-order recruitment times are repaired, not clobbered", {
  df <- data.frame(id = c("s","a","b","c"), recruiter.id = c("0","s","a","b"),
                   network.size = c(4,3,5,2), t = c(5,1,2,3),
                   stringsAsFactors = FALSE)
  d <- suppressWarnings(suppressMessages(as.rds.data.frame(df, id = "id",
         recruiter.id = "recruiter.id", network.size = "network.size", time = "t")))
  expect_equal(unname(get.recruitment.time(d)), c(5, 5, 5, 5))
})

# 3.7 part 2: single-index forms were misrouted (nargs() lost in the wrapper) --
# `x["X"]` dispatched as row selection and `x["newcol"] <- v` recycled across a
# row, both erroring. Fix: branch on the index count (nargs adjusted for the
# extra warn/drop args). Single index => column selection.
test_that("single-index selection and assignment work (part 2)", {
  data(faux)
  # single index -> column selection (plain data.frame, base semantics)
  z <- faux["X"]
  expect_false(is.rds.data.frame(z))
  expect_equal(ncol(z), 1L)
  expect_equal(nrow(z), nrow(faux))
  expect_equal(ncol(faux[c("X", "Y")]), 2L)

  # single-index column assignment adds a column, keeps rds.data.frame
  d <- faux
  d["newz"] <- seq_len(nrow(d))
  expect_s3_class(d, "rds.data.frame")
  expect_true("newz" %in% names(d))

  # two-index forms unaffected: row subset reconstructs, [i,j] drops to a vector
  z2 <- suppressWarnings(faux[c(1, 5), ])
  expect_s3_class(z2, "rds.data.frame")
  expect_equal(nrow(z2), 2L)
  expect_equal(length(faux[1:3, "X"]), 3L)
})
