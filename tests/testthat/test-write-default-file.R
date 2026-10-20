context("write.netdraw / write.rdsat default file (3.5)")

# Finding 3.5. The default file=NULL branch used substitute(rds.data) on a
# *local* variable, which substitutes the data frame's value, so
# paste(getwd(), <data.frame>) produced one garbage path per column and file()
# errored ("invalid 'description' argument"). Fix: capture the object name from
# the promise at entry via deparse(substitute(x)); the default now writes files
# named after the passed object in the working directory.

test_that("write.netdraw with default file writes object-named files", {
  data(faux)
  wd <- tempfile("wd"); dir.create(wd)
  old <- setwd(wd); on.exit(setwd(old), add = TRUE)
  expect_error(write.netdraw(faux), NA)              # previously: invalid 'description'
  expect_true(file.exists(file.path(wd, "faux.DL")))
})

test_that("write.rdsat with default file writes an object-named file", {
  data(faux)
  wd <- tempfile("wd"); dir.create(wd)
  old <- setwd(wd); on.exit(setwd(old), add = TRUE)
  expect_error(write.rdsat(faux), NA)
  expect_true(file.exists(file.path(wd, "faux.rdsat")))
})
