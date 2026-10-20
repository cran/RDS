context("as.rds.data.frame preserves the notes attribute (3.6)")

# Finding 3.6. The `notes` default's else branch returned attr(df,"time")
# (copy-pasted from the time default), so every reconstruction -- including
# every subset via the now-registered [.rds.data.frame -- replaced the notes
# attribute with the time-variable name ("MY NOTES" -> "wave").

test_that("notes survives reconstruction and subsetting", {
  data(faux)
  faux$wave <- get.wave(faux)
  d <- as.rds.data.frame(faux, notes = "MY NOTES", time = "wave")
  expect_equal(attr(d, "notes"), "MY NOTES")
  expect_equal(attr(as.rds.data.frame(d), "notes"), "MY NOTES")   # reconstruct
  expect_equal(attr(d[1:5, ], "notes"), "MY NOTES")               # subset
})
