context("bootstrap.incidence: outcomes align with sorted rows (5.4)")

# Finding 5.4. bootstrap.incidence sorts the data (rds <- rds.data[order(time,
# wave),]) but extracted hiv/recent from the UNSORTED rds.data and then attached
# the derived combination outcome to the sorted `rds` -- so outcomes were
# scrambled onto the wrong subjects for any input not already time-sorted
# (fauxtime: 502/511 rows out of order; estimate 0.146 vs correct 0.100). The
# point estimate ($Incidence) comes from inc.func(rds, FALSE, FALSE) with no
# bootstrap randomness, so it must be exactly invariant to input row order.

test_that("incidence point estimate is invariant to input row order", {
  data(fauxtime)
  d <- fauxtime
  d$hiv    <- (seq_len(nrow(d)) %% 3 == 0)
  d$recent <- (seq_len(nrow(d)) %% 7 == 0) & d$hiv
  tm <- get.recruitment.time(d)
  w  <- get.wave(d)

  inc <- function(x) suppressWarnings(bootstrap.incidence(
    x, "recent", "hiv", weight.type = "RDS-II",
    number.of.bootstrap.samples = 2, verbose = FALSE)$Incidence)

  # fauxtime is not time-sorted; pre-fix the unsorted and pre-sorted runs
  # disagreed (0.146 vs 0.100). They must now match.
  expect_equal(inc(d), inc(d[order(tm, w), ]), tolerance = 1e-8)
})
