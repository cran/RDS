

library(testthat)
library(RDS)

context("homophily.R")

test_that("homophily",{
  data(fauxmadrona)
  res <- homophily.estimates(fauxmadrona,outcome.variable="disease",
                             N=1000,weight.type="RDS-II")
  expect_equal(structure(1.60397806961622, class = "table"),res@estimate)
})

# 4.13: recruitment homophily on an all-seeds sample (no wave-1 recruits) used
# to crash with "NA/NaN argument" because wave.one.start <- match(1, wave) is NA
# and 1:(NA-1) errors. It now stops with a clear message.
test_that("recruitment homophily on an all-seeds sample errors cleanly (4.13)",{
  n <- 20
  df <- data.frame(id = as.character(1:n),
                   recruiter.id = rep("seed", n),
                   net.size = rep(3L, n),
                   grp = factor(rep(c("a","b"), length.out = n)))
  rds <- as.rds.data.frame(df, id = "id", recruiter.id = "recruiter.id",
                           network.size = "net.size", population.size = 1000)
  expect_equal(sort(unique(get.wave(rds))), 0)   # no wave-1 recruits
  expect_error(homophily.estimates(rds, "grp", recruitment = TRUE),
               "No recruitment")
})


