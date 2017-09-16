context("incidence.R")

test_that("inc",{
  hiv <- rep(FALSE, 5000)
  hiv[1:1000] <- TRUE
  recent <- rep(NA, 5000)
  recent[1:1000] <- FALSE
  recent[1:50] <- TRUE
  expect_true(round(RDS:::incidence(recent,hiv),5) == 0.01895)
  
  

})
  