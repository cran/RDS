
library(testthat)
library(RDS)

context("rds.data.frame.R")

test_that("rds.data.frame",{
  dat <- data.frame(id = c(1,2,3,4,5),
  recruiter.id = c(0,1,1,2,0),
  network.size.variable = c(2,4,3,5,4))
  rds <- as.rds.data.frame(dat)
  dat$time <- c(1,2,3,4,5)
  rds <- as.rds.data.frame(dat,time="time")
  expect_true(length(get.recruitment.time(rds)) == 5)
  expect_equal(class(get.recruitment.time(rds)), "numeric")
  
  dat$time <- c("1970-01-02", "1970-01-03", "1970-01-04", "1970-01-05", "1970-01-06")
  rds <- as.rds.data.frame(dat,time="time")
  expect_true(length(get.recruitment.time(rds)) == 5)
  expect_equal(class(get.recruitment.time(rds, to.numeric = FALSE)), c("POSIXct", "POSIXt")) 
  expect_equal(class(get.recruitment.time(rds)), "numeric") 
  assert.valid.rds.data.frame(rds)
})