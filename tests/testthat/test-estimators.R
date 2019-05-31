# TODO: Add comment
# 
# Author: ianfellows
###############################################################################



library(testthat)
library(RDS)

context("rds-freq.R")

df <- structure(list(id = 1:10, 
                     recruiter.id = 0:9, 
                     degree = c(16, 13, 13, 7, 7, 7, 13, 8, 12, 15), x = c(-0.406194954884142, -2.31045703106648, 
                         0.492457818020673, 0.520175024870627, 0.921941578183436, -0.158956546457245, 
                         0.202097220140023, -0.778286503628266, 1.32036321437896, 0.264852950094362
                     ), 
                     y = structure(c(1L, 2L, 1L, 2L, 2L, 2L, 2L, 1L, 2L, 1L), 
                                   .Label = c("FALSE",  "TRUE"), class = "factor"), 
                     z = structure(c(1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L), 
                                   .Label = c("a", "b"), class = "factor")), 
                .Names = c("id", "recruiter.id", "degree", "x", "y", "z"), 
                row.names = c(NA, -10L), 
                class = "data.frame")
df$time <- 1:10
rdsdat <- as.rds.data.frame(df, network.size = "degree", population.size = c(NA, 100, NA), time="time")
rdsdat$m <- c(NA, "t","t","f",NA,NA,NA,NA,NA,NA)

test_that("rds-MA",{
  v <- MA.estimates(rdsdat,trait.variable = "y", N=100,verbose=FALSE,number.of.coupons = 1)
  expect_true(v$estimate[1] > .26 && v$estimate[1] < .38)
  v <- MA.estimates(rdsdat,trait.variable = "z", N=100,verbose=FALSE,number.of.coupons = 1)
  expect_true(v$estimate[1] < .05)
  expect_warning(v <- MA.estimates(rdsdat,trait.variable = "m", N=100,verbose=FALSE,number.of.coupons = 1))
  expect_true(v$estimate[1] > .8)
})


test_that("rds-hcg",{
  est <- RDS.HCG.estimates(rds.data=rdsdat,outcome.variable='x')
  expect_true(round(est$estimate,2) < 0.1) 
  est <- RDS.HCG.estimates(rds.data=rdsdat,outcome.variable='y')
  expect_true(round(est$estimate[1],1) <= .4)   
  est <- RDS.HCG.estimates(rds.data=rdsdat,outcome.variable='z')
  
  rdsdat2 <- rdsdat
  rdsdat2$degree[4:6] <- NA
  
  expect_warning(RDS.HCG.estimates(rds.data=rdsdat2,outcome.variable='x'))
  expect_warning(RDS.HCG.estimates(rds.data=rdsdat2,outcome.variable='y'))
  expect_warning(RDS.HCG.estimates(rds.data=rdsdat2,outcome.variable='z'))

})




