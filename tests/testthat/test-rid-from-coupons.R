context("rid.from.coupons: mandatory subject.coupon + tree reconstruction (3.11)")

# 3.11. subject.coupon defaulted to NULL, so data[[subject.coupon]] crashed
# ("attempt to select less than one element") -- it is mandatory in practice and
# is now a required argument.

test_that("subject.coupon is required and the recruiter tree is reconstructed", {
  d <- data.frame(own = c("s1", "c1", "c2", "c3"),
                  cpn.1 = c("c1", "c3", "c4", NA),
                  cpn.2 = c("c2", NA, NA, NA),
                  stringsAsFactors = FALSE)
  expect_error(rid.from.coupons(d, coupon.variables = c("cpn.1", "cpn.2")),
               "subject.coupon")
  rid <- rid.from.coupons(d, subject.coupon = "own",
                          coupon.variables = c("cpn.1", "cpn.2"))
  expect_equal(rid, c("seed", "1", "1", "2"))
})
