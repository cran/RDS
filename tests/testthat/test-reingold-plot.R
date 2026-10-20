context("reingold.tilford.plot edge cases (3.9)")

# 3.9 part 1: el[rid!=sid,] lacked drop=FALSE, so a 2-node tree (exactly one
# non-seed edge) collapsed el to a vector and crashed at order(el[,3]) with
# "incorrect number of dimensions".
test_that("a 2-node tree does not crash the reingold-tilford plot", {
  skip_if_not_installed("igraph")
  d <- as.rds.data.frame(
    data.frame(id = c("1", "2"), recruiter.id = c("0", "1"),
               network.size = c(3, 4)),
    id = "id", recruiter.id = "recruiter.id", network.size = "network.size")
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  expect_error(suppressWarnings(reingold.tilford.plot(d)), NA)
})
