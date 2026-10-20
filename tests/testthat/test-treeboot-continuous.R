context("Treeboot uncertainty rejects continuous outcomes (5.12)")

# Finding 5.12. treeboot.bootstrap.estimates assumes a categorical outcome (one
# column per unique value), so a numeric/continuous outcome produced a garbled
# matrix and crashed downstream ("length of 'dimnames' ..."). Rather than
# support it, Treeboot now errors clearly for continuous outcomes; categorical
# Treeboot is unaffected.

test_that("continuous outcome under Treeboot errors clearly", {
  data(fauxmadrona)
  d <- fauxmadrona
  set.seed(1)
  d$cts <- as.numeric(d$disease) + rnorm(nrow(d), 0, 0.01)
  expect_error(
    RDS.bootstrap.intervals(d, "cts", weight.type = "RDS-II",
      uncertainty = "Treeboot", number.of.bootstrap.samples = 10, verbose = FALSE),
    "does not support continuous")
})

test_that("categorical outcome under Treeboot still works", {
  data(faux)
  iv <- RDS.bootstrap.intervals(faux, "X", weight.type = "RDS-II",
    uncertainty = "Treeboot", number.of.bootstrap.samples = 20, verbose = FALSE)
  expect_true(all(is.finite(iv$interval)))
})

# Finding 5.13. treeboot generated seed ids "0","1",...,"11" in the same flat
# namespace as child ids paste0(parent, i), so the 12th seed ("11") collided
# with the first child of the 2nd seed ("1"+"1"). With check.valid=FALSE the
# duplicate ids passed and match(recruiter.id, id) miswired the tree. Ids are
# now namespaced: seeds "s0","s1",... and children paste0(parent, ".", i).
test_that("treeboot bootstrap ids are unique with many small seed trees (5.13)", {
  nseed <- 20
  ids  <- c(paste0("s", 1:nseed), paste0("c", 1:nseed))
  rids <- c(rep("seed", nseed), paste0("s", 1:nseed))   # each c_k recruited by s_k
  df <- data.frame(id = ids, recruiter.id = rids, net.size = rep(3L, 2 * nseed),
                   X = factor(rep(c("a", "b"), length.out = 2 * nseed)),
                   stringsAsFactors = FALSE)
  rds <- as.rds.data.frame(df, id = "id", recruiter.id = "recruiter.id",
                           network.size = "net.size", population.size = 1000)
  id <- get.id(rds); rid <- get.rid(rds)
  rrow <- match(rid, id); rrow[is.na(rrow)] <- 0
  n <- length(id); seed.rows <- which(get.wave(rds) == 0)
  recruits <- lapply(seq_len(n), function(i) which(rrow == i))
  outcomes <- data.frame(rds)["X"]; netsize <- rds[["net.size"]]

  treeboot <- getFromNamespace("treeboot", "RDS")
  # pre-fix this produced duplicate id "11" on the very first seed
  for (s in 1:30) {
    set.seed(s)
    tb <- treeboot(seed.rows, recruits, netsize, outcomes,
                   get.population.size(rds), fixed.size = TRUE)
    bids <- get.id(tb)
    expect_false(any(duplicated(bids)),
                 info = paste("duplicate bootstrap ids at seed", s))
  }
})
