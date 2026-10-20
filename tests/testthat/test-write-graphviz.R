context("write.graphviz: correct edges on unsorted data (3.3)")

# Finding 3.3. write.graphviz reordered the edge list by wave and *then* filtered
# it with the un-reordered `rid != sid` mask. On data whose row order differs
# from wave order the masks misaligned, so seed-placeholder edges leaked through
# (bogus "<seed>" -> x, a node that should emit nothing) and real recruitment
# edges were dropped. Fix: drop the (irrelevant, for a digraph) wave reorder so
# the mask stays aligned, and use drop=FALSE so a single-edge dataset doesn't
# collapse `el` to a vector and crash `el[,1]`.

edges_of <- function(rds) {
  f <- tempfile(fileext = ".dot"); on.exit(unlink(f))
  write.graphviz(rds, f)
  txt <- paste(readLines(f, warn = FALSE), collapse = "\n")
  e <- regmatches(txt, gregexpr('"[^"]+" -> "[^"]+"', txt))[[1]]
  sort(gsub('"', '', gsub(' -> ', ' ', e)))
}

test_that("edges are correct and seed-placeholder-free on row-shuffled data", {
  data(faux)
  set.seed(1)
  sh  <- faux[sample(nrow(faux)), ]                 # row order != wave order
  sid <- as.character(get.seed.rid(sh))
  rid <- as.character(get.rid(sh)); id <- as.character(get.id(sh))
  correct <- sort(paste(rid[rid != sid], id[rid != sid]))   # real recruiter -> recruit edges

  got <- edges_of(sh)
  expect_false(any(startsWith(got, paste0(sid, " "))))      # no "<seed>" -> x edges
  expect_equal(got, correct)                                # exactly the real tree edges
})

test_that("a single-edge dataset does not error (drop=FALSE)", {
  d <- as.rds.data.frame(
    data.frame(id = c("1", "2"), recruiter.id = c("0", "1"),
               network.size = c(3, 4)),
    id = "id", recruiter.id = "recruiter.id", network.size = "network.size")
  expect_error(write.graphviz(d, tempfile(fileext = ".dot")), NA)
})
