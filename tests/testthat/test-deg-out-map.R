context("deg.out.map keeps every subject in its pool (7.2)")

# Finding 7.2. The continuous-outcome imputation pool deg.out.map[[i]][[d]] was
# built with `length(lis) <= d[j]` guarding the "initialize vs append" branch.
# Because lis[[d]] grows the list to length d, length(lis) equals the running-max
# degree, so a subject whose degree equals that max took the *initialize* branch
# and OVERWROTE the existing pool instead of appending -- silently dropping
# every earlier outcome at that degree (constant-degree worst case: 60 outcomes
# collapse to {20,40,60}). Fixed to `<` in gile-bootstrap.R and hcg-bootstrap.R.
#
# Direct guard: capture the real deg.out.map via trace and assert the pooled
# outcomes account for every subject (sum of pool sizes == n). Pre-fix this is
# short by the number of overwrites (one per group max on fauxmadrona: 497/500).

# Sum of all pool sizes in the deg.out.map built inside `fn` during `run()`.
dom_total <- function(fn, run) {
  if (exists(".dom_capture", envir = globalenv()))
    rm(".dom_capture", envir = globalenv())
  suppressMessages(trace(fn, where = asNamespace("RDS"), print = FALSE,
        exit = quote(if (exists("deg.out.map"))
          assign(".dom_capture", deg.out.map, globalenv()))))
  on.exit({
    suppressMessages(untrace(fn, where = asNamespace("RDS")))
    if (exists(".dom_capture", envir = globalenv()))
      rm(".dom_capture", envir = globalenv())
  })
  suppressWarnings(run())
  dom <- get(".dom_capture", envir = globalenv())
  sum(vapply(dom, function(l) sum(lengths(l)), numeric(1)))
}

test_that("Gile continuous bootstrap keeps every subject in deg.out.map (7.2)", {
  data(fauxmadrona)
  d <- fauxmadrona
  set.seed(1); d$cts <- round(rnorm(nrow(d)), 3)   # continuous; varied degrees
  total <- dom_total("sppsboot4ppsreal", function() {
    set.seed(7)
    RDS.bootstrap.intervals(d, "cts", weight.type = "RDS-II", uncertainty = "Gile",
      N = 1000, number.of.bootstrap.samples = 2,
      number.ss.samples.per.iteration = 100, verbose = FALSE)
  })
  expect_equal(total, nrow(d))
})

test_that("HCG continuous bootstrap keeps every subject in deg.out.map (7.2)", {
  data(fauxmadrona)
  d <- fauxmadrona
  set.seed(1); d$cts <- round(rnorm(nrow(d)), 3)
  total <- dom_total("HCG.bootstrap", function() {
    set.seed(7)
    RDS.bootstrap.intervals(d, "cts", weight.type = "HCG", uncertainty = "HCG",
      N = 1000, number.of.bootstrap.samples = 2, verbose = FALSE)
  })
  expect_equal(total, nrow(d))   # no NA degrees in fauxmadrona
})
