context("Recruits-per-subject plot keeps every recruiter (3.8)")

# Finding 3.8. The non-stratified "Recruits per subject" plot did
# table(recruiter.id[!seed])[-1] -- but !seed already excludes the seed
# placeholder, so [-1] deleted the first *real* recruiter (fauxmadrona: 253
# recruiters, 252 plotted). The stratified branch omits [-1] and is correct.

test_that("all recruiters appear in the recruits-per-subject plot", {
  skip_if_not_installed("ggplot2")
  data(fauxmadrona)
  rid <- get.rid(fauxmadrona)
  n.recruiters <- length(table(rid[rid != get.seed.rid(fauxmadrona)]))
  p  <- suppressWarnings(plot(fauxmadrona, plot.type = "Recruits per subject"))
  plotted <- sum(suppressWarnings(ggplot2::ggplot_build(p))$data[[1]]$count)
  expect_equal(plotted, n.recruiters)
})

# Finding 3.14(e). plot.rds.data.frame used the deprecated dot-dot notation
# (..n../..y..) and qplot(), which emit ggplot2 deprecation warnings and are on
# a path to removal (would eventually break R CMD check examples). Replaced with
# after_stat()/geom_bar(). Also added match.arg(plot.type) so typos error instead
# of silently returning NULL.
test_that("no ggplot2 deprecation warnings across plot types (3.14e)", {
  skip_if_not_installed("ggplot2")
  data(fauxmadrona)
  types <- c("Network size by wave", "Recruits by wave",
             "Recruits per seed", "Recruits per subject")
  build_warnings <- function(t, strat) {
    ws <- character(0)
    withCallingHandlers(
      ggplot2::ggplot_build(plot(fauxmadrona, plot.type = t, stratify.by = strat)),
      warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") }
    )
    ws
  }
  for (t in types) {
    expect_false(any(grepl("deprecated", build_warnings(t, NULL), ignore.case = TRUE)),
                 info = paste("unstratified:", t))
    expect_false(any(grepl("deprecated", build_warnings(t, "disease"), ignore.case = TRUE)),
                 info = paste("stratified:", t))
  }
})

test_that("plot.rds.data.frame rejects an unknown plot.type (3.14e)", {
  data(fauxmadrona)
  expect_error(plot(fauxmadrona, plot.type = "nonsense"), "should be one of")
})
