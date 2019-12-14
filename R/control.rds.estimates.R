utils::globalVariables(c(".control.rds.estimates"))
#' Auxiliary for Controlling RDS.bootstrap.intervals
#'
#' Auxiliary function as user interface for fine-tuning RDS.bootstrap.intervals algorithm,
#' which computes interval estimates for via bootstrapping.
#'
#' This function is only used within a call to the \code{\link{RDS.bootstrap.intervals}}
#' function.
#'
#' Some of the arguments are not yet fully implemented. It will evolve slower to incorporate more
#' arguments as the package develops.
#'
#' @param confidence.level The confidence level for the confidence intervals.
#' The default is 0.95 for 95\%.
#' @param SS.infinity The sample proportion, \code{n/N}, below which the computation
#' of the SS weights should simplify to that of the \code{RDS-II} weights.
#' @param lowprevalence Standard confidence interval procedures can be inaccurate when the
#' outcome expected count is close to zero. This sets conditions where alternatives to the
#' standard are used for the \code{ci.type="hmg"} option. See Details for its use.
#' @param discrete.cutoff The minimum proportion of the values of the outcome variable that
#' need to be unique before the variable is judged to be continuous.
#' @param useC Use a C-level implementation of Gile's bootstrap (rather than
#' the R level). The implementations should be computational
#' equivalent (except for speed).
#' @param number.of.bootstrap.samples The number of bootstrap samples to take
#' in estimating the uncertainty of the estimator. If \code{NULL} it defaults
#' to the number necessary to compute the standard error to accuracy 0.001.
#' @param hcg.reltol Relative convergence tolerance for the HCG estimator.  The algorithm stops if
#' it is unable to reduce the log-likelihood by a factor of \code{reltol * (abs(log-likelihood) + reltol)}
#' at a step. Defaults to \code{sqrt(.Machine$double.eps)}, typically about \code{1e-8}.
#' @param hcg.BS.reltol Relative convergence tolerance for the bootstrap of the HCG estimator. 
#' It has the same interpretation as \code{hcg.reltol} except it is applied to each bootstrap sample.
#' It is typically the same or larger than \code{hcg.reltol}.
#' @param hcg.max.optim The number of iterations on the likelihood optimization for the HCG estimator.
#' @param seed Seed value (integer) for the random number generator.  See
#' \code{\link[base]{set.seed}}
#' @return A list with arguments as components.
#' @details Standard confidence interval procedures can be inaccurate when the
#' outcome expected count is close to zero. In these cases
#' the combined Agresti-Coull and the bootstrap-t interval of
#' Mantalos and Zografos (2008) can be used.
#' The \code{lowprevalence} argument is a
#' two vector parameter setting the conditions under which the approximation is used.
#' The first is the penalty term on the differential activity. If the observed number
#' of the rare group minus the product of the first parameter and the differential
#' activity is lower than the
#' second parameter, the low prevalence approximation is used.
#' @seealso \code{\link{RDS.bootstrap.intervals}}
#' @keywords models
#' @export
control.rds.estimates <- function(confidence.level = 0.95,
                                  SS.infinity = 0.01,
                                  lowprevalence = c(8, 14),
                                  discrete.cutoff = 0.8,
                                  useC = TRUE,
                                  number.of.bootstrap.samples = NULL,
                                  hcg.reltol=sqrt(.Machine$double.eps),
                                  hcg.BS.reltol=100000*sqrt(.Machine$double.eps),
                                  hcg.max.optim=500,
                                  seed = NULL) {
  formal.args <- formals(sys.function())
  if (!exists(".control.rds.estimates")) {
    control <- list()
    for (arg in names(formal.args))
      control[arg] <- list(get(arg))
  } else{
    control <- .control.rds.estimates
    if (!missing(confidence.level)) {
      control[["confidence.level"]] <- confidence.level
    }
    if (!missing(SS.infinity)) {
      control[["SS.infinity"]] <- SS.infinity
    }
    if (!missing(lowprevalence)) {
      control[["lowprevalence"]] <- lowprevalence
    }
    if (!missing(discrete.cutoff)) {
      control[["discrete.cutoff"]] <- discrete.cutoff
    }
    if (!missing(useC)) {
      control[["useC"]] <- useC
    }
    if (!missing(number.of.bootstrap.samples)) {
      control[["number.of.bootstrap.samples"]] <-
        number.of.bootstrap.samples
    }
    if (!missing(hcg.reltol)) {
      control[["hcg.reltol"]] <- hcg.reltol
    }
    if (!missing(hcg.BS.reltol)) {
      control[["hcg.BS.reltol"]] <- hcg.BS.reltol
    }
    if (!missing(hcg.max.optim)) {
      control[["hcg.max.optim"]] <- hcg.max.optim
    }
    if (!missing(seed)) {
      control[["seed"]] <- seed
    }
  }
  
  RDS::set.control.class("control.rds.estimates")
}
