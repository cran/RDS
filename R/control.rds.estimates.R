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
#' @param SS.infinity The sample proportion, \code{n/N}, below which the computation of the SS weights should simplify to that of the \code{RDS-II} weights.
#' @param useC Use a C-level implementation of Gile's bootstrap (rather than
#' the R level). The implementations should be computational 
#' equivalent (except for speed).
#' @param number.of.bootstrap.samples The number of bootstrap samples to take
#' in estimating the uncertainty of the estimator. If \code{NULL} it defaults
#' to the number necessary to compute the standard error to accuracy 0.001.
#' @param seed Seed value (integer) for the random number generator.  See
#' \code{\link[base]{set.seed}}
#' @return A list with arguments as components.
#' @seealso \code{\link{RDS.bootstrap.intervals}}
#' @keywords models
control.rds.estimates<-function(confidence.level=0.95,
                              SS.infinity=0.01,
                              useC=TRUE,
                              number.of.bootstrap.samples = NULL,
                              seed=NULL){

  control<-list()
  formal.args<-formals(sys.function())
  for(arg in names(formal.args))
    control[arg]<-list(get(arg))

  set.control.class()
}
