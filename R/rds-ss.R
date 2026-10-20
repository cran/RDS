

#' Gile's SS Estimates
#' 
#' This function computes the sequential sampling (SS) estimates for a
#' categorical variable or numeric variable.
#' 
#' 
#' @param rds.data An \code{rds.data.frame} that indicates recruitment patterns
#' by a pair of attributes named ``id'' and ``recruiter.id''.
#' @param outcome.variable A string giving the name of the variable in the
#' \code{rds.data} that contains a categorical or numeric variable to be
#' analyzed.
#' @param N An estimate of the number of members of the population being
#' sampled. If \code{NULL} it is read as the \code{population.size.mid} attribute of
#' the \code{rds.data} frame. If that is missing no finite population correction is
#' applied and the estimate is the same as RDS-II.
#' @param subset An optional criterion to subset \code{rds.data} by. It is
#' an R expression which, when evaluated, subset the
#' data. In plain English, it can be something like \code{subset = seed > 0} to
#' exclude seeds. It can also be the name of a logical vector of the same length of
#' the outcome variable where TRUE means include it in the analysis. If
#' \code{NULL} then no subsetting is done.
#' @param number.ss.samples.per.iteration The number of samples to take in
#' estimating the inclusion probabilities in each iteration of the sequential
#' sampling algorithm. If \code{NULL} it is read as the
#' eponymous attribute of \code{rds.data}. If that
#' is missing it defaults to 5000.
#' @param number.ss.iterations The number of iterations of the sequential
#' sampling algorithm. If that is missing it defaults to 5.
#' @param control A list of control parameters for algorithm
#' tuning. Constructed using\cr
#' \code{\link{control.rds.estimates}}.
#' @param hajek logical; Use the standard Hajek-type estimator of Gile (2011)
#' or the standard Hortitz-Thompson. The default is TRUE.
#' @param to.factor force variable to be a factor
#' @return 
#' An object of class \code{rds.weighted.estimate}.
#' 
#' @author Krista J. Gile with help from Mark S. Handcock
#' @seealso \code{\link{RDS.I.estimates}}, \code{\link{RDS.II.estimates}}
#' @references
#' 
#' Gile, Krista J. 2011 \emph{Improved Inference for Respondent-Driven Sampling
#' Data with Application to HIV Prevalence Estimation}, \emph{Journal of the
#' American Statistical Association}, 106, 135-146.
#' 
#' Gile, Krista J., Handcock, Mark S., 2010. Respondent-driven Sampling:
#' An Assessment of Current Methodology, Sociological Methodology, 40,
#' 285-327. <doi:10.1111/j.1467-9531.2010.01223.x>
#'
#' Gile, Krista J., Beaudry, Isabelle S. and Handcock, Mark S., 2018 
#' Methods for Inference from Respondent-Driven Sampling Data,
#' Annual Review of Statistics and Its Application
#' <doi:10.1146/annurev-statistics-031017-100704>.
#'
#' Gile, Krista J., Handcock, Mark S., 2015 \emph{Network Model-Assisted
#' Inference from Respondent-Driven Sampling Data}, Journal of the Royal Statistical Society, A. <doi:10.1111/rssa.12091>.
#' 
#' Salganik, M., Heckathorn, D. D., 2004. \emph{Sampling and estimation in
#' hidden populations using respondent-driven sampling}. Sociological
#' Methodology 34, 193-239.
#' 
#' Volz, E., Heckathorn, D., 2008. \emph{Probability based estimation theory
#' for Respondent Driven Sampling}. The Journal of Official Statistics 24 (1),
#' 79-97.
#' @keywords survey manip
#' @examples
#' 
#' data(fauxmadrona)
#' RDS.SS.estimates(rds.data=fauxmadrona,outcome.variable="disease",N=1000)
#' 
#' @export RDS.SS.estimates
RDS.SS.estimates <- function(rds.data,outcome.variable,N=NULL,subset=NULL,
           number.ss.samples.per.iteration=500,number.ss.iterations=5,
           control=control.rds.estimates(),
           hajek=TRUE, to.factor=FALSE){
    se <- substitute(subset)
    subset <- eval(se,rds.data,parent.frame())
    if(length(outcome.variable) == 1){
      result <- RDS.estimates.local(
        rds.data=rds.data,
        outcome.variable=outcome.variable,
        N=N,
        subset=subset,
        number.ss.samples.per.iteration=number.ss.samples.per.iteration,
        number.ss.iterations=number.ss.iterations,
        control=control,
        hajek=hajek,
        weight.type="Gile's SS",
        to.factor=to.factor)
    }
    else {
      result <- lapply(outcome.variable,function(g){
        RDS.estimates.local(
          rds.data=rds.data,
          outcome.variable=g,
          N=N,
          subset=subset,
          number.ss.samples.per.iteration=number.ss.samples.per.iteration,
          number.ss.iterations=number.ss.iterations,
          control=control,
          hajek=hajek,
          weight.type="Gile's SS",
          to.factor=to.factor)
      })
      names(result) <- outcome.variable
    }
    return(result)
  }


