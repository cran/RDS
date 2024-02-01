
#' RDS-II Estimates
#' 
#' This function computes the RDS-II estimates for a categorical variable or
#' the RDS-II estimate for a numeric variable.
#' 
#' 
#' @param rds.data An \code{rds.data.frame} that indicates recruitment patterns
#' by a pair of attributes named ``id'' and ``recruiter.id''.
#' @param outcome.variable A string giving the name of the variable in the
#' \code{rds.data} that contains a categorical or numeric variable to be
#' analyzed.
#' @param empir.lik If true, and outcome.variable is numeric, standard errors
#' based on empirical likelihood will be given.
#' @param subset An optional criterion to subset \code{rds.data} by. It is
#' an R expression which, when evaluated, subset the
#' data. In plain English, it can be something like \code{subset = seed > 0} to
#' exclude seeds. It can also be the name of a logical vector of the same length of
#' the outcome variable where TRUE means include it in the analysis. If
#' \code{NULL} then no subsetting is done.
#' @param N Population size to be used to calculate the empirical likelihood interval. If NULL, this value is
#' taken to be the population.size.mid attribute of the data and if that is not set, no finite population
#' correction is used.
#' @param to.factor force variable to be a factor
#' @return 
#' If \code{outcome.variable} is numeric then the RDS-II estimate of the mean is returned, otherwise a vector of proportion estimates is returned.
#' If the \code{empir.lik} is true, an object of class \code{rds.interval.estimate} is returned. This is a list with components
#' \itemize{ \item\code{estimate}: The numerical point estimate of proportion
#' of the \code{trait.variable}.  \item\code{interval}: A matrix with six
#' columns and one row per category of \code{trait.variable}: \itemize{
#' \item\code{point estimate}: The HT estimate of the population mean.
#' \item\code{95\% Lower Bound}: Lower 95\% confidence bound.  \item\code{95\%
#' Upper Bound}: Upper 95\% confidence bound.  \item\code{Design Effect}: The
#' design effect of the RDS.  \item\code{s.e.}: Standard error.  \item\code{n}:
#' Count of the number of sample values with that value of the trait.  } }
#' 
#' Otherwise, an object of class \code{rds.II.estimate} is returned.
#' 
#' @author Mark S. Handcock and W. Whipple Neely
#' @seealso \code{\link{RDS.I.estimates}}, \code{\link{RDS.SS.estimates}}
#' @references Gile, Krista J., Handcock, Mark S., 2010. Respondent-driven Sampling:
#' An Assessment of Current Methodology, Sociological Methodology, 40,
#' 285-327. <doi:10.1111/j.1467-9531.2010.01223.x>
#'
#' Gile, Krista J., Beaudry, Isabelle S. and Handcock, Mark S., 2018 
#' Methods for Inference from Respondent-Driven Sampling Data,
#' Annual Review of Statistics and Its Application
#' <doi:10.1146/annurev-statistics-031017-100704>.
#' 
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
#' data(faux)
#' RDS.II.estimates(rds.data=faux,outcome.variable='X')
#' RDS.II.estimates(rds.data=faux,outcome.variable='X',subset= Y!="blue")
#' 
#' @export RDS.II.estimates
RDS.II.estimates <- function(rds.data,outcome.variable, N=NULL,subset=NULL,empir.lik=TRUE,to.factor=FALSE){
  se <- substitute(subset)
  subset <- eval(se,rds.data,parent.frame())
  if(length(outcome.variable) == 1){
    result <- RDS.estimates.local(rds.data,outcome.variable,subset=subset,
                                  empir.lik=empir.lik, weight.type="RDS-II", N=N, to.factor=to.factor)
  }
  else{
    result <- lapply(outcome.variable,function(g){
      RDS.estimates.local(rds.data, g, subset, empir.lik=empir.lik, 
                          weight.type="RDS-II", N=N, to.factor=to.factor)
    })
    names(result) <- outcome.variable
  }
  return(result)
}


