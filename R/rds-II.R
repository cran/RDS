
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
#' @param subset An optional criterion to subset \code{rds.data} by. It is
#' an R expression which, when evaluated, subset the
#' data. In plain English, it can be something like \code{subset = seed > 0} to
#' exclude seeds. It can also be the name of a logical vector of the same length of
#' the outcome variable where TRUE means include it in the analysis. If
#' \code{NULL} then no subsetting is done.
#' @param to.factor force variable to be a factor
#' @return 
#' An object of class \code{rds.weighted.estimate}.
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
RDS.II.estimates <- function(rds.data,outcome.variable, subset=NULL,to.factor=FALSE){
  se <- substitute(subset)
  subset <- eval(se,rds.data,parent.frame())
  if(length(outcome.variable) == 1){
    result <- RDS.estimates.local(rds.data,outcome.variable,subset=subset,
                                  weight.type="RDS-II", to.factor=to.factor)
  }
  else{
    result <- lapply(outcome.variable,function(g){
      RDS.estimates.local(rds.data, g, subset,  
                          weight.type="RDS-II", to.factor=to.factor)
    })
    names(result) <- outcome.variable
  }
  return(result)
}


