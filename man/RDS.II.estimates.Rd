\name{RDS.II.estimates}
\alias{RDS.II.estimates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ RDS-II Estimates}
\description{
This function computes the RDS-II estimates for a categorical variable or the RDS-II estimate for a 
numeric variable.
}
\usage{
RDS.II.estimates(rds.data, trait.variable, network.variable)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rds.data}{A data frame.  This data frame must identify recruitment patters by a pair of fields named ``recruitment.id'' and ``recruiter.id''.}
  \item{trait.variable}{A categorical or numeric variable to be analyzed.}
  \item{network.variable}{A string giving the name of the variable in the \code{rds.data} that contains the network sizes of 
    survey respondents.}
 }
%\details{}
\value{If \code{trait.variable} is numeric then the RDS-II estimate of the mean is returned, otherwise a vector of proportion estimates is returned.}
\references{
Gile, K. J., Handcock, M. S., 2009b. \emph{Respondent-driven sampling: An assessment 
of current methodology}. Under review, Nuffield College, University of Oxford.

Neely, W. W., 2009. \emph{Bayesian methods for data from respondent driven sampling}. Dissertation in-progress, Department of Statistics, University of Wisconsin, Madison.

Salganik, M., Heckathorn, D. D., 2004. \emph{Sampling and estimation in hidden populations using respondent-driven sampling}. Sociological Methodology 34, 193-239.

Volz, E., Heckathorn, D., 2008. \emph{Probability based estimation theory for Respondent Driven Sampling}. The Journal of Official Statistics 24 (1), 79-97.
}
\author{W. Whipple Neely}
%\note{}
\seealso{
\itemize{
   \item{RDS.I.estimates}{\code{\link{RDS.I.estimates}}}
   \item{RDS.I.estimates}{\code{\link{RDS.I.DS.estimates}}}
   }
}
\examples{
data(simulated)
RDS.II.estimates(rds.data=simulated.data,trait.variable='X',network.variable='network.size')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survey}
\keyword{manip}
