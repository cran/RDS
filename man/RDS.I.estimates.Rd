\name{RDS.I.estimates}
\alias{RDS.I.estimates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute RDS-I Estimates}
\description{
This function computes the basic RDS-I estimates for a categorical variable.
}
\usage{
RDS.I.estimates(rds.data, group.variable, network.variable)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rds.data}{A data frame.  This data frame must identify recruitment patters by a pair of fields named ``recruitment.id'' and ``recruiter.id''.}
  \item{group.variable}{A categorical variable to be analyzed.}
  \item{network.variable}{A string giving the name of the variable in the \code{rds.data} that contains the network sizes of 
    survey respondents.}
}
\details{}
\value{A vector of proportion estimates.}
\references{
Gile, K. J., Handcock, M. S., 2009b. \emph{Respondent-driven sampling: An assessment 
of current methodology}. Under review, Nuffield College, University of Oxford.

Neely, W. W., 2009. \emph{Bayesian methods for data from respondent driven sampling}. Dissertation in-progress, Department of Statistics, University of Wisconsin, Madison.

Salganik, M., Heckathorn, D. D., 2004. \emph{Sampling and estimation in hidden populations using respondent-driven sampling}. Sociological Methodology 34, 193-239.

Volz, E., Heckathorn, D., 2008. \emph{Probability based estimation theory for Respondent Driven Sampling}. The Journal of Official Statistics 24 (1), 79-97.
}
\author{W. Whipple Neely}
\note{This is an initial version of the software.  A more developed version of the package will be released during Summer of 2009.}
\seealso{
\itemize{
   \item{RDS.I.DS.estimates}{\code{\link{RDS.I.DS.estimates}}}
   \item{RDS.II.estimates}{\code{\link{RDS.II.estimates}}}
   }
}
\examples{
data(simulated)
RDS.I.estimates(rds.data=simulated.data,group.variable='X',network.variable='network.size')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey}
\keyword{manip}
