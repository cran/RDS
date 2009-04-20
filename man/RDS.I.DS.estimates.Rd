\name{RDS.I.DS.estimates}
\alias{RDS.I.DS.estimates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{RDS-I/DS Estimates}
\description{
This function computes the basic RDS-I/DS estimates for a categorical variable.  This is the ``data smoothed'' version of RDS where it is assumed that the observed Markov process is reversible.
}
\usage{
RDS.I.DS.estimates(rds.data, group.variable, network.variable)
}
%- maybe also 'usage' for other objects documented here.
\arguments{  
       \item{rds.data}{A data frame.  This data frame must identify recruitment patters by a pair of fields named          ``recruitment.id'' and ``recruiter.id''.}
  \item{group.variable}{A categorical variable to be analyzed.}
  \item{network.variable}{A string giving the name of the variable in the \code{rds.data} that contains the network sizes of 
    survey respondents.}
}
\details{}
\value{A vector of proportion estimates constructed using the ``smoothed'' transition matrix.}
\references{
Gile, K. J., Handcock, M. S., 2009b. \emph{Respondent-driven sampling: An assessment 
of current methodology}. Under review, Nuffield College, University of Oxford.

Neely, W. W., 2009. \emph{Bayesian methods for data from respondent driven sampling}. Dissertation in-progress, Department of Statistics, University of Wisconsin, Madison.

Salganik, M., Heckathorn, D. D., 2004. \emph{Sampling and estimation in hidden populations using respondent-driven sampling}. Sociological Methodology 34, 193-239.

Volz, E., Heckathorn, D., 2008. \emph{Probability based estimation theory for Respondent Driven Sampling}. The Journal of Official Statistics 24 (1), 79-97.
}
\author{W. Whipple Neely}
\note{}
\seealso{
\itemize{
   \item{RDS.I.estimates}{\code{\link{RDS.I.estimates}}}
   \item{RDS.II.estimates}{\code{\link{RDS.II.estimates}}}
   }
}

\examples{
data(simulated)
RDS.I.DS.estimates(rds.data=simulated.data,group.variable='X',network.variable='network.size')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survey}
\keyword{manip}
