\name{RDS.bootstrap.intervals}
\alias{RDS.bootstrap.intervals}
\title{RDS Bootstrap Interval Estimates}
\usage{
  RDS.bootstrap.intervals(rds.data, outcome.variable,
    weight.type = NULL, uncertainty = NULL, N = NULL,
    subset = NULL, confidence.level = 0.95,
    ties.to.trait = NULL,
    number.of.bootstrap.samples = 500, ...)
}
\arguments{
  \item{rds.data}{An \code{rds.data.frame} that indicates
  recruitment patterns by a pair of attributes named ``id''
  and ``recruiter.id''.}

  \item{outcome.variable}{A string giving the name of the
  variable in the \code{rds.data} that contains a
  categorical or numeric variable to be analyzed.}

  \item{weight.type}{A string giving the type of estimator
  to use. The options are \code{"Gile's SS"},
  \code{"RDS-I"}, \code{"RDS-II"}, \code{"RDS-I (DS)"}, and
  \code{"Arithemic Mean"}. If \code{NULL} it defaults to
  \code{"Gile's SS"}.}

  \item{uncertainty}{A string giving the type of
  uncertainty estimator to use. The options are
  \code{"SRS"}, \code{"Gile"} and \code{"Salganik"}. This
  is usually determined by \code{weight.type} to be
  consistent with the estimator's origins. The estimators
  RDS-I, RDS-I (DS), and RDS-II default to
  \code{"Salganik"}, "Arithmetic Mean" defaults to
  \code{"SRS"} and "Gile's SS" defaults to the
  \code{"Gile"} bootstrap.}

  \item{N}{An estimate of the number of members of the
  population being sampled. If \code{NULL} it is read as
  the \code{population.size.mid} attribute of the
  \code{rds.data} frame. If that is missing it defaults to
  1000.}

  \item{subset}{An optional criterion to subset
  \code{rds.data} by. It is a character string giving an R
  expression which, when evaluated, subset the data. In
  plain English, it can be something like \code{"seed > 0"}
  to exclude seeds. It can be the name of a logical vector
  of the same length of the outcome variable where TRUE
  means include it in the analysis. If \code{NULL} then no
  subsetting is done.}

  \item{confidence.level}{The confidence level for the
  confidence intervals. The default is 0.95 for 95\%.}

  \item{number.of.bootstrap.samples}{The number of
  bootstrap samples to take in estimating the uncertainty
  of the estimator. If \code{NULL} it defaults to 500.}

  \item{ties.to.trait}{A vector of strings giving the names
  of the variables in the \code{rds.data} that contains the
  respondent's report of the number of ties they have to
  population members in each of the trait values. If
  missing, it is computed from the sample as the number of
  ties to recruitees who have each of the trait values.}

  \item{...}{Additional arguments for RDS.*.estimates.}
}
\value{
  An object of class \code{rds.interval.estimate}
  summarizing the inference. The confidence interval and
  standard error are based on the bootstrap procedure. In
  additon, the object has attribute \code{bsresult} which
  provides details of the bootstrap procedure. The contents
  of the \code{bsresult} attribute depends on the
  \code{uncertainty} used. If
  \code{uncertainty=="Salganik"} then \code{bsresult} is a
  vector of standard deviations of the bootstrap samples.
  If \code{uncertainty=="Gile's SS"} then \code{bsresult}
  is a list with components for the bootstrap point
  estimate, the bootstrap samples themselves and the
  standard deviations of the bootstrap samples. If
  \code{uncertainty=="SRS"} then \code{bsresult} is NULL.
}
\description{
  This function computes an interval estimate for one or
  more categorical variables. It optionally uses attributes
  of the RDS data set to determine the type of estimator
  and type of uncertainty estimate to use.
}
\examples{
data(fauxmadrona)
RDS.bootstrap.intervals(rds.data=fauxmadrona,weight.type="RDS-II",uncertainty="Salganik",
	outcome.variable="disease",N=1000,number.of.bootstrap.samples=50)
}
\references{
  Gile, Krista J. 2011 \emph{Improved Inference for
  Respondent-Driven Sampling Data with Application to HIV
  Prevalence Estimation}, \emph{Journal of the American
  Statistical Association}, 106, 135-146.

  Gile, Krista J., Handcock, Mark S., 2010
  \emph{Respondent-driven Sampling: An Assessment of
  Current Methodology}. Sociological Methodology 40,
  285-327.
}
\keyword{manip}
\keyword{survey}

