\name{cumulative.estimate}
\alias{cumulative.estimate}
\title{Calculates estimates at each successive wave of the sampling process}
\usage{
  cumulative.estimate(rds.data, outcome.variable,
    est.func = RDS.II.estimates, ...)
}
\arguments{
  \item{rds.data}{An rds.data.frame}

  \item{outcome.variable}{The outcome}

  \item{est.func}{A function taking rds.data and
  outcome.variable as parameters and returning an
  rds.weighted.estimate object}

  \item{...}{additional parameters for est.func}
}
\description{
  Calculates estimates at each successive wave of the
  sampling process
}

