\name{RDS.compare.proportions}
\alias{RDS.compare.proportions}
\title{Compares the rates of two variables against one another.}
\usage{
  RDS.compare.proportions(first.interval, second.interval,
    M = 10000)
}
\arguments{
  \item{first.interval}{An \code{rds.interval.estimate}
  object fit with either "Gile" or "Salganik" uncertainty.}

  \item{second.interval}{An \code{rds.interval.estimate}
  object fit with either "Gile" or "Salganik" uncertainty.}

  \item{M}{The number of bootstrap resamplings to use}
}
\description{
  Compares the rates of two variables against one another.
}
\details{
  This function preforms a bootstrap test comparing the
  data(faux) int1 <- RDS.bootstrap.intervals(faux,
  outcome.variable=c("X"), weight.type="RDS-II",
  uncertainty="Salganik", N=1000,
  number.ss.samples.per.iteration=1000,
  confidence.level=0.95, number.of.bootstrap.samples=500)
  int2 <- RDS.bootstrap.intervals(faux,
  outcome.variable=c("Y"), weight.type="RDS-II",
  uncertainty="Salganik", N=1000,
  number.ss.samples.per.iteration=1000,
  confidence.level=0.95, number.of.bootstrap.samples=500)
  RDS.compare.proportions(int1,int2)
}

