\name{compute.weights}
\alias{compute.weights}
\title{Compute estimates of the sampling weights of the respondent's observations
based on various estimators}
\usage{
  compute.weights(rds.data,
    weight.type = c("Gile's SS", "RDS-I", "RDS-I (DS)", "RDS-II", "Arithmetic Mean"),
    N = NULL, subset = NULL, ...)
}
\arguments{
  \item{rds.data}{An \code{rds.data.frame} that indicates
  recruitment patterns by a pair of attributes named ``id''
  and ``recruiter.id''.}

  \item{weight.type}{A string giving the type of estimator
  to use. The options are \code{"Gile's SS"},
  \code{"RDS-I"}, \code{"RDS-II"}, \code{"RDS-I/DS"}, and
  \code{"Arithemic Mean"}. It defaults to \code{"Gile's
  SS"}.}

  \item{N}{An estimate of the number of members of the
  population being sampled. If \code{NULL} it is read as
  the \code{population.size.mid} attribute of the
  \code{rds.data} frame. If that is missing, the weights
  will sum to 1. Note that this parameter is required for
  Gile's SS.}

  \item{subset}{A logical expression subsetting rds.data.}

  \item{...}{Additional parameters passed to the individual
  weighting algorithms.}
}
\value{
  A vector of weighte for each of the respondents. It is of
  the same size as the number of rows in \code{rds.data}.
}
\description{
  Compute estimates of the sampling weights of the
  respondent's observations based on various estimators
}
\seealso{
  \code{\link{rds.I.weights}},
  \code{\link{gile.ss.weights}}, \code{\link{vh.weights}}
}

