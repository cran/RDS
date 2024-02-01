#' MA Estimates
#'
#' This function computes the sequential sampling (MA) estimates for a
#' categorical variable or numeric variable.
#'
#'
#' @param rds.data An \code{rds.data.frame} that indicates recruitment patterns
#' by a pair of attributes named ``id'' and ``recruiter.id''.
#' @param trait.variable A string giving the name of the variable in the
#' \code{rds.data} that contains a categorical or numeric variable to be
#' analyzed.
#' @param seed.selection An estimate of the mechanism guiding the choice of
#' seeds.  The choices are \describe{ \item{"allwithtrait"}{indicating that all
#' the seeds had the trait;} \item{"random"}{meaning they were, as if, a simple
#' random sample of individuals from the population;}
#' \item{"sample"}{indicating that the seeds are taken as those in the sample
#' (and resampled for the population with that composition if necessary);}
#' \item{"degree"}{is proportional to the degree of the individual;}
#' \item{"allwithtraitdegree"}{indicating that all the seeds had the trait and
#' the probability of being a seed is proportional to the degree of the
#' respondent.} }
#' @param number.of.seeds The number of seeds chosen to initiate the sampling.
#' @param number.of.coupons The number of coupons given to each respondent.
#' @param number.of.iterations The number of iterations used at the core of the
#' algorithm.
#' @param N An estimate of the number of members of the population being
#' sampled. If \code{NULL} it is read as the \code{pop.size.mid} attribute of
#' the \code{rds.data} frame. If that is missing it defaults to 1000.
#' @param M1 The number of networked populations generated at each iteration.
#' @param M2 The number of (full) RDS samples generated for each networked
#' population at each iteration.
#' @param seed The random number seed used to initiate the computations.
#' @param initial.sampling.probabilities Initialize sampling probabilities for
#' the algorithm. If missing, they are taken as proportional to degree, and
#' this is almost always the best starting values.
#' @param MPLE.samplesize Number of samples to take in the computation of the
#' maximum pseudolikelihood estimator (MPLE) of the working model parameter.
#' The default is almost always sufficient.
#' @param SAN.nsteps Number of MCMC proposals for all the annealing runs combined.
#' @param sim.interval Number of MCMC steps between each of the M1 sampled
#' networks per iteration.
#' @param number.of.cross.ties The expected number of ties between those with
#' the trait and those without. If missing, it is computed based on the
#' respondent's reports of the number of ties they have to population members
#' who have the trait (i.e. \code{ties.to.trait.variable}) and do not have the
#' trait (i.e. \code{ties.not.to.trait.variable}).
#' @param parallel Number of processors to use in the computations. The default
#' is 1, that is no parallel processing.
#' @param parallel.type The type of cluster to start. e.g. 'PSOCK', 'MPI', etc.
#' @param max.degree Impose ceiling on degree size.
#' @param verbose Should verbose diagnostics be printed while the algorithm is
#' running.
#' @param full.output More verbose output
#' @param SAN.maxit A ceiling on the number of simulated annealing iterations.
#' @return If \code{trait.variable} is numeric then the model-assisted estimate
#' of the mean is returned, otherwise a vector of proportion estimates is
#' returned. If \code{full.output=TRUE} this leads to:
#'
#' If \code{full.output=FALSE} this leads to an object of class
#' \code{rds.interval.estimate} which is a list with component
#' \describe{
#'   \item{estimate}{the numerical point estimate of proportion of the\code{trait.variable}.}
#'   \item{interval}{a matrix with size columns and one row per category of \code{trait.variable}:
#'     \describe{
#'       \item{point estimate}{ The HT estimate of the population mean.}
#'       \item{95\% Lower Bound}{Lower 95\% confidence bound}
#'       \item{95\% Upper Bound}{Upper 95\% confidence bound}
#'       \item{Design Effect}{The design effect of the RDS}
#'       \item{s.e.}{standard error}
#'       \item{n}{count of the number of sample values with that value of the trait}
#'     }
#'   }
#'   \item{rds.data}{an \code{rds.data.frame} that indicates recruitment
#'      patterns by a pair of attributes named ``id'' and ``recruiter.id''.}
#'   \item{N}{an estimate of the number of members of the population being
#'     sampled. If \code{NULL} it is read as the \code{pop.size.mid} attribute of
#'     the \code{rds.data} frame. If that is missing it defaults to 1000.}
#'   \item{M1}{the number of networked populations generated at each iteration.}
#'   \item{M2}{the number of (full) RDS populations generated for each
#'     networked population at each iteration.}
#'   \item{seed}{the random number seed used to initiate the computations.}
#'   \item{seed.selection}{an estimate of the mechanism guiding the choice
#'     of seeds.  The choices are 
#'     \describe{
#'       \item{"allwithtrait"}{indicating that
#'         all the seeds had the trait;} \item{"random"}{meaning they were, as if, a
#'         simple random sample of individuals from the population;}
#'       \item{"sample"}{indicating that the seeds are taken as those in the sample
#'         (and resampled for the population with that composition if necessary);}
#'       \item{"degree"}{is proportional to the degree of the individual;}
#'       \item{"allwithtraitdegree"}{indicating that all the seeds had the trait and
#'         the probability of being a seed is proportional to the degree of the
#'         respondent.}
#'     }
#'   }
#'   \item{number.of.seeds}{The number of seeds chosen to initiate the sampling.}
#'   \item{number.of.coupons}{The number of coupons given to each respondent.}
#'   \item{number.of.iterations}{The number of iterations used at the core of the algorithm.}
#'   \item{outcome.variable}{The name of the outcome variable}
#'   \item{weight.type}{The type of weighting used (i.e. MA)}
#'   \item{uncertainty}{The type of weighting used (i.e. MA)}
#'   \item{details}{A list of other diagnostic output from the computations.}
#'   \item{varestBS}{Output from the bootstrap procedure. A list with two
#'     elements: \code{var} is the bootstrap variance, and \code{BSest} is the
#'     vector of bootstrap estimates themselves.}
#'   \item{coefficient}{estimate of the parameter of the ERGM for the
#'     network.}
#' }
#' 
#' @author Krista J. Gile with help from Mark S. Handcock
#' @seealso \code{\link{RDS.I.estimates}}, \code{\link{RDS.I.estimates}}
#' @references Gile, Krista J. 2011 Improved Inference for
#' Respondent-Driven Sampling Data with Application to HIV Prevalence
#' Estimation, Journal of the American Statistical Association, 106,
#' 135-146.
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
#' %Gile, Krista J., Handcock, Mark S., 2011 Network Model-Assisted
#' %Inference from Respondent-Driven Sampling Data, ArXiv Preprint.
#'
#' %%Neely, W. W., 2009. \emph{Bayesian methods for data from respondent driven
#' %sampling}. Dissertation in-progress, Department of Statistics, University of
#' %Wisconsin, Madison.
#'
#' %Salganik, M., Heckathorn, D. D., 2004. Sampling and estimation in
#' %hidden populations using respondent-driven sampling. Sociological
#' %Methodology 34, 193-239.
#'
#' %Volz, E., Heckathorn, D., 2008. Probability based estimation theory
#' %for Respondent Driven Sampling. The Journal of Official Statistics 24 (1),
#' %79-97.
#' @keywords survey manip
#' @examples
#'
#' \dontrun{
#' data(faux)
#' MA.estimates(rds.data=faux,trait.variable='X')
#' }
#'
#' @export MA.estimates
#' @importFrom ergm summary_formula
MA.estimates <- function(rds.data,
                         trait.variable,
                         seed.selection = "degree",
                         number.of.seeds = NULL,
                         number.of.coupons = NULL,
                         number.of.iterations = 3,
                         N = NULL,
                         M1 = 25,
                         M2 = 20,
                         seed = 1,
                         initial.sampling.probabilities = NULL,
                         MPLE.samplesize = 50000,
                         SAN.maxit = 5,
                         SAN.nsteps = 2^19,
                         sim.interval = 10000,
                         number.of.cross.ties = NULL,
                         max.degree = NULL,
                         parallel = 1,
                         parallel.type = "PSOCK",
                         full.output = FALSE,
                         verbose = TRUE)
{
  if (!(trait.variable %in% names(rds.data))) {
    stop(sprintf("No variable called %s appears in the data.",
                 trait.variable))
  }
  
  network.size <- attr(rds.data, "network.size.variable")
  
  # Deal with missing values. This needs to be improved!
  
  remvalues <- !is.na(rds.data[[network.size]])
  if (sum(remvalues) < nrow(rds.data)) {
    warning(
      paste(
        nrow(rds.data) - sum(remvalues),
        "of",
        nrow(rds.data),
        "values were missing their network size values and were devalued"
      ),
      call. = FALSE
    )
    rds.data[[network.size]][!remvalues] <-
      max(rds.data[[network.size]], na.rm = TRUE)
  }
  
  rds.data[[trait.variable]] <-
    as.vector(rds.data[[trait.variable]])
  remvalues <- !is.na(rds.data[[trait.variable]])
  if (sum(remvalues) < nrow(rds.data)) {
    warning(
      paste(
        nrow(rds.data) - sum(remvalues),
        "of",
        nrow(rds.data),
        "outcome values were missing and were set to the null value."
      ),
      call. = FALSE
    )
    rds.data[[trait.variable]][!remvalues] <-
      min(rds.data[[network.size]], na.rm = TRUE)
  }
  
  if (is.null(N)) {
    N <- attr(rds.data, "population.size.mid")
    if (is.null(N)) {
      N <- ceiling(nrow(rds.data) / 0.04)
      warning(paste(
        "Parameter N missing, with no default for this data set. Using N =",
        N,
        "\n"
      ))
    } else{
      cat("\nNote: Using the data's mid population size estimate: N =",
          N,
          "\n")
    }
  }
  
  remvalues <-
    rds.data[[network.size]] == 0 | is.na(rds.data[[network.size]])
  if (any(remvalues)) {
    warning(
      paste(
        sum(remvalues),
        "of",
        nrow(rds.data),
        "network sizes were zero. The estimator will presume these are",
        max(rds.data[[network.size]], na.rm = TRUE)
      ),
      call. = FALSE
    )
    rds.data[[network.size]][remvalues] <-
      max(rds.data[[network.size]], na.rm = TRUE)
  }
  
  recruiter.id <-
    as.character(rds.data[[attr(rds.data, "recruiter.id")]])
  outcome <- as.vector(rds.data[[trait.variable]])
  toutcome <- table(outcome)
  maxa <- which.max(toutcome)
  toutcome[maxa] <- 0
  maxa <- c(maxa, which.max(toutcome))
  onames <- sort(names(toutcome)[maxa])
  # So only binary outcomes are dealth with
  outcome <- as.numeric(outcome == max(onames, na.rm = TRUE))
  outcome[is.na(outcome)] <- 0
  deg <- rds.data[[network.size]]
  if (is.null(max.degree)) {
    max.degree <- max(deg, na.rm = TRUE)
  }
  deg[deg > max.degree] <- max.degree
  recruits <- get.number.of.recruits(rds.data)
  
  if (is.null(number.of.coupons)) {
    number.of.coupons <- attr(rds.data, "max.coupons")
    if (is.null(number.of.coupons)) {
      number.of.coupons <- max(recruits)
      warning(
        paste(
          "Parameter number.of.coupons missing, with no
          default for this data set. Using number.of.coupons =",
          number.of.coupons,
          ", the maximum number in the data set\n"
        )
      )
    } else{
      cat(
        "\nNote: Using the data's maximum number of coupons as the number.of.coupons =",
        number.of.coupons,
        "\n"
      )
    }
  }
  

  # modify the vector of numbers of referrals to diseased nodes
  # to make an estimate of the number of ties to diseased in the pop.
  todis <- rep(0, nrow(rds.data))
  b = tapply(outcome, recruiter.id, sum, na.rm = TRUE)[-1]
  d = match(names(b), as.character(rds.data[[attr(rds.data, "id")]]))
  todis[d[!is.na(d)]] <- b[!is.na(d)]
  total.to.dis <- sum(todis)
  todis[recruits > 0] <-
    deg[recruits > 0] * todis[recruits > 0] / recruits[recruits > 0]
  todis[recruits == 0] <-
    deg[recruits == 0] * total.to.dis / sum(recruits)
  # modify the vector of number of referrals to undiseased
  # to make an estimate of the number of ties to undiseased in the pop.
  tonodis <- rep(0, nrow(rds.data))
  b = tapply(1 - outcome, recruiter.id, sum, na.rm = TRUE)[-1]
  d = match(names(b), as.character(rds.data[[attr(rds.data, "id")]]))
  tonodis[d[!is.na(d)]] <- b[!is.na(d)]
  total.not.to.dis <- sum(tonodis)
  tonodis[recruits > 0] <-
    deg[recruits > 0] * tonodis[recruits > 0] / recruits[recruits > 0]
  tonodis[recruits == 0] <-
    deg[recruits == 0] * total.not.to.dis / sum(recruits)

  
  
  wave <- get.wave(rds.data)
  
  if (is.null(number.of.seeds)) {
    number.of.seeds <- sum(wave == 0, na.rm = TRUE)
  }
  
  fixinitial <- switch(
    seed.selection,
    "sample" = -2,
    "random" = -3,
    "allwithtrait" = 1,
    "degree" = 3,
    "allwithtraitdegree" = 2,
    1
  )
  if (fixinitial == -2) {
    attr(fixinitial, "table.seeds") <- table(outcome[wave == 0])
  }
  if (N < length(outcome)) {
    stop(
      sprintf(
        "The population size, %d, is less than the sample
        size, %d. The population size must be at least as large as the
        sample size for the estimate to make sense.",
        N,
        length(outcome)
      )
    )
  }
  if (is.null(max.degree)) {
    max.degree <- max(deg)
  }
  out <- getnetest(
    dissample = outcome,
    degsample = deg,
    N = N,
    todisall = todis,
    tonodisall = tonodis,
    nit = number.of.iterations,
    M1 = M1,
    M2 = M2,
    seed = seed,
    short = FALSE,
    theta0 = -1,
    parallel = parallel,
    parallel.type = parallel.type,
    Qstart = initial.sampling.probabilities,
    MPLEsamplesize = MPLE.samplesize,
    maxdeg = max.degree,
    SAN.maxit = SAN.maxit,
    SAN.nsteps = SAN.nsteps,
    fixinitial = fixinitial,
    nsamp0 = number.of.seeds,
    coupons = number.of.coupons,
    crossStat = number.of.cross.ties,
    interval = sim.interval,
    net = NULL,
    rds.data = rds.data,
    trait.variable = trait.variable,
    verbose = verbose
  )

  # The next line is bias correction
  outprev <- out$prev
  estimate <- c(out$prev, 0, 0)
  a <- out$prev + c(-1.96, 1.96) * sqrt(out$varest)
  #
  a[a < sum(outcome) / N] <- sum(outcome) / N
  a[a > (1 - sum(1 - outcome) / N)] <- 1 - sum(1 - outcome) / N
  estimate[2:3] <- a
  estimate <- as.numeric(estimate)
  names(estimate)[1] <- trait.variable
  nsamplesbyoutcome <- table(as.numeric(outcome))
  #
  # Design effect
  #
  varoutcome <- estimate[1] * (1 - estimate[1])
  #    Note the finite sample correction factor
  varsrs <-
    (((N - length(outcome)) / (N - 1)) * varoutcome / length(outcome))
  de <- ((estimate[3] - estimate[2]) / (1.96 * 2)) ^ 2 / varsrs
  estimate <-
    c(estimate,
      de,
      (estimate[3] - estimate[2]) / (1.96 * 2),
      nsamplesbyoutcome[2])
  names(estimate)[c(4, 5, 6)] <- c("Design Effect", "s.e.", "n")
  estimate <-
    rbind(c((1 - estimate)[c(1, 3, 2)], estimate[4:6]), estimate)
  estimate[1, 6] <- nsamplesbyoutcome[1]
  rownames(estimate) <- onames
  estimate <- as.numeric(estimate)
  names(estimate) <- NULL
  names(estimate)[1:2] <- onames
  estimate <-
    rds.interval.estimate(
      estimate,
      trait.variable,
      weight.type = "MA",
      uncertainty = "MA",
      weights = out$weights
    )
  
  if (full.output) {
    result <- list(
      rds.data,
      N = N,
      M1 = M1,
      M2 = M2,
      seed = seed,
      seed.selection = seed.selection,
      number.of.seeds = number.of.seeds,
      number.of.coupons = number.of.coupons,
      number.of.iterations = number.of.iterations,
      estimate = estimate,
      details = out$istuff,
      varestBS = out$varestBS,
      varest = out$varest,
      coefficient = out$coef
    )
  } else{
    result <- estimate
  }
  return(result)
}
