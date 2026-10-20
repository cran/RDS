
.bootstrap.ci <- function(estimate, observed.estimate, bsests, ci.type, confidence.level) {
  # estimate: point/lower/upper matrix, rows aligned to observed.estimate (outclasses order)
  # bsests:   replicates x groups, columns aligned the same way
  qs <- function() t(apply(bsests, 2, stats::quantile,
                           c((confidence.level + 1)/2, (1 - confidence.level)/2), na.rm = TRUE))
  estimate[, 1] <- observed.estimate
  switch(ci.type,
         pivotal    = { estimate[, 2:3] <- 2 * observed.estimate - qs(); estimate },
         quantile   = { estimate[, 2:3] <- qs()[, 2:1];                   estimate },
         proportion = {
           if (all(observed.estimate >= 0) &&
               (any(observed.estimate < 0.15) || any(estimate[, 2:3] < 0 | estimate[, 2:3] > 1)))
             estimate[, 2:3] <- qs()[, 2:1]
           estimate[, 2:3][estimate[, 2:3] < 0] <- 0
           estimate[, 2:3][estimate[, 2:3] > 1] <- 1
           estimate
         },
         estimate  # "t" default: keep the SE-based lower/upper
  )
}



RDS.bootstrap.intervals.local <- function(rds.data,
                                          outcome.variable,
                                          weight.type,
                                          uncertainty,
                                          N,
                                          subset,
                                          number.of.bootstrap.samples,
                                          confidence.level = .95,
                                          control = control.rds.estimates(),
                                          fast = FALSE,
                                          useC = FALSE,
                                          csubset = "",
                                          ci.type = "t",
                                          to.factor=FALSE, 
                                          cont.breaks = 3,
                                          ...) {
  if (is(rds.data, "rds.data.frame")) {
    if (!(outcome.variable %in% names(rds.data))) {
      stop(sprintf(
        "No variable called %s appears in the data.",
        outcome.variable
      ))
    }
    network.size <- attr(rds.data, "network.size.variable")
    if (!(network.size %in% names(rds.data)))
      stop("invalid network size variable")
    rds.data[[network.size]] <- as.numeric(rds.data[[network.size]])
  }
  else {
    stop("rds.data must be of type rds.data.frame")
  }
  
  if (is.null(weight.type)) {
    weight.type <- "Gile's SS"
  }
  weight.type <- match.arg(
    weight.type,
    c(
      "Gile's SS",
      "RDS-I",
      "RDS-II",
      "RDS-I (DS)",
      "HCG",
      "Arithmetic Mean"
    )
  )
  if (is.na(weight.type)) {
    # User typed an unrecognizable name
    stop(
      paste(
        'You must specify a valid weight.type. The valid types are "Gile\'s SS","RDS-I", "RDS-II", "RDS-I (DS)", "HCG", and "Arithmetic Mean"'
      ),
      call. = FALSE
    )
  }
  
  if (is.null(uncertainty)) {
    if (weight.type == "Gile's SS") {
      uncertainty <- "Gile"
    } else if (weight.type %in% c("RDS-I", "RDS-I (DS)", "RDS-II")) {
      uncertainty <- "Salganik"
    } else if (weight.type == "Arithmetic Mean") {
      uncertainty <- "SRS"
    } else if (weight.type == "HCG"){
      uncertainty <- "HCG"
    }
  }
  uncertainty <-
    match.arg(uncertainty, c("Treeboot", "Gile", "Salganik", "SRS", "HCG"))
  if (is.null(N)) {
    N <- attr(rds.data, "population.size.mid")
    if (is.null(N) && uncertainty == "Gile") {
      stop(
        "N not specified with no default in the data. N is required for Gile's SS estimator and bootstrap"
      )
    }
    if (!is.null(N))
      cat("\nNote: Using the data's mid population size estimate: N =",
          N,
          "\n")
  }
  
  ###########################################################################################
  # Check for missing values and warn the user if any are removed.   This should really taken
  # care of elsewhere.  NB: It is also worth considering the semantics of the message
  # "values were missing and were removed".
  remvalues <-
    rds.data[[network.size]] == 0 | is.na(rds.data[[network.size]])
  if (any(remvalues)) {
    warning(
      paste(
        sum(remvalues),
        "of",
        nrow(rds.data),
        "network sizes were missing or zero. The estimator will presume these are",
        max(rds.data[[network.size]], na.rm = TRUE)
      ),
      call. = FALSE
    )
    
    rds.data[[network.size]][remvalues] <-
      max(rds.data[[network.size]], na.rm = TRUE)
  }
  rds.data.nomiss <- rds.data

  se <- substitute(subset)
  subset <- eval(se, rds.data, parent.frame())
  if (is.null(se) | is.null(subset)) {
    subset <- rep(TRUE, length = nrow(rds.data.nomiss))
  } else{
    subset[is.na(subset)] <- FALSE
  }
  # Sub-population size for reporting and the SRS variance; the weights and
  # the bootstrap engines are fit on the full data with the full N.
  N.report <- N
  if (!all(subset) && !is.null(N)) {
    #use VH estimator to adjust population size to sub-population
    # defined by subset
    tmp.wts <- vh.weights(rds.data.nomiss[[network.size]])
    tmp.wts <- tmp.wts / sum(tmp.wts)
    prop <- sum(tmp.wts * subset)
    N.report <- N * prop
    if (N.report < sum(subset))
      stop(sprintf("Estimated sub-population size, %f, smaller than subsets size, %f.", N.report, sum(subset)))
  }

  # Boundary: retain excluded subjects (missing outcome / outside subset) as
  # an internal group so recruitment chains stay intact; the group is removed
  # from all outputs. A continuous outcome cannot carry the internal level in
  # the outcome vector itself: excluded subjects' values are set to NA, and
  # the engines carry them as the internal level of the *cut* factor - each
  # (group, degree) cell's members stand in for one another, so an excluded
  # subject contributes its degree and linkage but draws an NA value, which
  # the estimates mask (see NA_SUBSET_HANDLING.md).
  outcome0 <- rds.data.nomiss[[outcome.variable]]
  is.cts <- is.numeric(outcome0) & !to.factor
  excl <- .excluded.mask(outcome0, subset)
  if (!is.cts) {
    if (any(excl)) {
      rds.data.nomiss[[outcome.variable]] <-
        .sort.levels(.mint.excluded(as.factor(outcome0), excl))
    }
  } else if (any(excl & !is.na(outcome0))) {
    # Subjects outside the subset must not contribute outcome information.
    rds.data.nomiss[[outcome.variable]][excl] <- NA
  }

  #if (mean(duplicated(rds.data.nomiss[[outcome.variable]], na.rm = TRUE)) < control$discrete.cutoff) {
  #  is.cts <- TRUE
  #}
  outcome <- rds.data.nomiss[[outcome.variable]]
  #only categorical estimates are supported
  if(!is.cts){
    outcome <- factor(rds.data.nomiss[[outcome.variable]])
           #Make sure the factor labels are alphabetic!
    outcome <- .sort.levels(outcome)
    rds.data.nomiss[[outcome.variable]] <- outcome
    outclasses <- setdiff(levels(outcome), .rds.excluded)
    if(any(excl)){
      nsamplesbyoutcome <- table(factor(as.character(outcome[!excl]), levels = outclasses))
    }else{
      nsamplesbyoutcome <- table(as.numeric(outcome))
    }
  }else{
    nsamplesbyoutcome <- sum(!excl & !is.na(outcome))
  }
  g <- length(setdiff(stats::na.omit(unique(as.character(rds.data.nomiss[[outcome.variable]]))), .rds.excluded))
  if (g == 1) {
    warning(paste(outcome.variable, "has only one level. Skipping..."))
    return(invisible(NULL))
  }
  nsamples <- sum(!excl & !is.na(as.vector(outcome)))

  rds <- RDS.estimates.local(
    rds.data = rds.data.nomiss,
    outcome.variable = outcome.variable,
    N = N,
    weight.type = weight.type,
    control = control,
    to.factor=to.factor,
    cont.breaks=cont.breaks,
    hcg.wave.fallback = TRUE,
    ...
  )
  
  if (is.rds.interval.estimate(rds)) {
    observed.estimate <- rds$estimate
    weights.all <- rds$weights
  } else{
    observed.estimate <- rds@estimate
    weights.all <- rds@weights
  }
  
  # outcome variance 
  if(!is.cts){
    varoutcome <- observed.estimate * (1 - observed.estimate) 
  }else{
    varoutcome <- wtd.var(rds.data.nomiss[[outcome.variable]], weights = weights.all)
  }
  
  if (is.null(number.of.bootstrap.samples) ||
      number.of.bootstrap.samples < 10) {
    if (is.rds.interval.estimate(rds)) {
      EL.se <- matrix(rds$interval, ncol = 6, byrow = FALSE)[1, 5]
    } else{
      EL.se <-
        sqrt(varoutcome / length(weights.all))
    }
    number.of.bootstrap.samples <-
      min(500, max(20, ceiling(0.5 * (EL.se / 0.002) ^ 2 + 1)))
    if (!is.numeric(number.of.bootstrap.samples)) {
      number.of.bootstrap.samples <- 500
    }
    cat(sprintf(
      "Using %d bootstrap samples.\n",
      number.of.bootstrap.samples
    ))
  }
  #
  #  Confidence intervals
  #
  crit <- stats::qnorm((1 - confidence.level) / 2, lower.tail = FALSE)
  if (uncertainty == "Salganik") {
    bs <- salganik.bootstrap.se(
      rds.data = rds.data.nomiss,
      group.variable = outcome.variable,
      number.of.bootstrap.samples = number.of.bootstrap.samples,
      estimator.name = weight.type,
      N = N,
      to.factor=to.factor,
      ...
    )
    attr(bs, "is.cts") <- is.cts 
    if (!is.cts) {
      ord    <- match(outclasses, names(bs))        # <- FIX: correct direction (was names(bs), outclasses)
      bso    <- bs[ord]
      bsests <- attr(bs, "bsresult")[, ord, drop = FALSE]
    } else {
      bso    <- bs                                  # <- cts: single column, NO reorder (fixes the crash)
      bsests <- as.matrix(attr(bs, "bsresult"))
      outclasses <- outcome.variable                # normalize the cts label here for all branches
    }
    estimate <- cbind(observed.estimate, observed.estimate - crit * bso, observed.estimate + crit * bso)
    colnames(estimate) <- c("point", "lower", "upper")
    estimate <- .bootstrap.ci(estimate, observed.estimate, bsests, ci.type, confidence.level)
    
  } else if (uncertainty == "HCG") {

    bs <- HCG.bootstrap.se(
      rds.data = rds.data.nomiss,
      group.variable = outcome.variable,
      number.of.bootstrap.samples = number.of.bootstrap.samples,
      estimator.name = weight.type,
      N = N,
      cont.breaks=cont.breaks,
      to.factor=to.factor,
      control = control,
      ...
    )
    attr(bs, "is.cts") <- is.cts 
    if (!is.cts) {
      ord    <- match(outclasses, names(bs))        # <- FIX: correct direction (was names(bs), outclasses)
      bso    <- bs[ord]
      bsests <- attr(bs, "bsresult")[, ord, drop = FALSE]
    } else {
      bso    <- bs                                  # <- cts: single column, NO reorder (fixes the crash)
      bsests <- as.matrix(attr(bs, "bsresult"))
      outclasses <- outcome.variable                # normalize the cts label here for all branches
    }
    estimate <- cbind(observed.estimate, observed.estimate - crit * bso, observed.estimate + crit * bso)
    colnames(estimate) <- c("point", "lower", "upper")
    estimate <- .bootstrap.ci(estimate, observed.estimate, bsests, ci.type, confidence.level)
    
  } else if (uncertainty == "Treeboot") {
    if(is.cts)
      stop("uncertainty = \"Treeboot\" does not support continuous outcomes; use \"Salganik\", \"HCG\", or \"Gile\".")
    bs <- treeboot.bootstrap.se(
      rds.data = rds.data.nomiss,
      group.variable = outcome.variable,
      number.of.bootstrap.samples = number.of.bootstrap.samples,
      estimator.name = weight.type,
      N = N,
      ...
    )
    attr(bs, "is.cts") <- is.cts 
    if (!is.cts) {
      ord    <- match(outclasses, names(bs))        # <- FIX: correct direction (was names(bs), outclasses)
      bso    <- bs[ord]
      bsests <- attr(bs, "bsresult")[, ord, drop = FALSE]
    } else {
      bso    <- bs                                  # <- cts: single column, NO reorder (fixes the crash)
      bsests <- as.matrix(attr(bs, "bsresult"))
      outclasses <- outcome.variable                # normalize the cts label here for all branches
    }
    estimate <- cbind(observed.estimate, observed.estimate - crit * bso, observed.estimate + crit * bso)
    colnames(estimate) <- c("point", "lower", "upper")
    estimate <- .bootstrap.ci(estimate, observed.estimate, bsests, ci.type, confidence.level)
  } else if (uncertainty == "SRS") {
    if(is.cts){
      varsrs <- varoutcome / nsamples
    }else{
      varsrs <- observed.estimate * (1 - observed.estimate) / nsamples
    }
    if (!is.null(N))
      varsrs <- ((N.report - nsamples) / (N.report - 1)) * varsrs

    estimate <- cbind(
      observed.estimate,
      observed.estimate -
        crit * sqrt(varsrs),
      observed.estimate + crit * sqrt(varsrs)
    )
    colnames(estimate) <- c("point", "lower", "upper")
    bs <- varsrs
    attr(bs, "is.cts") <- is.cts
    #attr(bs, "is.quantile") <- is.quantile
    attr(bs, "mu") <- observed.estimate
    
  } else if (uncertainty == "Gile") {
    bs <- SSBS.estimates(
      rds.data.nomiss,
      outcome.variable,
      confidence.level = confidence.level,
      N = N,
      refs.to.trait=NULL,
      weight.type = weight.type,
      number.of.bootstrap.samples = number.of.bootstrap.samples,
      control = control,
      to.factor=to.factor,
      cont.breaks=cont.breaks,
      fast = fast,
      useC = useC,
      ...
    )
    if (attr(bs, "is.cts")) {
      observed.estimate <- bs[1]; outclasses <- outcome.variable
      estimate <- matrix(bs, nrow = 1); bsests <- as.matrix(attr(bs, "bsresult")$bsests)
    } else {
      # If the internal excluded group was simulated, remove it (by name)
      # from each replicate and renormalize before taking SDs/quantiles.
      bsr <- attr(bs, "bsresult")
      if (.rds.excluded %in% colnames(bsr$bsests)) {
        keep   <- colnames(bsr$bsests) != .rds.excluded
        be     <- bsr$bsests[, keep, drop = FALSE]
        be     <- be / rowSums(be)
        pt     <- bsr$point_estimate[names(bsr$point_estimate) != .rds.excluded]
        pt     <- pt / sum(pt)
        sds    <- apply(be, 2, stats::sd, na.rm = TRUE)
        bs2    <- cbind(pt, pt - crit * sds, pt + crit * sds)
        rownames(bs2) <- names(pt)
        bsr$bsests <- be
        bsr$point_estimate <- pt
        bsr$se_estimate <- sds
        attributes(bs2) <- c(attributes(bs2),
                             attributes(bs)[setdiff(names(attributes(bs)),
                                                    c("dim","dimnames","bsresult"))])
        attr(bs2, "bsresult") <- bsr
        bs <- bs2
      }
      ord      <- match(outclasses, rownames(bs))   # <- FIX: correct direction
      estimate <- bs[ord, , drop = FALSE]
      bsests   <- attr(bs, "bsresult")$bsests[, ord, drop = FALSE]
    }
    colnames(estimate) <- c("point", "lower", "upper")
    estimate <- .bootstrap.ci(estimate, observed.estimate, bsests, ci.type, confidence.level)
  }
  #
  # Design effect (for factor outcomes)
  #
  #if (attr(bs, "is.cts")) {
    #if (attr(bs, "is.quantile")) {
    #  varoutcome <-
    #    continuous * (1 - continuous) / (stats::dnorm(
    #      observed.estimate,
    #      mean = attr(bs, "mu"),
    #      sd = sqrt(attr(bs, "sigma2"))
    #    )) ^ 2
    #} else{
    #varoutcome <- attr(bs, "sigma2")
    #}
  #} else{
    #varoutcome <- estimate[, 1] * (1 - estimate[, 1])
  #}
  #       Note the finite sample correction factor
  if (!is.null(N)) {
    varsrs <- (((N.report - nsamples) / (N.report - 1)) * varoutcome / nsamples)
  } else{
    varsrs <- varoutcome / nsamples
  }
  de <- ((estimate[, 3] - estimate[, 2]) / (crit * 2)) ^ 2 / varsrs
  estimate <-
    cbind(estimate,
          de,
          (estimate[, 3] - estimate[, 2]) / (crit * 2),
          nsamplesbyoutcome)
  colnames(estimate)[c(4, 5, 6)] <- c("Design Effect", "s.e.", "n")
  if (attr(bs, "is.cts")) {
    rownames(estimate) <- outcome.variable
  } else{
    outclasses <-
      setdiff(sort(unique(as.vector(rds.data.nomiss[[outcome.variable]]))), .rds.excluded)
    outclasses[outclasses == "NA.NA"] <- "NA"
    rownames(estimate) <- outclasses[1:nrow(estimate)]
    estimate <- as.numeric(estimate)
    names(estimate)[1:g] <- outclasses[1:g]
  }
  if (exists("bs")){
    attr(estimate, "bsresult") <- attr(bs, "bsresult")
    attr(estimate, "is.cts") <- attr(bs, "is.cts")
  }
  
  result <- rds.interval.estimate(
    estimate,
    outcome.variable,
    weight.type,
    uncertainty,
    weights.all,
    N = N.report,
    csubset = csubset,
    conf.level = confidence.level
  )
  attr(result, "bsresult") <- attr(estimate, "bsresult")
  attr(result, "is.cts") <- attr(estimate, "is.cts")
  
  return(result)
}






#' RDS Bootstrap Interval Estimates
#'
#' This function computes an interval estimate for one or more categorical
#' variables. It optionally uses attributes of the RDS data set to determine
#' the type of estimator and type of uncertainty estimate to use.
#'
# Note that the current estimator only works on binary outcome variables. (I think this cane be removed now--non-binary categorical outcome variables work.)
#'
#' @param rds.data An \code{rds.data.frame} that indicates recruitment patterns
#' by a pair of attributes named ``id'' and ``recruiter.id''.
#' @param outcome.variable A string giving the name of the variable in the
#' \code{rds.data} that contains a categorical or numeric variable to be
#' analyzed.
#' @param weight.type A string giving the type of estimator to use. The options
#' are \code{"Gile's SS"}, \code{"RDS-I"}, \code{"RDS-II"}, \code{"RDS-I (DS)"},
#' and \code{"Arithemic Mean"}. If \code{NULL} it defaults to \code{"Gile's
#' SS"}.
#' @param uncertainty A string giving the type of uncertainty estimator to use.
#' The options are \code{"SRS"}, \code{"Salganik"}, \code{"Treeboot"},
#' \code{"Gile"} and \code{"HCG"}; see \code{\link{RDS.bootstrap.mechanics}}
#' for how each resampling procedure works and what it assumes. This is usually
#' determined by \code{weight.type} to be consistent with the estimator's
#' origins. The estimators RDS-I, RDS-I (DS), and RDS-II default to \code{"Salganik"},  "Arithmetic
#' Mean" defaults to \code{"SRS"}, "Gile's SS" defaults to the \code{"Gile"}
#' bootstrap and "HCG" defaults to the \code{"HCG"} bootstrap.
#' @param N An estimate of the number of members of the population being
#' sampled. If \code{NULL} it is read as the \code{population.size.mid} attribute of
#' the \code{rds.data} frame. If that is missing it defaults to 1000.
#' @param subset An optional criterion to subset \code{rds.data} by. It is a
#' character string giving an R expression which, when evaluated, subset the
#' data. In plain English, it can be something like \code{"seed > 0"} to
#' exclude seeds. It can be the name of a logical vector of the same length of
#' the outcome variable where TRUE means include it in the analysis. If
#' \code{NULL} then no subsetting is done.
#' @param confidence.level The confidence level for the confidence intervals.
#' The default is 0.95 for 95\%.
#' @param number.of.bootstrap.samples The number of bootstrap samples to take
#' in estimating the uncertainty of the estimator. If \code{NULL} it defaults
#' to the number necessary to compute the standard error to accuracy 0.001.
#' \code{outcome.variable}. Otherwise it will compute the population frequencies of each value of the \code{outcome.variable}.
#' @param fast Use a fast bootstrap where the weights are reused from the
#' estimator rather than being recomputed for each bootstrap sample.
#' @param useC Use a C-level implementation of Gile's bootstrap (rather than
#' the R level). The implementations should be a computational
#' equivalent estimator (except for speed).
#' @param ci.type Type of confidence interval to use, if possible. If "t", use lower and upper confidence
#' interval values based on the standard deviation of the bootstrapped values and a t multiplier.
#' If "pivotal", use lower and upper confidence interval values based on the
#' basic bootstrap (also called the pivotal confidence interval).
#' If "quantile", use lower and upper confidence interval values based on the
#' quantiles of the bootstrap sample.
#' If "proportion", use the "t" unless the estimated proportion is less than 0.15 or the bounds are
#' outside [0,1 . In this case, try the "quantile" and constrain the bounds to be compatible with [0,1].
#' @param control A list of control parameters for algorithm
#' tuning. Constructed using\cr
#' \code{\link{control.rds.estimates}}.
#' @param to.factor force variable to be a factor
#' @param cont.breaks For continuous variates, some bootstrap proceedures require categorical data.
#' In these cases, in order to contruct each bootstrap replicate, the outcome variable is split into
#' cont.breaks categories.
#' @param ... Additional arguments for RDS.*.estimates.
#' @return An object of class \code{rds.interval.estimate} summarizing the inference.
#' The confidence interval and standard error are based on the bootstrap procedure.
#' In additon, the object has attribute \code{bsresult} which provides details of the
#' bootstrap procedure. The contents of the \code{bsresult} attribute depends on the
#' \code{uncertainty} used. If \code{uncertainty=="Salganik"} then \code{bsresult} is a
#' vector of standard deviations of the bootstrap samples.
#' If \code{uncertainty=="Gile's SS"} then
#' \code{bsresult} is a list with components for the bootstrap point estimate,
#' the bootstrap
#' samples themselves and the standard deviations of the bootstrap samples.
#' If \code{uncertainty=="SRS"} then \code{bsresult} is NULL.
#' @references Gile, Krista J. 2011 \emph{Improved Inference for
#' Respondent-Driven Sampling Data with Application to HIV Prevalence
#' Estimation}, \emph{Journal of the American Statistical Association}, 106,
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
#' @seealso \code{\link{RDS.bootstrap.mechanics}} for a description of each
#' bootstrap procedure, and \code{\link{control.rds.estimates}} for tuning
#' parameters.
#' @keywords survey manip
#' @examples
#'
#' \dontrun{
#' data(fauxmadrona)
#' RDS.bootstrap.intervals(rds.data=fauxmadrona,weight.type="RDS-II",
#'      uncertainty="Salganik",
#' 	outcome.variable="disease",N=1000,number.of.bootstrap.samples=50)
#'
#' data(fauxtime)
#' RDS.bootstrap.intervals(rds.data=fauxtime,weight.type="HCG",
#'      uncertainty="HCG",
#' 	outcome.variable="var1",N=1000,number.of.bootstrap.samples=10)
#' }
#'
#' @export
RDS.bootstrap.intervals <- function(rds.data,
                                    outcome.variable,
                                    weight.type = NULL,
                                    uncertainty = NULL,
                                    N = NULL,
                                    subset = NULL,
                                    confidence.level = .95,
                                    number.of.bootstrap.samples = NULL,
                                    fast = TRUE,
                                    useC = TRUE,
                                    ci.type = "t",
                                    control = control.rds.estimates(),
                                    to.factor=FALSE, 
                                    cont.breaks = 3,
                                    ...) {
  se <- substitute(subset)
  if (is.null(se)) {
    csubset <- ""
  } else{
    csubset <- as.character(enquote(substitute(subset)))[2]
  }
  subset <- eval(se, rds.data, parent.frame())
  if (length(outcome.variable) == 1) {
    result <- RDS.bootstrap.intervals.local(
      rds.data,
      outcome.variable,
      weight.type,
      uncertainty,
      N,
      subset,
      confidence.level,
      number.of.bootstrap.samples = number.of.bootstrap.samples,
      control = control,
      fast = fast,
      useC = useC,
      ci.type = ci.type,
      csubset = csubset,
      to.factor = to.factor,
      cont.breaks = cont.breaks,
      ...
    )
  }
  else {
    result <- lapply(
      X = outcome.variable,
      FUN = function(g) {
        RDS.bootstrap.intervals.local(
          rds.data,
          g,
          weight.type,
          uncertainty,
          N,
          subset,
          confidence.level,
          number.of.bootstrap.samples = number.of.bootstrap.samples,
          control = control,
          fast = fast,
          useC = useC,
          ci.type = ci.type,
          csubset = csubset,
          to.factor = to.factor,
          cont.breaks = cont.breaks,
          ...
        )
      }
    )
    names(result) <- outcome.variable
  }
  return(result)
}
