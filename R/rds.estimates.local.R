RDS.estimates.local <- function(rds.data,
           outcome.variable,
           subset = NULL,
           control = control.rds.estimates(),
           weight.type,
           N = NULL,
           to.factor = FALSE,
           cont.breaks = 3,
           hcg.wave.fallback = FALSE,
           ...) {
    if (is(rds.data, "rds.data.frame")) {
      if (!(outcome.variable %in% names(rds.data))) {
        stop(sprintf(
          "No variable called %s appears in the data.",
          outcome.variable
        ))
      }
      network.size <- attr(rds.data, "network.size.variable")
      
    }
    else{
      stop("rds.data must be of type rds.data.frame")
    }
    
    if (to.factor) {
      outcome <- factor(rds.data[[outcome.variable]])
      # 		Make sure the factor labels are alphabetic!
      outcome <- factor(outcome, levels = levels(outcome)[order(levels(outcome))])
      rds.data[[outcome.variable]] <- outcome
    }
    
    
    if (is.null(N)) {
      N <- get.population.size(rds.data)[2]
      if(is.na(N))
        N <- NULL
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
    
    #########################################################################################
    # The simple data management tasks have been taken care of, so it's now time to compute
    # the RDS estimates.   The cases for numeric and categorical outcomes are handled
    # separately.
    
    se <- substitute(subset)
    subset <- eval(se, rds.data, parent.frame())
    if (is.null(se) | is.null(subset)) {
      csubset <- ""
      subset <- rep(TRUE, length = nrow(rds.data.nomiss))
    } else{
      csubset <- as.character(enquote(substitute(subset)))[2]
      if (nchar(csubset) > 20)
        csubset <- "subsetted"
      subset[is.na(subset)] <- FALSE
    }

    # Boundary: subjects excluded by the subset expression (and, for
    # categorical outcomes, missing outcomes) are retained internally as
    # their own group so that recruitment chains through them stay intact.
    # They receive weight 0 and are removed from the estimate below.
    outcome0 <- rds.data.nomiss[[outcome.variable]]
    is.cts.outcome <- is.numeric(outcome0)
    excl <- .excluded.mask(outcome0, subset)
    if (!is.cts.outcome && any(excl)) {
      rds.data.nomiss[[outcome.variable]] <-
        .sort.levels(.mint.excluded(as.factor(outcome0), excl))
    }

    force.discrete <- is.cts.outcome && weight.type %in% c("HCG","RDS-I","RDS-I (DS)")
    if(force.discrete){
      #out <- cut(rds.data.nomiss[[outcome.variable]],
      #               breaks = unique(quantile(rds.data.nomiss[[outcome.variable]],
      #                                        probs = seq(0, 1, length.out=cont.breaks+1))), include.lowest=TRUE)
      out <- suppressWarnings(.cut2.excl(outcome0, g=cont.breaks, mask=excl))
      tmp.out <- outcome0
      rds.data.nomiss[[outcome.variable]] <- out
    }
    weights.all <- compute.weights(
      rds.data.nomiss,
      weight.type = weight.type,
      outcome.variable = outcome.variable,
      N = N,
      control = control,
      hcg.wave.fallback = hcg.wave.fallback,
      ...
    )
    if(force.discrete){
      rds.data.nomiss[[outcome.variable]] <- tmp.out
    }
    weights.all[excl] <- 0

    #########################################################################################
    # The simple data management tasks have been taken care of, so it's now time to compute
    # the RDS-II estimates.   The cases for numeric and categorical outcomes are handled
    # separately.
    outcome <- rds.data.nomiss[[outcome.variable]]

    if (is.numeric(outcome)) {
      estimate <- HT.estimate(weights = weights.all[!excl], outcome = outcome[!excl])
    } else{
      outcome <- as.factor(outcome)
      # Restrict to retained subjects and observed levels (the internal
      # excluded level is removed by name, never by position).
      o <- outcome[!excl]
      keep.lev <- setdiff(levels(o), .rds.excluded)
      # 		Make sure the factor labels are alphabetic!
      o <- factor(o, levels = keep.lev[order(keep.lev)])

      # Calculate the within group totals of the *inverse* of the weights, again we exclude the seeds.
      wtd.degree.totals <- stats::xtabs(weights.all[!excl] ~ o)

      estimate.xtabs <- wtd.degree.totals / sum(wtd.degree.totals)
      # Now convert to a numeric vector ('estimate.xtabs' is currently of class 'xtabs')
      estimate <- as.numeric(estimate.xtabs)
      names(estimate) <- names(estimate.xtabs)
    }
    
    result <-
      rds.weighted.estimate(
        estimate,
        weights.all,
        outcome.variable,
        weight.type,
        subset = subset,
        csubset = csubset
      )    
    return(result)
  }
