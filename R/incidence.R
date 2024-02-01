
# Calculates assay based incidence estimate
# 
# @param recent logical indicator of recent infection
# @param hiv logical indicator of hiv infection
# @param weights sampling weights
# @param mean.duration Estimated mean duration of recent infection (MDRI) (days)
# @param frr Estimated false-recent rate (FRR)
# @param post.infection.cutoff Post-infection time cut-off T, separating "true-recent" from "false-recent" results (days)
# @examples 
# hiv <- rep(FALSE, 5000)
# hiv[1:1000] <- TRUE
# recent <- rep(NA, 5000)
# recent[1:1000] <- FALSE
# recent[1:50] <- TRUE
# incidence(recent,hiv)
incidence <- function(recent, hiv, weights=NULL, mean.duration=200, frr=.01, post.infection.cutoff=730){
  mdrihat <- mean.duration / 365.25
  post.infection.cutoff <- post.infection.cutoff / 365.25
  p_s <- wtd.mean(!hiv, weights=weights, na.rm=TRUE)
  p_pos <- 1 - p_s
  p_RgivenTested <- wtd.mean(recent, weights=weights, na.rm=TRUE)
  p_NRgivenTested <- 1 - p_RgivenTested
  #incest_den <- p_s*(mdrihat-T*frrhat)
  max(0, (p_RgivenTested*p_pos-frr*p_pos)/(p_s*(mdrihat-frr*post.infection.cutoff)))
}


#' Calculates incidence and bootstrap confidence intervals for immunoassay data collected with RDS
#' @param rds.data an rds.data.frame
#' @param recent.variable The name of the variable indicating recent infection 
#' @param hiv.variable The name of the variable indicating of hiv infection
#' @param N Population size
#' @param weight.type A string giving the type of estimator to use. The options
#' are \code{"Gile's SS"}, \code{"RDS-I"}, \code{"RDS-II"}, \code{"RDS-I/DS"},
#' and \code{"Arithemic Mean"}. It defaults to \code{"Gile's
#' SS"}.
#' @param mean.duration Estimated mean duration of recent infection (MDRI) (days)
#' @param frr Estimated false-recent rate (FRR)
#' @param post.infection.cutoff Post-infection time cut-off T, separating "true-recent" from "false-recent" results (days)
#' @param number.of.bootstrap.samples The number of bootstrap samples used to construct the interval.
#' @param se.mean.duration The standard error of the mean.duration estimate
#' @param se.frr The standard error of the false recency estimate
#' @param confidence.level The level of confidence for the interval
#' @param verbose verbosity control
#' @param ... additional arguments to compute.weights
#' @details 
#' The recent.variable and hiv should be the names of logical variables. Otherwise they are converted to logical using as.numeric(x) > 0.5.
#' 
#' This function estimates incidence using RDS sampling wieghts. Confidence intervals are constucted using HCG bootstraps.
#' See http://www.incidence-estimation.org/ for additional information on (non-RDS) incidence estimation.
#' @examples 
#' data(faux)
#' faux$hiv <- faux$X == "blue"
#' faux$recent <- NA
#' faux$recent[faux$hiv] <- runif(sum(faux$hiv)) < .2
#' faux$recent[runif(nrow(faux)) > .5] <- NA
#' faux$hiv[is.na(faux$recent)][c(1,6,10,21)] <- NA
#' attr(faux,"time") <- "wave"
#' bootstrap.incidence(faux,"recent","hiv",weight.type="RDS-II", number.of.bootstrap.samples=100)
#' @export
bootstrap.incidence <- function(rds.data, recent.variable, hiv.variable,  N=NULL, weight.type=c("Gile's SS","RDS-I","RDS-I (DS)","RDS-II","Arithmetic Mean","HCG"),
                                    mean.duration=200, frr=.01, post.infection.cutoff=730, 
                                    number.of.bootstrap.samples = 1000,
                                    se.mean.duration=0, se.frr=0,confidence.level = .95,
                                verbose=TRUE, ...){
  
  weight.type <- match.arg(weight.type,
                           c("Gile's SS","RDS-I", "RDS-II", "RDS-I (DS)","Arithmetic Mean","HCG"))
  
  network.size <- attr(rds.data, "network.size.variable")
  remvalues <- rds.data[[network.size]] == 0 | is.na(rds.data[[network.size]])
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
    rds.data[[network.size]][remvalues] <- max(rds.data[[network.size]], na.rm = TRUE)
  }
  
  #Set up
  n <- nrow(rds.data)
  if(is.null(N))
    N <- n * 100000
  wave <- get.wave(rds.data)
  if (!has.recruitment.time(rds.data)) {
    time <- rep(1, nrow(rds.data))
  } else{
    time <- get.recruitment.time(rds.data)
  }
  rds <- rds.data[order(time, wave),]
  if (!has.recruitment.time(rds)) {
    time <- rep(1, nrow(rds.data))
  } else{
    time <- get.recruitment.time(rds)
  }
  
  # Create combination variable of both row and column
  hiv <- as.numeric(rds.data[[hiv.variable]]) > .5
  recent <- as.numeric(rds.data[[recent.variable]]) > .5
  v <- rep(NA, n)
  v[!hiv] <- 0
  v[hiv & is.na(recent)] <- 1
  v[hiv & !recent] <- 2
  v[hiv & recent] <- 3
  v <- as.factor(v)
  varname <- .make.unique(names(rds), "variable")
  rds[[varname]] <- v
  
  inc.func <- function(x, boot.frr=TRUE, boot.mean.duration=TRUE){
    wts <- compute.weights(
      x,
      outcome.variable = varname,
      weight.type = weight.type,
      N = N,
      ...
    )
    wts <- wts * nrow(x) / sum(wts)
    
    #bootstrap these estimates
    if(boot.frr)
      frr.boot <- max(0, stats::rnorm(1, frr, se.frr))
    else
      frr.boot <- frr
    if(boot.mean.duration)
      mean.duration.boot <- max(1, stats::rnorm(1, mean.duration, se.mean.duration))
    else
      mean.duration.boot <- mean.duration
    
    v <- x[[varname]]
    hiv <- v %in% c("1","2","3")
    hiv[is.na(v)] <- NA
    recent <- v == "3"
    recent[v %in% c("0","1")] <- NA
    incidence(recent, hiv, weights=wts, mean.duration=mean.duration.boot, 
              frr=frr.boot, post.infection.cutoff=post.infection.cutoff)
  }
  
  estimate <- inc.func(rds, FALSE, FALSE)
  
  # Perform bootstrap
  boot.stats <-
    unlist(
      HCG.bootstrap(
        rds,
        varname,
        number.of.bootstrap.samples,
        N = N,
        fun = inc.func,
        small.fraction = TRUE,
        verbose = verbose
      )
    )
  boot.var <- stats::var(boot.stats)
  
  mult <- -stats::qnorm((1 - confidence.level) / 2)
  result <- data.frame(Incidence=estimate, SE=sqrt(boot.var),
                       Lower=max(0, estimate - mult*sqrt(boot.var)),
                       Upper=estimate + mult*sqrt(boot.var))
  rownames(result) <- recent.variable
  result
}
