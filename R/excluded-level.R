################################################################################
# Internal machinery for representing excluded subjects (missing outcome or
# outside a subset= expression) as an explicit factor level, so that the
# recruitment tree stays intact through the weighting and bootstrap machinery
# and the excluded group is removed only at the output boundary.
#
# Boundary functions (compute.weights, RDS.estimates.local,
# RDS.bootstrap.intervals and the secondary boundaries) own all recoding and
# stripping. Estimation engines below the boundary receive complete data in
# which the excluded group is an ordinary factor level; they must never see
# missing outcome values (see .assert.complete.outcome).
#
# The reserved level is always located by name, never by position.
# See NA_SUBSET_HANDLING.md for the full design.
################################################################################

.rds.excluded <- "~~~excluded~~~"

# TRUE where a subject is excluded from analysis: missing outcome, already
# carrying the reserved level (data arriving from an internal caller such as
# a bootstrap replicate), or outside the subset.
.excluded.mask <- function(outcome, subset = NULL){
  m <- is.na(outcome)
  m <- m | (as.character(outcome) %in% .rds.excluded)
  if(!is.null(subset))
    m <- m | !subset
  m
}

# Recode a categorical outcome so that excluded subjects form their own
# level, appended after the observed levels. Returns the input unchanged
# when nothing is excluded, so NA-free full-data analyses are untouched.
# Empty observed levels are dropped (matching the historical "drop 0 count
# levels" behavior under subsetting); existing level order is preserved.
.mint.excluded <- function(outcome, mask){
  if(!any(mask))
    return(outcome)
  x <- as.character(outcome)
  x[mask] <- .rds.excluded
  if(is.factor(outcome)){
    lev <- levels(outcome)
    lev <- lev[lev %in% x[!mask] & lev != .rds.excluded]
  }else{
    lev <- sort(unique(x[!mask]))
  }
  factor(x, levels = c(lev, .rds.excluded))
}

# Alphabetize observed levels; the reserved level, if present, is placed
# last for display stability. No correctness logic may depend on its
# position -- strip operations match by name.
.sort.levels <- function(f){
  lev <- levels(f)
  obs <- sort(lev[lev != .rds.excluded])
  if(.rds.excluded %in% lev)
    factor(f, levels = c(obs, .rds.excluded))
  else
    factor(f, levels = obs)
}

# Cut a continuous outcome into g groups using only retained values (so
# excluded subjects' values never shape the breakpoints); excluded subjects
# get the reserved level.
.cut2.excl <- function(x, g, mask){
  if(!any(mask))
    return(.cut2(x, g))
  f <- rep(NA_character_, length(x))
  cf <- .cut2(x[!mask], g = g)
  f[!mask] <- as.character(cf)
  f[mask] <- .rds.excluded
  factor(f, levels = c(levels(cf), .rds.excluded))
}

# Drop the excluded entry (matched by name) from a named estimate vector,
# renormalizing when it is a proportion vector.
.drop.excluded <- function(est, renormalize = TRUE){
  if(is.null(names(est)))
    return(est)
  keep <- names(est) != .rds.excluded
  if(all(keep))
    return(est)
  e <- est[keep]
  if(renormalize && isTRUE(sum(e) > 0))
    e <- e / sum(e)
  e
}

# Engine-side guard: estimation engines are only ever handed complete
# outcomes (the boundary recodes missing values into the reserved level).
.assert.complete.outcome <- function(outcome){
  if(anyNA(outcome))
    stop("missing outcome values reached an RDS estimation engine. ",
         "Missing-data and subset handling is performed by compute.weights, ",
         "the RDS.*.estimates functions and RDS.bootstrap.intervals; ",
         "please use those entry points.", call. = FALSE)
  invisible(TRUE)
}

# Impute missing degrees with the within-group median of the observed
# degrees (global median where a group has none). Shared by the HCG
# weighting and replicate-weight code paths.
.impute.degree.group.median <- function(degree, out){
  mdeg <- is.na(degree)
  if(!any(mdeg))
    return(degree)
  warning(paste(sum(mdeg), "subject(s) missing degree. Imputing within group median."),
          call. = FALSE)
  med <- stats::median(degree, na.rm = TRUE)
  levs <- unique(out[!is.na(out)])
  for(l in levs){
    grp <- !is.na(out) & out == l
    if(any(grp & mdeg)){
      gm <- stats::median(degree[grp], na.rm = TRUE)
      if(is.na(gm))
        gm <- med
      degree[grp & mdeg] <- gm
    }
  }
  degree[is.na(degree)] <- med
  degree
}
