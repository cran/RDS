#' Bootstrap uncertainty procedures for respondent-driven sampling
#'
#' @description An overview of the resampling procedures that
#' \code{\link{RDS.bootstrap.intervals}} uses to attach standard errors and
#' confidence intervals to respondent-driven sampling (RDS) estimates, and of
#' what each one assumes. The procedure is selected with the \code{uncertainty}
#' argument.
#'
#' @details
#' \strong{Why a special bootstrap is needed.} In RDS, participants recruit
#' people they already know, so the sample is not a simple random sample. It is
#' a branching recruitment tree; participants are reached with a probability
#' related to their number of social contacts (their degree); neighbours tend
#' to share traits (homophily); some groups recruit more actively than others
#' (differential activity); and when a study captures a large share of the
#' hidden population, later recruits deplete the pool of people still available
#' to be recruited. Any of these features can make an ordinary standard error
#' badly wrong -- usually too small. Each procedure below quantifies uncertainty
#' by generating many artificial datasets that reproduce a particular model of
#' the RDS process, recomputing the estimate on each dataset, and using the
#' spread of those estimates as the standard error (and, for the chosen interval
#' type, the confidence limits). The procedures differ in which features of the
#' design they reproduce. The default is chosen to match the estimator in use;
#' it can be overridden with the \code{uncertainty} argument.
#'
#' \strong{SRS (simple random sampling).} The naive baseline. It treats the
#' sample as if it were an ordinary simple random sample and reports the
#' design-weighted standard deviation of the outcome divided by the square root
#' of the sample size, with a finite-population correction when the population
#' size is known. It ignores recruitment dependence and degree bias, so it
#' generally understates the true uncertainty. It is useful mainly as a
#' reference or lower bound, and is the default for the arithmetic-mean
#' estimator.
#'
#' Steps (a direct calculation rather than a resampling procedure):
#'
#' 1. Compute the design-weighted variance of the outcome.
#'
#' 2. Divide by the sample size, and apply a finite-population correction when
#' the population size is known, to obtain the squared standard error.
#'
#' 3. Form the interval from the point estimate and this standard error.
#'
#' \strong{Salganik bootstrap.} Reproduces recruitment as a first-order Markov
#' chain on the trait of interest. From the observed recruiter-to-recruit pairs
#' it estimates the chance that a participant in one trait group recruits a
#' participant in each other group. Each artificial sample is built by starting
#' at a seed and repeatedly drawing the next recruit's group from those
#' transition probabilities until the original sample size is reached; within a
#' group, an individual's degree is drawn so as to reproduce the degree profile
#' of recruits in that group. This captures the trait-to-trait dependence and
#' the degree structure, but it assumes the recruitment chain behaves like a
#' stationary Markov chain and that the sample is a small fraction of the
#' population (so that no depletion occurs). It is the standard companion to the
#' RDS-I and RDS-II (Volz-Heckathorn) estimators.
#' 
#' Bootstrap steps (to generate one bootstrap sample):
#'
#' 0. Count the recruitment ties between each pair of groups to estimate the
#' probability \eqn{\theta_{ij}}{theta_ij} that a participant in group
#' \eqn{i}{i} recruits a participant in group \eqn{j}{j}.
#'
#' 1. Select a starting participant (seed) at random from the sample as the
#' active recruiter.
#'
#' 2. Draw a group for the recruit from the rate that people in the active
#' recruiter's group recruited people of each group (\eqn{\theta}{theta});
#' select a participant at random from the sample members of that group. This
#' participant becomes the new active recruiter.
#'
#' 3. If the sample size has not been met, go to step 2.
#'
#' \strong{Tree bootstrap ("Treeboot").} A nonparametric procedure that
#' resamples the observed recruitment tree itself rather than fitting a model to
#' it. New trees are grown by resampling the seeds and, at each recruiter,
#' resampling with replacement from that recruiter's own observed recruits, so
#' that the branching pattern and wave structure of the real study are
#' preserved; each simulated participant carries a real participant's trait and
#' degree. Because it assumes no particular model for how traits mix or how
#' degree is distributed, it reflects the actual clustering present in the data,
#' but it can only reproduce structure that is present in the observed tree and
#' needs a reasonably large sample to resample from. It applies to categorical
#' outcomes.
#'
#' Bootstrap steps (to generate one bootstrap sample):
#'
#' 0. Record the observed recruitment tree (who recruited whom), together with
#' each participant's trait and degree.
#'
#' 1. Resample the seeds, with replacement, from the observed seeds.
#'
#' 2. At each active recruiter, resample its recruits, with replacement, from
#' that recruiter's own observed recruits (so the number of recruits per
#' recruiter is preserved); each resampled participant carries a real
#' participant's trait and degree.
#'
#' 3. Continue down the tree, starting further trees from resampled seeds as
#' needed, until the sample size is reached, trimming the final wave if it
#' overshoots.
#' 
#' Because sampling in the bootstrap is done with replacement, this method tends to
#' have wider confidence intervals. RDS recruitment is done without replacement
#' (i.e. subjects are only recruited once) in practice.
#'
#' \strong{Gile's successive sampling (SS) bootstrap.} A model-based procedure
#' designed for studies that sample a non-negligible fraction of the hidden
#' population. From the sample it reconstructs an estimated population of the
#' specified size, then repeatedly draws fresh RDS samples from that synthetic
#' population by successive sampling: individuals are selected with probability
#' proportional to their degree and without replacement, so the available pool
#' is depleted as the sample grows -- exactly as in a real without-replacement
#' RDS. Branching is matched to the observed number of recruits per participant,
#' and trait mixing to the observed cross-group tie structure. Because it models
#' the finite population and the without-replacement depletion, its intervals
#' correctly widen as the sample fraction increases. It is the companion to the
#' Gile SS estimator.
#'
#' Bootstrap steps (to generate one bootstrap sample):
#'
#' 0. Using the successive-sampling inclusion probabilities, reconstruct an
#' estimated population of the specified size -- how many people fall in each
#' degree-by-trait class -- and estimate the cross-group tie (mixing) structure.
#'
#' 1. Draw the seeds from this population with probability proportional to
#' degree, without replacement.
#'
#' 2. Give each recruiter a number of recruits drawn (with replacement) from the
#' observed counts of returned coupons.
#'
#' 3. For each recruit, draw its group from the mixing structure given its
#' recruiter's group, select an individual of that group with probability
#' proportional to degree, and remove that individual from the population
#' (depletion, i.e. without replacement).
#'
#' 4. Repeat until the sample size is reached.
#'
#' \strong{Homophily configuration graph (HCG) bootstrap.} A model-based
#' procedure built on a network model that adds homophily and differential
#' activity. It estimates the population's group-to-group tie mixing and each
#' group's degree distribution, then regenerates recruitment: seeds are assigned
#' traits according to the estimated population composition, each subsequent
#' recruit is assigned a trait from the estimated mixing given its recruiter's
#' trait, and each participant is assigned a degree from the relevant group's
#' distribution (higher-degree people more likely). When the sample is a
#' non-negligible fraction of the population, degrees are drawn without
#' replacement so the pool depletes. The resulting intervals reflect
#' homophily, differing activity levels between groups, short recruitment
#' chains, and seed bias. When the sample is a small fraction of the population
#' this procedure coincides with the Salganik bootstrap. It is the companion to
#' the HCG estimator.
#'
#' Bootstrap steps (to generate one bootstrap sample):
#'
#' 0. Fit the sample to a homophily configuration graph model: estimate the group-to-group tie mixing
#' (\eqn{\theta}{theta}), the estimated population composition of traits, and
#' the degree distribution within each group.
#'
#' 1. Assign each seed a trait drawn from the estimated population composition,
#' and a degree drawn in proportion to degree from that group.
#'
#' 2. Following the observed recruitment tree, assign each recruit a trait drawn
#' from \eqn{\theta}{theta} given its recruiter's trait, and a degree drawn in
#' proportion to degree from that group. When the sample is a non-negligible
#' fraction of the population, remove that individual from the pool (depletion,
#' i.e. without replacement).
#'
#' 3. Continue until every position in the observed tree is filled.
#'
#' \strong{Choosing a procedure.} In most analyses the default tied to the
#' estimator is appropriate. Broadly: SRS is a naive reference; the Gile SS bootstrap 
#' is preferred when that
#' fraction is large; and the HCG bootstrap is preferred when homophily,
#' differential activity, or short recruitment chains are a concern.
#'
#' @references
#' Salganik, M. J. (2006). Variance estimation, design effects, and sample size
#' calculations for respondent-driven sampling. \emph{Journal of Urban Health},
#' 83(S1), 98-112.
#'
#' Baraff, A. J., McCormick, T. H., and Raftery, A. E. (2016). Estimating
#' uncertainty in respondent-driven sampling using a tree bootstrap method.
#' \emph{Proceedings of the National Academy of Sciences}, 113(51),
#' 14668-14673. <doi:10.1073/pnas.1617258113>
#'
#' Gile, K. J. (2011). Improved inference for respondent-driven sampling data
#' with application to HIV prevalence estimation. \emph{Journal of the American
#' Statistical Association}, 106, 135-146.
#'
#' Gile, K. J. and Handcock, M. S. (2010). Respondent-driven sampling: An
#' assessment of current methodology. \emph{Sociological Methodology}, 40,
#' 285-327. <doi:10.1111/j.1467-9531.2010.01223.x>
#'
#' Gile, K. J., Beaudry, I. S., Handcock, M. S. and Ott, M. Q. (2018). Methods
#' for inference from respondent-driven sampling data. \emph{Annual Review of
#' Statistics and Its Application}, 5, 65-93.
#' <doi:10.1146/annurev-statistics-031017-100704>
#'
#' Fellows, I. E. (2019). Respondent-driven sampling and the homophily
#' configuration graph. \emph{Statistics in Medicine}, 38, 131-150.
#'
#' @seealso \code{\link{RDS.bootstrap.intervals}},
#' \code{\link{control.rds.estimates}}
#' @name RDS.bootstrap.mechanics
NULL
