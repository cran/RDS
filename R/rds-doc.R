
#' This package provides functionality for carrying out estimation
#'     with data collected using Respondent-Driven Sampling. This includes
#'     Heckathorn's RDS-I and RDS-II estimators as well as Gile's Sequential
#'     Sampler estimator.
#' @import ggplot2 scales reshape2 gridExtra methods
#' @docType package
#' @name RDS
#' @useDynLib RDS
NULL


#' A Simulated RDS Data Set
#' @description This is a faux set used to demonstrate RDS functions and analysis.
#' @docType data
#' @keywords datasets
#' @format An rds.data.frame object
#' @references Gile, Krista J., Handcock, Mark S., 2010 \emph{Respondent-driven Sampling: An Assessment of Current Methodology},  \emph{Sociological Methodology}, 40, 285-327. 
#' @seealso \code{\link{fauxsycamore}}, \code{\link{fauxmadrona}}
#' @examples 
#' data(faux)
#' RDS.I.estimates(rds.data=faux,outcome.variable='X')
#' @name faux
NULL


#' A Simulated RDS Data Set with no seed dependency
#' 
#' This is a faux set used to illustrate how the estimators perform under
#' different populations and RDS schemes.
#' 
#' The population had N=1000 nodes.  In this case, the sample size is 500 so
#' that there is a relatively small sample fraction (50\%). There is homophily
#' on disease status (R=5) and there is differential activity by disease status
#' whereby the infected nodes have mean degree twice that of the uninfected
#' (w=1.8).
#' 
#' In the sampling, the seeds are chosen randomly from the full population, so
#' there is no dependency induced by seed selection.
#' 
#' Each sample member is given 2 uniquely identified coupons to distribute to
#' other members of the target population in their acquaintance.  Further each
#' respondent distributes their coupons completely at random from among those
#' they are connected to.
#' 
#' 
#' 
#' Here are the results for this data set and the sister \code{fauxsycamore}
#' data set:
#' 
#' \tabular{rlllllll}{ \bold{Name} \tab \bold{City} \tab \bold{Type} \tab
#' \bold{Mean} \tab \bold{RDS I (SH)} \tab \bold{RDS II (VH)} \tab \bold{SS}
#' \tab \bold{MA}\cr fauxsycamore \tab Oxford\tab seed dependency, 70\% \tab
#' 0.2408 \tab 0.1087 \tab 0.1372 \tab 0.1814 \tab 0.1843\cr fauxmadrona \tab
#' Seattle\tab no seed dependency, 50\% \tab 0.2592 \tab 0.1592 \tab 0.1644
#' \tab 0.1941 \tab 0.1978\cr fauxbanksia \tab Perth\tab ? \tab ? \tab ? \tab ?
#' \tab ? \tab ?  }
#' 
#' Even with only 50\% sample, the VH is substantially biased , and the SS and
#' MA do much better.  We expect the MA to perform about as well as the SS (and
#' it does).
#' 
#' 
#' @name fauxmadrona
#' @aliases fauxmadrona fauxmadrona.network
#' @docType data
#' @format An \code{rds.data.frame}
#' @seealso \code{\link{fauxsycamore}}, \code{\link{faux}}
#' @references Gile, Krista J., Handcock, Mark S., 2010 \emph{Respondent-driven
#' Sampling: An Assessment of Current Methodology}, \emph{Sociological
#' Methodology}, 40, 285-327.
#' @source It is
#' \code{/net/proj/rdsworkinggroup/kgile/sims/rdsnetsamps/rdssamp1000_5_2__-1.RData}
#' \cr With networks at: \cr
#' \code{/net/proj/rdsworkinggroup/kgile/sims/rdsnetsims/nets1000_5_21.RData}
#' \cr It is network 1 and the extraction code is in the directory for YesYes
#' on mosix. Look for "extract1.R" \cr The original network is included as
#' \code{fauxmadrona.network} as a \code{network} object.  \cr The data set
#' also includes the \code{data.frame} of the RDS data set as
#' \code{fauxmadrona}.  \cr Use \code{data(package="RDS")} to get a full list
#' of datasets.
#' 
#' @keywords datasets
NULL




#' A Simulated RDS Data Set with extreme seed dependency
#' 
#' This is a faux set used to demonstrate RDS functions and analysis.  The
#' population had N=715 nodes.  In this case, the sample size is 500 so that
#' there is a relatively large sample fraction (70\%). There is homophily on
#' disease status (R=5) and there is differential activity by disease status
#' whereby the infected nodes have mean degree twice that of the uninfected
#' (w=1.8).
#' 
#' In the sampling the seeds are chosen randomly from the infected population,
#' so there is extreme dependency induced by seed selection.
#' 
#' Each sample member is given 2 uniquely identified coupons to distribute to
#' other members of the target population in their acquaintance.  Further each
#' respondent distributes their coupons completely at random from among those
#' they are connected to.
#' 
#' With 70\% sample, the VH is substantially biased, so the SS (and presumably
#' MA) do much better.  We expect the MA to perform a bit better than the SS.
#' 
#' It is network 702 and its sample from YesYes on mosix. Look for
#' "extract702.R" \cr The original network is included as
#' \code{fauxsycamore.network} as a \code{network} object.  \cr The data set
#' also includes the \code{data.frame} of the RDS data set as
#' \code{fauxsycamore}.  \cr Use \code{data(package="RDS")} to get a full list
#' of datasets.
#' 
#' @name fauxsycamore
#' @aliases fauxsycamore fauxsycamore.network
#' @docType data
#' @format An rds.data.frame plus the original network as a network object
#' @seealso \code{\link{faux}}, \code{\link{fauxmadrona}}
#' @references Gile, Krista J., Handcock, Mark S., 2009.
#' \emph{Respondent-driven Sampling: An Assessment of Current Methodology},
#' \emph{Sociological Methodology}, 40, 285-327.
#' @keywords datasets
NULL




