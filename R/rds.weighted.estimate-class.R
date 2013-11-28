################################################################################
# Class definition for "rds.wieghted.estimate".   This class will be returned as the
# result of computing RDS estimates with no confidence interval.    
################################################################################

setClass("rds.weighted.estimate",
		representation(estimate="numeric",
				weights="numeric",
				outcome.variable="character",
				mean.group.weights="numeric",
				weight.type="character",
				subset="logical"))


setMethod("initialize",
		"rds.weighted.estimate",
		function(.Object,
				estimate,
				weights,
				outcome.variable,
				weight.type,
				subset) {
			
			.Object@estimate <- estimate
			.Object@weights <- weights
			.Object@outcome.variable <- outcome.variable
			.Object@weight.type <- weight.type
			if(is.null(subset))
				subset <- rep(TRUE,length(weights))
			.Object@subset <- subset
			return(.Object)
		})

# A "friendly" constructor.  

rds.weighted.estimate <- function(estimate, weights,
		outcome.variable,weight.type,subset){
	new("rds.weighted.estimate", estimate, weights,
			outcome.variable,weight.type,subset)
}


setMethod("print",signature="rds.weighted.estimate",
		definition=function(x){
			cat(x@weight.type,"Estimate for",x@outcome.variable,"\n")
			if(!is.null(attr(x@estimate,"se"))){
				se <- attr(x@estimate,"se")
				attr(x@estimate,"se") <- NULL
				print(x@estimate)
				if(se > 0){
					cat("Estimated standard error","\n")
					print(se)
				}
			}else{ 
				print(x@estimate)
			}
		}
)

setMethod("show",signature="rds.weighted.estimate",
		definition=function(object){
			cat(object@weight.type,"Estimate for",object@outcome.variable,"\n")
			if(!is.null(attr(object@estimate,"se"))){
				se <- attr(object@estimate,"se")
				attr(object@estimate,"se") <- NULL
				print(object@estimate)
				if(se > 0){
					cat("Estimated standard error","\n")
					print(se)
				}
			}else{ 
				print(object@estimate)
			}
		}
)




