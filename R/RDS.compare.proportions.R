
#' Compares the rates of two variables against one another.
#' @param first.interval An \code{rds.interval.estimate} object fit with either "Gile" or "Salganik" uncertainty.
#' @param second.interval An \code{rds.interval.estimate} object fit with either "Gile" or "Salganik" uncertainty.
#' @param M The number of bootstrap resamplings to use
#' @details This function preforms a bootstrap test comparing the 
#' data(faux)
#' int1 <- RDS.bootstrap.intervals(faux, outcome.variable=c("X"), weight.type="RDS-II", 
#' 		uncertainty="Salganik", N=1000, number.ss.samples.per.iteration=1000, 
#' 		confidence.level=0.95, number.of.bootstrap.samples=500)
#' int2 <- RDS.bootstrap.intervals(faux, outcome.variable=c("Y"), 
#' 		weight.type="RDS-II", uncertainty="Salganik", N=1000, number.ss.samples.per.iteration=1000, confidence.level=0.95, number.of.bootstrap.samples=500)
#' RDS.compare.proportions(int1,int2)
#' @export
RDS.compare.proportions <- function(first.interval,second.interval,M=10000){
	
	bsresult1 <- attr(first.interval,"bsresult") 
	boots1 <- if(is.matrix(bsresult1)) bsresult1 else bsresult1$bsests
	levs1 <- colnames(boots1)
	nlev1 <- length(levs1)
	bs1 <- boots1[sample(1:nrow(boots1),replace=TRUE,size=M),]
	tmp <- first.interval$estimate
	obs <- tmp[match(names(tmp),levs1)]
	bs1 <- sweep(bs1,2,colMeans(bs1)-obs)

	bsresult2 <- attr(second.interval,"bsresult") 
	boots2 <- if(is.matrix(bsresult2)) bsresult2 else bsresult2$bsests
	levs2 <- colnames(boots2)
	nlev2 <- length(levs2)
	bs2 <- boots2[sample(1:nrow(boots2),replace=TRUE,size=M),]
	tmp <- second.interval$estimate
	obs <- tmp[match(names(tmp),levs2)]
	bs2 <- sweep(bs2,2,colMeans(bs2)-obs)

	res <- matrix(NA,ncol=nlev2,nrow=nlev1)
	colnames(res) <- levs2
	rownames(res) <- levs1
	for(i in 1:nlev1){
		for(j in 1:nlev2){
			bs <- bs1[,i] - bs2[,j]
			pval <- 2*min(c(mean(bs<=0),mean(bs>=0)))
			res[i,j] <- pval
		}
	}
	res <- as.data.frame(res)
	class(res) <- c("pvalue.table","data.frame")
	res
}

#' Displays a pvalue.table
#' @param x a pvalue.table object
#' @param ... additional parameters passed to print.data.frame.
#' @export
#' @method print pvalue.table
print.pvalue.table <- function(x,...){
	cat("P-Value table comparing rates:\n")
	print.data.frame(x,...)
}
