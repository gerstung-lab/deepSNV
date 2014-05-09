# The cormorant algorithm computes variants using a precomputed reference.
# 
# Author: mg14
###############################################################################

#' Test against normals
#' 
#' No strand at this stage.
#' Added 29/11/2013
#' @param tumorCounts
#' @param errRates
#' @return array
#' @noRd
#' @export 
#' @author mg14
bbNormal <- function(tumorCounts, errRates, ...){	
	ncol = dim(tumorCounts)[3]/2
	
	countForward = tumorCounts[,,1:ncol, drop=FALSE]
	countBackward = tumorCounts[,,1:ncol + ncol, drop=FALSE]
	
	coverageForward = rep(rowSums(countForward, dims=2), dim(countForward)[3])
	coverageBackward = rep(rowSums(countBackward, dims=2), dim(countBackward)[3])
	
	pForward <- pbetabinom(countForward, coverageForward, errRates$nu.fw, (1-errRates$rho)/errRates$rho)
	pBackward <- pbetabinom(countBackward, coverageBackward, errRates$nu.bw, (1-errRates$rho)/errRates$rho)
	
	
}

#' Cormorant algorithm
#' @param tumorFiles 
#' @param normalFiles 
#' @param regions 
#' @param nucleotides 
#' @param mc.cores 
#' @param q 
#' @param ... 
#' @return A data.frame with p-values
#' 
#' @author mg14
#' @export
cormorant <- function(tumorFiles, normalFiles, regions, nucleotides, mc.cores=1, q=25, ...){
	tumorCounts <- loadAllData(tumorFiles, regions, mc.cores=mc.cores, q=q)
	normalCounts <- loadAllData(normalFiles, regions, mc.cores=mc.cores, q=q)
	errRates <- bbb(normalCounts, return.value = "err", truncate = 1 + .Machine$double.eps, ...)
	p <- bbNormal(tumorCounts, errRates, ...)
	ix <- sapply( c("A","T","C","G","-"), `==`,nucleotides)
	p <- t(matrix(p[rep(ix, each=nrow(p))], nrow=nrow(p)))
	colnames(p) <- tumorFiles
	result <- cbind(p, mu=errRates$nu[ix], rho=errRates$rho[ix]) 
	w <- which(ix, arr.ind=TRUE)
	return(result[order(w[,1], w[,2]),])
}

#' Beta-binomial distribution
#' @param x Counts, integer
#' @param n Size, integer
#' @param mu Probability, nonnegative numeric
#' @param rho Dispersion, numeric (0,1)
#' @param log Logical. Should logarithmic probabilities be returned?
#' @return Probability
#' 
#' @author mg14
#' @export
dbbinom <- function(x,n,mu,rho,log=FALSE){
	disp <- (1-rho)/rho
	res <- logbb(x,n,mu*disp, disp) + lgamma(n+1) - lgamma(x+1) -lgamma(n-x+1)
	if(!log)
		res <- exp(res)
	return(res)
}

#' Beta-binomial cumulative density
#' 
#' @param q Counts, integer
#' @param n Size, integer
#' @param mu Probability, nonnegative numeric
#' @param rho Dispersion, numeric (0,1)
#' @param log.p Logical. Should logarithmic probabilities be returned?
#' @return Cumulative probability
#' 
#' @author mg14
#' @export
pbbinom <- function(q,n,mu,rho, lower.tail=FALSE, log.p=FALSE){
	p <- q
	p[] <- sapply(1:length(q), function(i) {
				if(lower.tail)
					w <- 0:q[i]
				else
					w <- n[i]:q[i]
				sum(exp(dbbinom(w, n[i], mu[i],rho[i], log=TRUE)))
			})
	return(p)
}