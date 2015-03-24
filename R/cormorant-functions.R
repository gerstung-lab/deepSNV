# The cormorant algorithm computes variants using a precomputed reference.
# 
# Author: mg14
###############################################################################

#' Test against normals
#' 
#' No strand at this stage.
#' Added 29/11/2013
#' @param tumorCounts An array (dim: pos x nucleotides x samples) for the tumor counts
#' @param errRates An list with error rates for matching positions (as returned by bbb)
#' @param combine.method The method for combining the P-values in forward and reverse orientation.
#' @return array
#' @noRd
#' @export 
#' @author mg14
bbNormal <- function(tumorCounts, errRates, combine.method = "fisher"){	
	ncol = dim(tumorCounts)[3]/2
	
	countForward = tumorCounts[,,1:ncol, drop=FALSE]
	countBackward = tumorCounts[,,1:ncol + ncol, drop=FALSE]
	
	coverageForward = rep(rowSums(countForward, dims=2), dim(countForward)[3])
	coverageBackward = rep(rowSums(countBackward, dims=2), dim(countBackward)[3])
	
	pForward <- pbetabinom(countForward, coverageForward, errRates$nu.fw, (1-errRates$rho)/errRates$rho)
	pBackward <- pbetabinom(countBackward, coverageBackward, errRates$nu.bw, (1-errRates$rho)/errRates$rho)
	
	pTotal <- p.combine(pForward, pBackward, method=combine.method)
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
