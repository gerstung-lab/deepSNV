#### C implementations of the beta-binomial distribution


#' Beta-binomial probability distribution
#' @param x Counts
#' @param size Size
#' @param prob Probability
#' @param dispersion Dispersion rho/(1-rho)
#' @param log Should the natural logarithm be returned? Default=FALSE
#' @return Probability
#' 
#' @author mg14
#' @export
dbetabinom = function(x, size, prob, dispersion, log=FALSE){
	l = max(length(x),length(size),length(prob), length(dispersion))
	d = numeric(l)
	result = .C("dbetabinom",
			d,
			as.integer(l),
			as.integer(x),
			as.integer(length(x)),
			as.integer(size),
			as.integer(length(size)),
			as.numeric(prob),
			as.integer(length(prob)),
			as.numeric(dispersion),
			as.integer(length(dispersion)),
			as.integer(log),
			PACKAGE="deepSNV"
	)[[1]]
	return(result)
}

#' Cumulative beta-binomial probability distribution
#' 
#' The underlying implementation is written in C
#' @param x Counts
#' @param size Sample size
#' @param prob Probability
#' @param dispersion Dispersion rho/(1-rho)
#' @param log Should the natural logarithm be returned? Default=FALSE
#' @return Probability
#' 
#' @author mg14
#' @export
pbetabinom = function(x, size, prob, dispersion, log=FALSE){
	l = max(length(x),length(size),length(prob), length(dispersion))
	p = numeric(l)
	result = .C("pbetabinom",
			p,
			as.integer(l),
			as.integer(x),
			as.integer(length(x)),
			as.integer(size),
			as.integer(length(size)),
			as.numeric(prob),
			as.integer(length(prob)),
			as.numeric(dispersion),
			as.integer(length(dispersion)),
			as.integer(log),
			PACKAGE="deepSNV"
	)[[1]]	
	return(result)
}