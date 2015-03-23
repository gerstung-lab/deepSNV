
#' Beta-binomial probability distribution
#' @param x Counts
#' @param n Size
#' @param mu Probability
#' @param rho Dispersion. rho in (0,1)
#' @param log Return logarithmic values
#' @return d
#' 
#' @author mg14
#' @export
dbetabinom = function(x, n, mu, rho, log=FALSE){
	disp <- (1-rho)/rho
	l = max(length(x),length(n),length(mu), length(disp))
	d = numeric(l)
	result = .C("dbetabinom",
			d,
			as.integer(l),
			as.integer(x),
			as.integer(length(x)),
			as.integer(n),
			as.integer(length(n)),
			as.numeric(mu),
			as.integer(length(mu)),
			as.numeric(disp),
			as.integer(length(disp)),
			as.integer(log),
			PACKAGE="deepSNV"
	)[[1]]
	return(result)
}

#' Cumulative beta-binomial probability distribution
#' @param x Counts
#' @param n Sample size
#' @param mu Probability
#' @param rho Dispersion. rho in (0,1)
#' @param log Return logarithmic values
#' @return Probability
#' 
#' @author mg14
#' @export
pbetabinom = function(x, n, mu, rho, log=FALSE){
	disp <- (1-rho)/rho
	l = max(length(x),length(n),length(mu), length(disp))
	p = numeric(l)
	result = .C("pbetabinom",
			p,
			as.integer(l),
			as.integer(x),
			as.integer(length(x)),
			as.integer(n),
			as.integer(length(n)),
			as.numeric(mu),
			as.integer(length(mu)),
			as.numeric(disp),
			as.integer(length(disp)),
			as.integer(log),
			PACKAGE="deepSNV"
	)[[1]]	
	return(result)
}