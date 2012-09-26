# Experimental functions in deepSNV. Functions in this file are not guaranteed
# to work! They are not exported from the namespace.
# Author: Moritz Gerstung
###############################################################################

#' Test for PCR errors.
#' 
#' A test for PCR errors based on the Luria-Delbruck model of error propagation in exponentially growing populations. The test makes use of the fact that in the LD model
#' the variance scales as \eqn{\sigma^2 = 2 \mu T}, where \eqn{\mu} is the mutation rate and \eqn{T} the number of cycles. Moreover, if there are intially \eqn{N_0 \gg 1} molecules, 
#' one may use a Normal approximation for the distribution of errors with mean \eqn{\mu} and variance \eqn{\sigma^2/N_0}. 
#' @param test Either an \code{\link{deepSNV-class}} object or a named matrix with nucleotide counts.
#' @param control Missing if test is an \code{link{deepSNV-class}} object, otherwise a matrix with nucleotide counts.
#' @param alpha A parameter controlling the variance. Defined as the number of cycles tims the initial number of molecules. Default alpha = 33*1000.
#' @return A matrix with p-values for test and control of class matrix, otherwise an deepSNV object.
#' @author Moritz Gerstung
#' @noRd
setGeneric("PCRTest",
		function(test, control, ...) {
			standardGeneric("PCRTest")
		})

#' Test for PCR errors.
#' @noRd
setMethod("PCRTest",
		signature(test="matrix", control="matrix"),
		function(test, control, alpha = 33*1000) {
			test.rf <- RF(test+.5)
			control.rf <- RF(control+.5)
			both.rf <- RF(test+control+.5)
			v <- 2  * both.rf / alpha #Variance in LD-model #V = f/cycles/N0
			2*pnorm(abs(test.rf-control.rf)/(sqrt(v*2)), lower.tail=FALSE)
		})

#' Test for PCR errors.
#' @noRd
setMethod("PCRTest",
		signature(test="deepSNV", control="missing"),
		function(test, control, ...) {
			p.val <- PCRTest(test(test, total=T), control(test, total=T), ...)
			result <- initialize(test, p.val = pmax(p.val, test@p.val))
			return(result)
		}
)




#' Estimate overdispersion.
#' 
#' This method fits the overdispersion by minimizing the squared distance of the empirical cumulative distribution of the p-values of the \code{\link{deepSNV}} test
#' to a uniform distribution.
#' @note This experimental feature works only if the majority of sites contains no low-frequency variants.
#' @param test Either an \code{\link{deepSNV-class}} object or a named matrix with nucleotide counts.
#' @param control Missing if test is an \code{link{deepSNV-class}} object, otherwise a matrix with nucleotide counts.
#' @param weight Character string denoting how weights should be computed.
#' @param ... Additinal arguments passed to \code{\link{deepSNV}}
#' @author Moritz Gerstung
#' @noRd
setGeneric("overDispersion",
		function(test, control, weight, ...) {
			standardGeneric("overDispersion")
		})

#' Estimate overdispersion.
#' @noRd
setMethod("overDispersion",
		signature(test="matrix", control="matrix"),
		function(test, control, weight = c("none","geometric"), ...) {
			weight = match.arg(weight)
			optimize(function(x) {
						p = na.omit(as.vector(deepSNV(test,control, over.dispersion = x, ...)@p.val)) + .Machine$double.eps
						l <- length(p)
						rsq = (log(p) -log(rank(p)/l))^2
						if(weight == "geometric")
							rsq = rsq/rank(p)	
						sum(rsq)
					}, 
					interval=c(1,10))
		})

#' Estimate overdispersion.
#' @noRd
setMethod("overDispersion", 
		signature(test="deepSNV", control="missing"),
		function(test, control, weight = c("none","geometric"), ...){
			weight = match.arg(weight)
			overDispersion(test@test, test@control, weight, alternative=test@alternative, combine.method=test@combine.method, dirichlet.prior=test@dirichlet.prior,...)
		}
)


