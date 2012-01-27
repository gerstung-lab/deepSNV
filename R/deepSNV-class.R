# deepSNV class definiton
# 
# Author: gemoritz
###############################################################################


#' deepSNV class.
#' 
#' This class stores the contents of the deepSNV test. It is typically initialized with \code{\link{deepSNV-methods}}.
#' @author gemoritz
#' @export
#' @slot p.val The P-values of the test.
#' @slot test A matrix with the nucleotide counts in the test experiment. The column names of the nucleotide counts are A, T, C, G, - for the positivie strand and a, t, c, g, _ for the reverse. 
#' @slot control A matrix with the nucleotide counts in the control experiment. The column names must be the same as for the test.
#' @slot coordinates A \code{\link{data.frame}} with the genomic coordinates chr and pos, and other columns, if desired.
#' @slot dirichlet.prior A matrix with the nucleotide-specific Dirichlet prior
#' @slot alternative A string with the alternative used in the test.
#' @slot nucleotides A character vector with the nucleotides tested.
#' @slot regions A \code{\link{data.frame}} with columns chr, start, and stop.
#' @slot files A list with two entries test and control storing the filenames (if the object was initialized from two bam-files).
#' @slot combine.method The method for combining p-values as a character string.
#' @slot model The statistical model, either bin for binomial, or betabin for beta-binomial
#' @slot over.dispersion If the model is beta-binomial, the first parameter for the beta-binomial model, which is shared across sites.
#' @slot call The last function call to deepSNV.
#' @slot log.lik The log likelihood of the data under the null hypothesis. (Excluding zeros on the opposite site under a one-sided test.)
#' @seealso \code{\link{deepSNV-methods}}
#' @aliases deepSNV-class
#' @example deepSNV/inst/example/deepSNV-example.R
setClass("deepSNV",
		representation=representation(
				p.val = "matrix",
				test = "matrix",
				control = "matrix",
				coordinates = "data.frame",
				dirichlet.prior = "matrix",
				alternative = "character",
				nucleotides = "vector",
				regions = "data.frame",
				files = "list",
				combine.method = "character",
				over.dispersion = "numeric",
				model = "character",
				call = "call",
				log.lik = "numeric"
	)
)

