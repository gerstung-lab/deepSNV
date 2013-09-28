# deepSNV class definiton
# 
# Author: gemoritz
###############################################################################


#' deepSNV class.
#' 
#' This class stores the contents of the deepSNV test. It is typically initialized with \code{\link{deepSNV}}.
#' This class has the following slots:
#' \describe{
#' \item{p.val}{The P-values of the test.}
#' \item{test}{A matrix with the nucleotide counts in the test experiment. The column names of the nucleotide counts are A, T, C, G, - for the positivie strand and a, t, c, g, _ for the reverse.} 
#' \item{control}{A matrix with the nucleotide counts in the control experiment. The column names must be the same as for the test.}
#' \item{coordinates}{A \code{\link{data.frame}} with the genomic coordinates chr and pos, and other columns, if desired.}
#' \item{dirichlet.prior}{A matrix with the nucleotide-specific Dirichlet prior}
#' \item{pseudo.count}{The pseudo count if used)}
#' \item{alternative}{A string with the alternative used in the test.}
#' \item{nucleotides}{A character vector with the nucleotides tested.}
#' \item{regions}{A \code{\link{data.frame}} with columns chr, start, and stop.}
#' \item{files}{A list with two entries test and control storing the filenames (if the object was initialized from two bam-files).}
#' \item{combine.method}{The method for combining p-values as a character string.}
#' \item{model}{The statistical model, either bin for binomial, or betabin for beta-binomial}
#' \item{over.dispersion}{If the model is beta-binomial, the first parameter for the beta-binomial model, which is shared across sites.}
#' \item{call}{The last function call to deepSNV.}
#' \item{log.lik}{The log likelihood of the data under the null hypothesis. (Excluding zeros on the opposite site under a one-sided test.)}
#' }
#' @seealso \code{\link{deepSNV}}
#' @rdname deepSNV-class
#' @name deepSNV-class
#' @example inst/example/deepSNV-example.R
#' @author Moritz Gerstung
#' @exportClass deepSNV
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
				log.lik = "numeric",
				pseudo.count="numeric"
	)
)

