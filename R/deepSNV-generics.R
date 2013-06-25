# Generic methods of the deepSNV package
# 
# Author: Moritz Gerstung
###############################################################################


#' Scatter plot of relative nucleotide frequencies.
#' 
#' This function plots the relative nucleotide frequencies of the test against the control experiment on a logarithmit scale. The color of the symbols
#' denotes the nucleotide, and the area of the circle is proportional to the \eqn{- \log}{-log} of the p-value.
#' @param x A deep SNV object. 
#' @param sig.level By default, p-values below sig.level are drawn as filled circles.
#' @param col Color of the nucleotides.
#' @param col.null Color of insignificant nucleotides.
#' @param cex.min The minimal size of the points.
#' @param xlab The x-axis label.
#' @param ylab The y-axis label.
#' @param pch The plotting symbol. Default = 16 (filled circle)
#' @param ... Additional arguments passed to plot.
#' @return NULL
#' @author Moritz Gerstung
#' @export
#' @method plot deepSNV
#' @example inst/example/deepSNV-example.R
plot.deepSNV <- function (x, sig.level = NULL, col = NULL, col.null = "grey", cex.min = 0.2, ylab="Relative Frequency in Test", xlab="Relative Frequency in Control", pch = 16, ...) 
		{	
			test = RF(test(x, total=T)+1)
			control = RF(control(x, total=T)+1)
			range = c(10^floor(log10(c(min(test, control)))),1)
			logp = log(x@p.val)/log(10)
			steps = pretty(-c(max(min(logp, na.rm=T),-50), max(logp, na.rm=T)),4)
			scale = max(steps) / 10
			if(is.null(col)) col = nt.col[,2]
			if(length(col.null)==1) col.null=rep(col.null, ncol(test))
			if(is.null(sig.level)) sig.level = 0.05/sum(!is.na(x@p.val))
			
			colors <- col[(1:length(x@p.val) -1) %/% nrow(x@p.val) +1  ]
			colors[(logp > log10(sig.level))& !is.na(logp)] <- col.null[(1:length(x@p.val) -1) %/% nrow(x@p.val) +1  ][(logp > log10(sig.level)) & !is.na(logp)]
			plot(control,
					test,
					log='xy', 
					xlim=range, 
					ylim=range, 
					ylab=ylab, 
					xlab=xlab,
					pch=pch,#16,#c(1,16)[(logp < log10(sig.level)) + 1], 
					col=colors, 
					cex=sqrt(pmin(-logp,50))/scale+cex.min)
			legend("bottomright",1, 
					colnames(slot(x,"test"))[1:5],
					pch=pch, 
					pt.cex=sqrt(2), 
					col=col[1:ncol(x@test)], 
					bty='n',
					title = "nt")
			usr = par("usr")
			legend("bottomright",0, steps,pch=pch, pt.cex=sqrt(steps)/scale + cex.min, bty='n', title="-log10 P", inset=c(0.1,0))
			
			abline(0,1, lty=2)
		}

#' Summary of a deepSNV object
#' 
#' Tabularize significant SNVs by evalutating the p-values of the \code{\link{deepSNV}} test. 
#' @param object A \code{\link{deepSNV-class}} object.
#' @param sig.level The desired significance level.
#' @param adjust.method The adjustment method for multiple testing corrections. See \code{\link{p.adjust}} for details. Set to NULL, for no adjustment. Default "bonferroni".
#' @param fold.change The minimal fold change required of the relative frequency. Default 1.
#' @param value String. The type of the returned object. Either "data.frame" for a \code{\link{data.frame}} (default) or "VCF" for an Extended\code{\link{VCF-class}} object.
#' @return If value="data.frame", a \code{\link{data.frame}} with the following columns:
#' \item{chr}{The chromosome}
#' \item{pos}{The position (1-based)}
#' \item{ref}{The reference (consensus) nucleotide}
#' \item{var}{The variant nucleotide}
#' \item{p.val}{The (corrected) p-value}
#' \item{freq.var}{The relative frequency of the SNV}
#' \item{sigma2.freq.var}{The estimated variance of the frequency}
#' \item{n.tst.fw}{The variant counts in the test experiment, forward strand}
#' \item{cov.tst.fw}{The coverage in the test experiment, forward strand}
#' \item{n.tst.bw}{The variant counts in the test experiment, backward strand}
#' \item{cov.tst.bw}{The coverage in the test experiment, backward strand}
#' \item{n.ctrl.fw}{The variant counts in the control experiment, forward strand}
#' \item{cov.ctrl.fw}{The coverage in the control experiment, forward strand}
#' \item{n.ctrl.bw}{The variant counts in the control experiment, backward strand}
#' \item{cov.ctrl.bw}{The coverage in the control experiment, backward strand}
#' \item{raw.p.val}{The raw p-value}
#' If value = "VCF", this functions returns a  \code{\link{VCF-class}} object with the following entries:
#' FIXED:
#' \item{REF}{Reference allele in control sample. Note that deletions in the control sample will be reported like insertions, e.g. if the consensus of the control is A,- at positions 1 and 2 
#' (relative to the reference) and the test was A,A, then this would be denoted as REF="A" and VAR="AA" with coordinate IRanges(1,2). This may cause ambiguities when the VCF object is written 
#' to text with writeVcf(), which discards the width of the coordinate, and this variant remains indistinguishable from an insertion to the _reference_ genome.}
#' \item{VAR}{Variant allele in test sample}
#' \item{QUAL}{-10*log10(raw.p.val)}
#' INFO:
#' \item{VF}{Variant frequency. Variant allele frequency in the test minus variant allele frequency in the control.}
#' \item{VFV}{Variant frequency variance. Variance of the variant frequency; can be thought of as confidence interval.}
#' GENO (one column for test and one column for control):
#' \item{FW}{Forward allele count}
#' \item{BW}{Backward allele count}
#' \item{DFW}{Forward read depth}
#' \item{DBW}{Backward read depth}
#' @author Moritz Gerstung
#' @exportMethod summary
#' @example inst/example/deepSNV-example.R
#' @rdname summary-methods
#' @docType methods
#' @name summary
if (!isGeneric("summary"))
	setGeneric("summary", function(object, ...)
				standardGeneric("summary"),
			package = "deepSNV")

#' Summary for deepSNV object
#' @rdname summary-methods
#' @docType methods
#' @aliases summary,deepSNV-method
setMethod("summary",
		signature = signature(object = "deepSNV"),
		function (object, sig.level = 0.05, adjust.method = "bonferroni", fold.change=1, value=c("data.frame","VCF")) 
		{
			value <- match.arg(value)
			.significantSNV(object, sig.level, adjust.method, fold.change, value=value)
		}
)


#' Subsetting for deepSNV objects.
#' @param x A \code{\link{deepSNV-class}} object.
#' @param i Row indeces.
#' @param j Column (nucleotide) indeces.
#' @return A \code{\link{deepSNV-class}} object.
#' @author Moritz Gerstung
#' @exportMethod `[`
#' @name Extract
#' @docType methods
#' @rdname Extract-methods
#' @name `[`
#' @aliases [,deepSNV,ANY,ANY-method
#' @examples data(HIVmix)
#' HIVmix[1:10,]
setMethod("[",
		signature = signature(x = "deepSNV"),
		function (x, i, j){
			if (is.character(i))
				stop("subscript 'i' must be numeric")
			if (any(is.na(i)))
				stop("subscript 'i' contains NA")
			#if(length(i) == nrow(x@test)){
			
			result = initialize(x, 
					p.val = slot(x, "p.val")[i,j, drop=FALSE],
					test = slot(x, "test")[i,j, drop=FALSE],
					control = slot(x,"control")[i,j, drop=FALSE],
					coordinates = slot(x,"coordinates")[i,,drop=FALSE],
					nucleotides = slot(x,"nucleotides")[j,drop=FALSE]
			)
			#}else{
			#	stop("undefined subset selected")
			#}
			return(result)
		}
)

#' Show method for deepSNV objects
#' @param object A \code{\link{deepSNV-class}} object.
#' @return NULL
#' @author Moritz Gerstung
#' @exportMethod show
#' @example inst/example/deepSNV-example.R
#' @docType methods
#' @aliases show,deepSNV-method
setMethod("show",
		signature = signature(object = "deepSNV"),
		function (object) 
		{
			cat("Data: ", nrow(object@test), "positions x ", ncol(object@test), "characters\n")
			cat("Model: ", object@model, "\n")
			cat("Alternative: ", object@alternative, "\n")
			cat("Combine Method: ", object@combine.method, "\n")
			cat("P-Values:\n")
			if(nrow(object@p.val) > 20){
				print(head(object@p.val))
				cat("...\n")
				print(tail(object@p.val))
			}
			else
				print(object@p.val)
		}
)


