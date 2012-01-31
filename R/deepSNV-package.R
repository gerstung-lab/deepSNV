# Package definitions for deepSNV
# 
# Author: gemoritz
###############################################################################

#' Estimation of single nucleotide variants from paired deep sequencing experiments
#' 
#' \tabular{ll}{
#' Package: \tab deepSNV\cr
#' Type: \tab Package\cr
#' Version: \tab 0.99.3\cr
#' Date: \tab 2012-01-20\cr
#' License: \tab GPL3\cr
#' LazyLoad: \tab yes\cr
#' }
#' 
#' This packages provides algorithms for estimating SNVs and their frequencies from matched ultra-deep sequencing experiments of heterogeneous samples. 
#' It retrieves the nucleotide counts at each position and each strand and tests for differences between the two experiments with a likelihood ratio test using
#' either a binomial or beta-binomial model. 
#' The statistic can be tuned across genomic sites by a shared Dirichlet prior.
#' @name deepSNV-package
#' @docType package
#' @title Estimation of single nucleotide variants from paired deep sequencing experiments
#' @author Moritz Gerstung, ETH Zurich, D-BSSE \email{gemoritz@@ethz.ch}
#' @references See our upcoming paper ...
#' @keywords package
#' @seealso \code{\link{deepSNV}}
#' 
#' @example deepSNV/inst/example/deepSNV-example.R
#' @import Rsamtools
#' @import methods
#' @import Biostrings
#' @importFrom graphics plot
#' @importFrom VGAM vglm
#' @importFrom IRanges as.data.frame
#' @useDynLib deepSNV
NA

.onLoad <- function(lib, pkg){
	library.dynam("deepSNV", pkg, lib)
	#require(Biostrings)
	#require(graphics, methods)
}

.onUnload <- function(libpath){
	library.dynam.unload("deepSNV", libpath)
}

#' Example .bam data and true SNVs.
#' 
#' Two .bam alignments as example data sets are downloaded remotely via http. Sequenced were a 1,512 nt fragment of the HIV genome and a mixture (90\% + 10\%) with another variants. 
#' The two sequences were confirmed by Sanger sequencing and stored in the table trueSNVs. 
#' @docType data
#' @examples data(HIVmix)
#' data(trueSNVs)
#' table(p.adjust(p.val(HIVmix), method="BH") < 0.05, trueSNVs)
#' @name trueSNVs
#' @aliases HIVmix
NA


#' Example phiX data
#' 
#' Data from two phiX experiments sequenced on a GAIIx.
#' @docType data
#' @name phiX
#' @examples data(phiX, package="deepSNV")
#' plot(phiX)
#' phiN <- normalize(phiX, round=TRUE)
#' plot(phiN)
NA


#' Example RCC data
#' 
#' Deep sequencing experiments of a renal cell carcinoma and healthy control tissue. 
#' @docType data
#' @name RCC
#' @examples data("RCC", package="deepSNV")
#' summary(RCC, adjust.method="bonferroni")[,1:6]
#' plot(RCC)
#' RCC.bb <- estimateDispersion(RCC, alternative="two.sided")
#' summary(RCC.bb, adjust.method="bonferroni")[,1:6]
#' plot(RCC.bb)
NA

