# Package definitions for deepSNV
# 
# Author: Moritz Gerstung
###############################################################################

#' Detection of subclonal SNVs in deep sequencing experiments
#' 
#' This packages provides algorithms for detecting subclonal single nucleotide variants (SNVs) and their frequencies from ultra-deep sequencing data. 
#' It retrieves the nucleotide counts at each position and each strand from two .bam files and tests for differences between the two experiments with a likelihood ratio test using
#' either a binomial or and overdispersed beta-binomial model. 
#' The statistic can be tuned across genomic sites by a shared Dirichlet prior and there package provides procedures for normalizing sequencing data from different runs.
#' @name deepSNV-package
#' @docType package
#' @title Detection of subclonal SNVs in deep sequencing experiments
#' @author Moritz Gerstung, Wellcome Trust Sanger Institute, \email{moritz.gerstung@@sanger.ac.uk}
#' @references Gerstung M, Beisel C, Rechsteiner M, Wild P, Schraml P, Moch H, and Beerenwinkel N. Reliable detection of subclonal single-nucleotide variants in tumour cell populations. Nat Commun 3:811 (2012). \href{http://dx.doi.org/10.1038/ncomms1814}{DOI:10.1038/ncomms1814}.
#' @keywords package
#' @seealso \code{\link{deepSNV}}
#' @example inst/example/deepSNV-example.R
#' @import Rhtslib
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

#' Example prior
#' 
#' Prior from COSMIC v63 for the TP53 gene
#' @docType data
#' @name pi
#' @examples data("pi", package="deepSNV")
#' plot(pi[,1], type="h")
NA

#' Example count table
#' 
#' A table with counts of the HIVmix data set. Used for minimal unit testing.
#' @docType data
#' @name counts
#' @examples data("counts", package="deepSNV")
#' countsFromBam <- bam2R(file = system.file("extdata", "test.bam", package="deepSNV"), chr="B.FR.83.HXB2_LAI_IIIB_BRU_K034", start = 3120, stop=3140, q = 10)
#' all(counts == countsFromBam)
NA