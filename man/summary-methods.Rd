\name{summary-methods}
\alias{summary-methods}
\title{Summary of a deepSNV object...}
\description{Summary of a deepSNV object}
\details{\code{summary-methods}: Tabularize significant SNVs by evalutating the p-values of the \code{\link{deepSNV}} test.

}
\value{\code{summary-methods}: A data.frame with the following columns:
\item{chr}{The chromosome}
\item{pos}{The position (1-based)}
\item{ref}{The reference (consensus) nucleotide}
\item{var}{The variant nucleotide}
\item{p.val}{The (corrected) p-value}
\item{freq.var}{The relative frequency of the SNV}
\item{sigma2.freq.var}{The estimated variance of the frequency}
\item{n.tst.fw}{The variant counts in the test experiment, forward strand}
\item{cov.tst.fw}{The coverage in the test experiment, forward strand}
\item{n.tst.bw}{The variant counts in the test experiment, backward strand}
\item{cov.tst.bw}{The coverage in the test experiment, backward strand}
\item{n.ctrl.fw}{The variant counts in the control experiment, forward strand}
\item{cov.ctrl.fw}{The coverage in the control experiment, forward strand}
\item{n.ctrl.bw}{The variant counts in the control experiment, backward strand}
\item{cov.ctrl.bw}{The coverage in the control experiment, backward strand}
\item{raw.p.val}{The raw p-value}

}
\author{gemoritz}
\arguments{\item{object}{A \code{\link{deepSNV-class}} object.}
\item{sig.level}{The desired significance level.}
\item{adjust.method}{The adjustment method for multiple testing corrections. See \code{\link{p.adjust}} for details. Set to NULL, for no adjustment. Default "bonferroni".}
\item{fold.change}{The minimal fold change required of the relative frequency. Default 1.}
}
\examples{## Short example with 2 SNVs at frequency ~10%
regions <- data.frame(chr="B.FR.83.HXB2_LAI_IIIB_BRU_K034", start = 3120, stop=3140)
ex <- deepSNV(test = system.file("extdata", "test.bam", package="deepSNV"), control = system.file("extdata", "control.bam", package="deepSNV"), regions=regions, q=10)
show(ex)   # show method
plot(ex)   # scatter plot
summary(ex)   # summary with significant SNVs
ex[1:3,]   # subsetting the first three genomic positions
tail(test(ex, total=TRUE))   # retrieve the test counts on both strands
tail(control(ex, total=TRUE))

## Not run: Full example with ~ 100 SNVs. Requires an internet connection, but try yourself.
# regions <- data.frame(chr="B.FR.83.HXB2_LAI_IIIB_BRU_K034", start = 2074, stop=3585)
# HIVmix <- deepSNV(test = "http://www.bsse.ethz.ch/cbg/software/deepSNV/data/test.bam", control = "http://www.bsse.ethz.ch/cbg/software/deepSNV/data/control.bam", regions=regions, q=10)
data(HIVmix) # attach data instead..
show(HIVmix)
plot(HIVmix)
head(summary(HIVmix))}
\alias{summary,deepSNV-method}
\alias{summary}
\usage{\S4method{summary}{deepSNV}(object, sig.level=0.05, adjust.method="bonferroni", fold.change=1)
}

