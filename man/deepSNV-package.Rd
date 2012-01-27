\name{deepSNV-package}
\alias{deepSNV-package}
\alias{deepSNV-package}
\title{Estimation of single nucleotide variants from paired deep sequencing experiments}
\description{Estimation of single nucleotide variants from paired deep sequencing experiments}
\details{\tabular{ll}{
Package: \tab deepSNV\cr
Type: \tab Package\cr
Version: \tab 0.99.2\cr
Date: \tab 2012-01-20\cr
License: \tab GPL3\cr
LazyLoad: \tab yes\cr
}

This packages provides algorithms for estimating SNVs and their frequencies from matched ultra-deep sequencing experiments of heterogeneous samples. 
It retrieves the nucleotide counts at each position and each strand and tests for differences between the two experiments with a likelihood ratio test using
either a binomial or beta-binomial model. 
The statistic can be tuned across genomic sites by a shared Dirichlet prior.}
\docType{package}
\author{Moritz Gerstung, ETH Zurich, D-BSSE \email{gemoritz@ethz.ch}}
\references{See our upcoming paper ...}
\keyword{package}
\seealso{\code{\link{deepSNV}}}
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

