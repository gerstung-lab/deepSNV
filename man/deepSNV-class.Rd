\name{deepSNV-class}
\alias{deepSNV-class}
\title{deepSNV class.}
\description{deepSNV class.}
\details{This class stores the contents of the deepSNV test. It is typically initialized with \code{\link{deepSNV-methods}}.}
\author{gemoritz}
\seealso{\code{\link{deepSNV-methods}}}
\alias{deepSNV-class}
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
\section{Slots}{\describe{\item{\code{p.val}:}{(\code{\link{matrix}}) The P-values of the test.}\item{\code{test}:}{(\code{\link{matrix}}) A matrix with the nucleotide counts in the test experiment. The column names of the nucleotide counts are A, T, C, G, - for the positivie strand and a, t, c, g, _ for the reverse.}\item{\code{control}:}{(\code{\link{matrix}}) A matrix with the nucleotide counts in the control experiment. The column names must be the same as for the test.}\item{\code{coordinates}:}{(\code{\link{data.frame}}) A \code{\link{data.frame}} with the genomic coordinates chr and pos, and other columns, if desired.}\item{\code{dirichlet.prior}:}{(\code{\link{matrix}}) A matrix with the nucleotide-specific Dirichlet prior}\item{\code{alternative}:}{(\code{\link{character}}) A string with the alternative used in the test.}\item{\code{nucleotides}:}{(\code{\link{vector}}) A character vector with the nucleotides tested.}\item{\code{regions}:}{(\code{\link{data.frame}}) A \code{\link{data.frame}} with columns chr, start, and stop.}\item{\code{files}:}{(\code{\link{list}}) A list with two entries test and control storing the filenames (if the object was initialized from two bam-files).}\item{\code{combine.method}:}{(\code{\link{character}}) The method for combining p-values as a character string.}\item{\code{model}:}{(\code{\link{character}}) The statistical model, either bin for binomial, or betabin for beta-binomial}\item{\code{over.dispersion}:}{(\code{\link{numeric}}) If the model is beta-binomial, the first parameter for the beta-binomial model, which is shared across sites.}\item{\code{call}:}{(\code{\link{call}}) The last function call to deepSNV.}\item{\code{log.lik}:}{(\code{\link{numeric}}) The log likelihood of the data under the null hypothesis. (Excluding zeros on the opposite site under a one-sided test.)}}}
\section{Methods}{\describe{\item{\code{\link[=summary,deepSNV-method]{summary}}}{\code{signature(object = "deepSNV")}: Summary for deepSNV object}\item{\code{\link[=[,deepSNV-method]{[}}}{\code{signature(x = "deepSNV")}: Subsetting for deepSNV objects.}\item{\code{\link[=show,deepSNV-method]{show}}}{\code{signature(object = "deepSNV")}: Show method for deepSNV objects}\item{\code{\link[=deepSNV,deepSNV,missing-method]{deepSNV}}}{\code{signature(test = "deepSNV", control = "missing")}: Method for signature deepSNV,missing}\item{\code{\link[=estimateDirichlet,deepSNV-method]{estimateDirichlet}}}{\code{signature(control = "deepSNV")}: Method for signature deepSNV}\item{\code{\link[=repeatMask,deepSNV-method]{repeatMask}}}{}\item{\code{\link[=consensusSequence,deepSNV-method]{consensusSequence}}}{}\item{\code{\link[=estimateDispersion,deepSNV,missing-method]{estimateDispersion}}}{}\item{\code{\link[=normalize,deepSNV,missing-method]{normalize}}}{\code{signature(test = "deepSNV", control = "missing")}: Normalize nucleotide counts.}\item{\code{\link[=control,deepSNV-method]{control}}}{\code{signature(deepSNV = "deepSNV")}: Get controls counts.}\item{\code{\link[=test,deepSNV-method]{test}}}{\code{signature(deepSNV = "deepSNV")}: Get test counts.}\item{\code{\link[=p.val,deepSNV-method]{p.val}}}{\code{signature(deepSNV = "deepSNV")}: Get p-values}\item{\code{\link[=coordinates,deepSNV-method]{coordinates}}}{\code{signature(deepSNV = "deepSNV")}: Get coordinates}}}

