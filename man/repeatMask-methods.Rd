\name{repeatMask-methods}
\alias{repeatMask-methods}
\alias{repeatMask}
\title{Mask homopolymeric repeats.}
\description{Mask homopolymeric repeats.}
\details{\code{repeatMask-methods}: This function masks homopolymeric repeats longer than a given width. These are hot-spots of sequencing error and can confound the analysis.

}
\value{\code{repeatMask-methods}: A boolean vector where TRUE indicates a non-homopolymeric region.

}
\author{gemoritz}
\arguments{\item{x}{An object. Either a \code{\link{deepSNV-class}} object or a \code{\link{DNAString}} with the nucleotide sequence.}
\item{flank}{Boolean. Indicates whether the sites adjacent to the repeat should also be masked.}
\item{w}{Integer. The minimal length at which repeats should be masked. Default \code{w=0}.}
}
\examples{data(HIVmix)
which(repeatMask(HIVmix))}
\alias{repeatMask,DNAString-method}
\alias{repeatMask}
\usage{\S4method{repeatMask}{DNAString}(x, w=5, flank=TRUE)
\S4method{repeatMask}{deepSNV}(x, w=5, flank=TRUE)
}
\alias{repeatMask,deepSNV-method}
\alias{repeatMask}

