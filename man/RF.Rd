\name{RF}
\alias{RF}
\title{Relative frequencies.}
\usage{RF(freq, total=FALSE)
}
\description{Relative frequencies.}
\details{Convenience function to compute the relative frequencies from a matrix with absolute counts.}
\value{A matrix with the relative frequencies.}
\author{gemoritz}
\arguments{\item{freq}{A matrix with nucleotide counts.}
\item{total}{If the nucleotide counts have columns for forward and reverse direction, return each strand sepratatelu (FALSE), or add the two (TRUE).}
}
\examples{data(HIVmix)
RF(test(HIVmix))[1:10,]
RF(test(HIVmix), total=TRUE)[1:10,]}

