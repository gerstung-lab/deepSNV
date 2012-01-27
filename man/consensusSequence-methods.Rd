\name{consensusSequence-methods}
\alias{consensusSequence-methods}
\alias{consensusSequence}
\title{Calculate the consensus sequence.}
\description{Calculate the consensus sequence.}
\details{\code{consensusSequence-methods}: This function computes the consensus sequence from a matrix of nucleotide counts, or the control slot of a deepSNV object.

}
\value{\code{consensusSequence-methods}: A \code{\link{DNAString}} with the consensus sequence, or if vector = TRUE, a character vector.

}
\author{gemoritz}
\arguments{\item{x}{An object. Either an \code{\link{deepSNV-class}} object, or a named matrix with nucleotide counts.}
\item{vector}{Boolean where TRUE indicates that a character vector should be returned.}
\item{haploid}{Should the consensus be called for a haploid control? Otherwise, also all bases larger than het.cut are rerported. Default haploid = TRUE.}
\item{het.cut}{Heterozygous cutoff. If haploid = FALSE, report all nucleotides with relative frequency larger than het.cut. Default = 0.333.}
}
\examples{data(HIVmix)
seq = consensusSequence(HIVmix)
consensusSequence(HIVmix, vector=TRUE)[1:10]}
\alias{consensusSequence,matrix-method}
\alias{consensusSequence}
\usage{\S4method{consensusSequence}{matrix}(x, vector=FALSE, haploid=TRUE, het.cut=0.333)
\S4method{consensusSequence}{deepSNV}(x, vector=FALSE, haploid=TRUE, het.cut=0.333)
}
\alias{consensusSequence,deepSNV-method}
\alias{consensusSequence}

