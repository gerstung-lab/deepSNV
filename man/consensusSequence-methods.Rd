\docType{methods}
\name{consensusSequence}
\alias{consensusSequence}
\alias{consensusSequence,deepSNV-method}
\alias{consensusSequence,matrix-method}
\title{Calculate the consensus sequence.}
\arguments{
  \item{x}{An object. Either an \code{\link{deepSNV-class}}
  object, or a named matrix with nucleotide counts.}

  \item{vector}{Boolean where TRUE indicates that a
  character vector should be returned.}

  \item{haploid}{Should the consensus be called for a
  haploid control? Otherwise, also all bases larger than
  het.cut are rerported. Default haploid = TRUE.}

  \item{het.cut}{Heterozygous cutoff. If haploid = FALSE,
  report all nucleotides with relative frequency larger
  than het.cut. Default = 0.333.}

  \item{...}{Additional arguments passed to methods.}
}
\value{
  A \code{\link{DNAString}} with the consensus sequence, or
  if vector = TRUE, a character vector.
}
\description{
  This function computes the consensus sequence from a
  matrix of nucleotide counts, or the control slot of a
  deepSNV object.
}
\examples{
data(HIVmix)
seq = consensusSequence(HIVmix)
consensusSequence(HIVmix, vector=TRUE)[1:10]
}
\author{
  Moritz Gerstung
}

