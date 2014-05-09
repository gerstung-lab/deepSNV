\docType{methods}
\name{repeatMask}
\alias{repeatMask}
\alias{repeatMask,deepSNV-method}
\alias{repeatMask,DNAString-method}
\title{Mask homopolymeric repeats.}
\arguments{
  \item{x}{An object. Either a \code{\link{deepSNV-class}}
  object or a \code{\link{DNAString}} with the nucleotide
  sequence.}

  \item{flank}{Boolean. Indicates whether the sites
  adjacent to the repeat should also be masked.}

  \item{w}{Integer. The minimal length at which repeats
  should be masked. Default \code{w=0}.}
}
\value{
  A boolean vector where TRUE indicates a non-homopolymeric
  region.
}
\description{
  This function masks homopolymeric repeats longer than a
  given width. These are hot-spots of sequencing error and
  can confound the analysis.
}
\examples{
data(HIVmix)
which(repeatMask(HIVmix))
}
\author{
  Moritz Gerstung
}

