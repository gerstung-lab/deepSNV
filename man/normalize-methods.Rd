\name{normalize}
\alias{normalize}
\alias{normalize-methods}
\alias{normalize,deepSNV,missing-method}
\alias{normalize,matrix,matrix-method}
\title{Normalize nucleotide counts.}
\usage{
  normalize(test, control, ...)

  \S4method{normalize}{matrix,matrix}(test, control,
    round=TRUE, ...)

  \S4method{normalize}{deepSNV,missing}(test, control, ...)
}
\arguments{
  \item{test}{Either an \code{\link{deepSNV-class}} object
  or a named matrix with nucleotide counts.}

  \item{control}{Missing if test is an
  \code{link{deepSNV-class}} object, otherwise a matrix
  with nucleotide counts.}

  \item{round}{Logical. Should normalized counts be rounded
  to integers? Default=TRUE}

  \item{...}{Parameters passed to \code{\link{loess}}.}
}
\value{
  A \code{\link{deepSNV-class}} object.
}
\description{
  This functions performs a \code{\link{loess}}
  normalization of the nucleotide. This experimental
  feature can be used to compare experiments from different
  libraries or sequencing runs that may have differing
  noise characteristics.

  Normalize nucleotide counts.

  Normalize nucleotide counts.
}
\note{
  This feature is somewhat experimental and the results
  should be treated with care. Sometimes it can be better
  to leave the data unnormalized and use a model with
  greater dispersion instead.
}
\examples{
data(phiX, package = "deepSNV")
plot(phiX)
phiN <- normalize(phiX, round = TRUE)
plot(phiN)
}
\author{
  Moritz Gerstung
}

