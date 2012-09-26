\name{Extract}
\alias{[,deepSNV,ANY,ANY-method}
\alias{Extract}
\title{Subsetting for deepSNV objects.}
\arguments{
  \item{x}{A \code{\link{deepSNV-class}} object.}

  \item{i}{Row indeces.}

  \item{j}{Column (nucleotide) indeces.}

  \item{drop}{For matrices and arrays.  If \code{TRUE} the
  result is coerced to the lowest possible dimension (see
  the examples).  This only works for extracting elements,
  not for the replacement.  See \code{\link{drop}} for
  further details.  }
}
\value{
  A \code{\link{deepSNV-class}} object.
}
\description{
  Subsetting for deepSNV objects.
}
\examples{
data(HIVmix)
HIVmix[1:10,]
}
\author{
  Moritz Gerstung
}

