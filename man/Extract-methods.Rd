\docType{methods}
\name{Extract}
\alias{[,deepSNV,ANY,ANY-method}
\alias{Extract}
\title{Subsetting for deepSNV objects.}
\arguments{
  \item{x}{A \code{\link{deepSNV-class}} object.}

  \item{i}{Row indeces.}

  \item{j}{Column (nucleotide) indeces.}
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

