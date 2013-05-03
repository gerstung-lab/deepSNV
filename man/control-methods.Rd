\docType{methods}
\name{control}
\alias{control}
\alias{control,deepSNV-method}
\title{Get control counts}
\arguments{
  \item{deepSNV}{a \code{\link{deepSNV-class}} object}

  \item{total}{Logical. If true the sum of both strands is
  returned}
}
\value{
  A matrix with the absolute frequencies summed over both
  strands.
}
\description{
  Convenience function to obtain the control counts from a
  deepSNV object.
}
\examples{
data(HIVmix)
control(HIVmix)[1:10,]
control(HIVmix, total=TRUE)[1:10,]
}

