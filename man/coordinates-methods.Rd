\docType{methods}
\name{coordinates}
\alias{coordinates}
\alias{coordinates,deepSNV-method}
\title{Get coordinates}
\arguments{
  \item{deepSNV}{a \code{\link{deepSNV-class}} object}
}
\value{
  A \code{\link{data.frame}} with columns "chrom(osome)"
  and "pos(ition)".
}
\description{
  Convenience function to get the coordinates from a
  deepSNV object.
}
\examples{
data(HIVmix)
coordinates(HIVmix)[1:10,]
}

