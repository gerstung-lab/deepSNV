\docType{methods}
\name{p.val}
\alias{p.val}
\alias{p.val,deepSNV-method}
\title{Get p-values}
\arguments{
  \item{deepSNV}{a \code{\link{deepSNV-class}} object}
}
\value{
  A matrix with the p-values.
}
\description{
  Convenience function to get the p-values from a deepSNV
  object.
}
\examples{
data(HIVmix)
p.val(HIVmix)[1:10,]
}

