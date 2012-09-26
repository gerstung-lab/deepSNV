\name{test}
\alias{test}
\alias{test-methods}
\alias{test,deepSNV-method}
\title{Get test counts}
\usage{
  \S4method{test}{deepSNV}(deepSNV, total = FALSE)
}
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
  Convenience function to obtain the test counts from a
  deepSNV object.
}
\examples{
data(HIVmix)
test(HIVmix)[1:10,]
test(HIVmix, total=TRUE)[1:10,]
}

