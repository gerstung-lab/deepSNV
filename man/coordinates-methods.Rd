\name{coordinates-methods}
\alias{coordinates-methods}
\alias{coordinates}
\title{Get coordinates...}
\description{Get coordinates}
\details{\code{coordinates-methods}: Convenience function to get the coordinates from a deepSNV object.

}
\value{\code{coordinates-methods}: A \code{\link{data.frame}} with columns "chrom(osome)" and "pos(ition)".

}
\arguments{\item{deepSNV}{a \code{\link{deepSNV-class}} object}
}
\examples{data(HIVmix)
coordinates(HIVmix)[1:10,]}
\alias{coordinates,deepSNV-method}
\alias{coordinates}
\usage{\S4method{coordinates}{deepSNV}(deepSNV)
}

