\name{test-methods}
\alias{test-methods}
\alias{test}
\title{Get test counts...}
\description{Get test counts}
\details{\code{test-methods}: Convenience function to obtain the test counts from a deepSNV object.

}
\value{\code{test-methods}: A matrix with the absolute frequencies summed over both strands.

}
\arguments{\item{deepSNV}{a \code{\link{deepSNV-class}} object}
\item{total}{Logical. If true the sum of both strands is returned}
}
\examples{data(HIVmix)
test(HIVmix)[1:10,]
test(HIVmix, total=TRUE)[1:10,]}
\alias{test,deepSNV-method}
\alias{test}
\usage{\S4method{test}{deepSNV}(deepSNV, total=F)
}

