\name{control-methods}
\alias{control-methods}
\alias{control}
\title{Get control counts...}
\description{Get control counts}
\details{\code{control-methods}: Convenience function to obtain the control counts from a deepSNV object.

}
\value{\code{control-methods}: A matrix with the absolute frequencies summed over both strands.

}
\arguments{\item{deepSNV}{a \code{\link{deepSNV-class}} object}
\item{total}{Logical. If true the sum of both strands is returned}
}
\examples{data(HIVmix)
control(HIVmix)[1:10,]
control(HIVmix, total=TRUE)[1:10,]}
\alias{control,deepSNV-method}
\alias{control}
\usage{\S4method{control}{deepSNV}(deepSNV, total=F)
}

