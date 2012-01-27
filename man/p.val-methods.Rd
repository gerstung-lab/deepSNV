\name{p.val-methods}
\alias{p.val-methods}
\alias{p.val}
\title{Get p-values...}
\description{Get p-values}
\details{\code{p.val-methods}: Convenience function to get the p-values from a deepSNV object.

}
\value{\code{p.val-methods}: A matrix with the p-values.

}
\arguments{\item{deepSNV}{a \code{\link{deepSNV-class}} object}
}
\examples{data(HIVmix)
p.val(HIVmix)[1:10,]}
\alias{p.val,deepSNV-method}
\alias{p.val}
\usage{\S4method{p.val}{deepSNV}(deepSNV)
}

