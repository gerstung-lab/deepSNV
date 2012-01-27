\name{manhattanPlot}
\alias{manhattanPlot}
\title{Manhattan plot.}
\usage{manhattanPlot(x, col=nt.col)
}
\description{Manhattan plot.}
\details{This functions performs a Manhattan plot of the p-values of a deepSNV test against the position}
\value{NULL.}
\author{gemoritz}
\arguments{\item{x}{An \code{\link{deepSNV}} object.}
\item{col}{An optional vector of colors for the nucleotides.}
}
\examples{data(HIVmix)
manhattanPlot(HIVmix)}

