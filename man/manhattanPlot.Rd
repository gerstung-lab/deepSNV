\name{manhattanPlot}
\alias{manhattanPlot}
\title{Manhattan plot.}
\usage{
  manhattanPlot(x, col = nt.col)
}
\arguments{
  \item{x}{An \code{\link{deepSNV}} object.}

  \item{col}{An optional vector of colors for the
  nucleotides.}
}
\value{
  NULL.
}
\description{
  This functions performs a Manhattan plot of the p-values
  of a deepSNV test against the position
}
\examples{
data(HIVmix)
manhattanPlot(HIVmix)
}
\author{
  Moritz Gerstung
}

