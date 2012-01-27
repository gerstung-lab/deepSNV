\name{Extract-methods}
\alias{[,deepSNV-method}
\alias{[}
\title{Subsetting for deepSNV objects.}
\usage{\S4method{[}{deepSNV}(x, i, j)
}
\description{Subsetting for deepSNV objects.}
\value{A \code{\link{deepSNV-class}} object.}
\author{gemoritz}
\alias{[,deepSNV,ANY,ANY-method}
\arguments{\item{x}{A \code{\link{deepSNV-class}} object.}
\item{i}{Row indeces.}
\item{j}{Column (nucleotide) indeces.}
}
\examples{data(HIVmix)
HIVmix[1:10,]}

