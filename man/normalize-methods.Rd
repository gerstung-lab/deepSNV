\name{normalize-methods}
\alias{normalize-methods}
\alias{normalize}
\title{Normalize nucleotide counts.}
\description{Normalize nucleotide counts.}
\details{\code{normalize-methods}: This functions performs a \code{\link{loess}} normalization of the nucleotide. This experimental feature can 
be used to compare experiments from different libraries or sequencing runs that may have differing noise characteristics.

}
\value{\code{normalize-methods}: A \code{\link{deepSNV-class}} object.

}
\author{gemoritz}
\arguments{\item{test}{Either an \code{\link{deepSNV-class}} object or a named matrix with nucleotide counts.}
\item{control}{Missing if test is an \code{link{deepSNV-class}} object, otherwise a matrix with nucleotide counts.}
\item{round}{Logical. Should normalized counts be round to integers? Default=TRUE}
\item{...}{Parameters passed to \code{\link{loess}}.}
}
\examples{data(phiX, package = "deepSNV")
plot(phiX)
phiN <- normalize(phiX, round = TRUE)
plot(phiN)}
\alias{normalize,matrix,matrix-method}
\alias{normalize}
\usage{\S4method{normalize}{matrix,matrix}(test, control, round=TRUE, ...)
\S4method{normalize}{deepSNV,missing}(test, control, ...)
}
\alias{normalize,deepSNV,missing-method}
\alias{normalize}

