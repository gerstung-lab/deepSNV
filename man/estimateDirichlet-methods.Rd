\name{estimateDirichlet-methods}
\alias{estimateDirichlet-methods}
\alias{estimateDirichlet}
\title{Learn a base-specific Dirichlet prior.}
\description{Learn a base-specific Dirichlet prior.}
\details{\code{estimateDirichlet-methods}: The prior learns the parameters of a Dirichlet distribution seperately for each consensus base. The expected value of the Dirichlet distributions 
is the base-substitution matrix, where rows correspond to the initial nucleotide and columns to the substituted nucleotide. The absolute values determine 
the higher moments of the Dirichlet distributions. After having learned the prior the \code{\link{deepSNV-class}} test is recomputed.

}
\value{\code{estimateDirichlet-methods}: An \code{\link{deepSNV-class}} object.

}
\author{gemoritz}
\arguments{\item{control}{Either a matrix with nucleotide counts or a \code{\link{deepSNV-class}} object.}
}
\examples{data(phiX)
estimateDirichlet(phiX)}
\alias{estimateDirichlet,matrix-method}
\alias{estimateDirichlet}
\usage{\S4method{estimateDirichlet}{matrix}(control)
\S4method{estimateDirichlet}{deepSNV}(control)
}
\alias{estimateDirichlet,deepSNV-method}
\alias{estimateDirichlet}

