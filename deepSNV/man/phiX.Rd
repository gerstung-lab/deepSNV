\docType{data}
\name{phiX}
\alias{phiX}
\title{Example phiX data}
\description{
  Data from two phiX experiments sequenced on a GAIIx.
}
\examples{
data(phiX, package="deepSNV")
plot(phiX)
phiN <- normalize(phiX, round=TRUE)
plot(phiN)
}
\keyword{datasets}

