\docType{methods}
\name{estimateDispersion}
\alias{estimateDispersion}
\alias{estimateDispersion,deepSNV,missing-method}
\alias{estimateDispersion,matrix,matrix-method}
\title{Estimate the Dispersion factor in a beta-binomial model.}
\arguments{
  \item{test}{Either a deepSNV object, or a matrix with the
  test counts.}

  \item{control}{Missing if test is a deepSNV object,
  otherwise missing.}

  \item{alternative}{The alternative to be tested. One of
  "greater", "less", "two-sided" (default). If test is a
  deepSNV object, automatically taken from the
  corresponding slot if unspecified.}

  \item{interval}{The interval to be screened for the
  overdispersion factor. Default (0,1000).}
}
\value{
  A \code{\link{deepSNV-class}} object if the input was a
  deepSNV object. Otherwise the loglikelihood and the
  estimated parameter.
}
\description{
  This function estimates the dispersion factor in a
  beta-binomial model of the nucleotide counts. This model
  assumes that the count for nucleotide j at position i is
  distributed after a beta-binomial \eqn{X_{i,j}\sim
  \mathrm{BB}(n_i; \alpha, \beta_{ij})}{X_ib ~ BB(n_i;
  alpha, beta_ij)}, where \eqn{n_i}{n_i} is the coverage.
  The base and nucleotide specific parameter
  \eqn{\beta_{ij}}{beta_ij} is estimated from the local
  mean by the method-of-moments estimate,
  \eqn{\alpha}{alpha} is a shared overdispersion parameter.
  It is estimated via a numerical optimization of the
  likelihood under the null-hypothesis.
}
\examples{
data("RCC", package="deepSNV")
plot(RCC)
summary(RCC)[,1:6]
RCC.bb = estimateDispersion(RCC, alternative = "two.sided")
summary(RCC.bb)
}
\author{
  Moritz Gerstung
}

