\name{pbetabinom}
\alias{pbetabinom}
\title{Cumulative beta-binomial probability distribution}
\usage{
  pbetabinom(x, size, prob, dispersion, log = FALSE)
}
\arguments{
  \item{x}{Counts}

  \item{size}{Sample size}

  \item{prob}{Probability}

  \item{dispersion}{Dispersion rho/(1-rho)}

  \item{log}{Should the natural logarithm be returned?
  Default=FALSE}
}
\value{
  Probability
}
\description{
  The underlying implementation is written in C
}
\author{
  mg14
}

