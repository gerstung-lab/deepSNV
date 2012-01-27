\name{p.combine}
\alias{p.combine}
\title{Combine two p-values...}
\usage{p.combine(p1, p2, method=c("fisher", "max", "average", "prod"))
}
\description{Combine two p-values}
\details{This function combines two P-values into a single one using a statistic defined by method. "fisher" uses the product of the two, in this case the logarithm of the product is
\eqn{\chi^2_4} distributed. If the method = "max", the resulting P-value is \eqn{\max\{P_1,P_2\}^2}. For method = "average" the mean is used, yielding a P-value of \eqn{2 x^2}{2 x^2} if \eqn{x=(P_1+P_2)/2 < .5}{x=(P_1+P_2)/2 < .5}
and  \eqn{1-2 x^2}{1-2 x^2} otherwise.}
\value{p-values}
\author{gemoritz}
\arguments{\item{p1}{P-value 1}
\item{p2}{P-value 2}
\item{method}{One of "fisher" (default), "max" or "average"}
}
\examples{p1 <- runif(1000)
p2 <- runif(1000)
hist(p1)
p.avg = p.combine(p1,p2, method="average")
hist(p.avg)
p.fish = p.combine(p1,p2, method="fisher")
hist(p.fish)
p.max = p.combine(p1,p2, method="max")
hist(p.max)
pairs(data.frame(p1,p2,p.fish,p.max,p.avg))}

