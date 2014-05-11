\docType{data}
\name{RCC}
\alias{RCC}
\title{Example RCC data}
\description{
  Deep sequencing experiments of a renal cell carcinoma and
  healthy control tissue.
}
\examples{
data("RCC", package="deepSNV")
summary(RCC, adjust.method="bonferroni")[,1:6]
plot(RCC)
RCC.bb <- estimateDispersion(RCC, alternative="two.sided")
summary(RCC.bb, adjust.method="bonferroni")[,1:6]
plot(RCC.bb)
}
\keyword{datasets}

