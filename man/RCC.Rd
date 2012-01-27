\name{RCC}
\alias{RCC}
\alias{RCC}
\title{Example RCC data...}
\description{Example RCC data}
\details{Deep sequencing experiments of a renal cell carcinoma and healthy control tissue.}
\docType{data}
\examples{data("RCC", package="deepSNV")
summary(RCC, adjust.method="bonferroni")[,1:6]
plot(RCC)
RCC.bb <- estimateDispersion(RCC, alternative="two.sided")
summary(RCC.bb, adjust.method="bonferroni")[,1:6]
plot(RCC.bb)}

