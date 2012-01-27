\name{trueSNVs}
\alias{trueSNVs}
\alias{trueSNVs}
\title{Example...}
\description{Example .bam data and true SNVs.}
\details{Two .bam alignments as example data sets are downloaded remotely via http. Sequenced were a 1,512 nt fragment of the HIV genome and a mixture (90\% + 10\%) with another variants. 
The two sequences were confirmed by Sanger sequencing and stored in the table trueSNVs.}
\docType{data}
\alias{HIVmix}
\examples{data(HIVmix)
data(trueSNVs)
table(p.adjust(p.val(HIVmix), method="BH") < 0.05, trueSNVs)}

