\name{makePrior}
\alias{makePrior}
\title{Compute a prior from a COSMIC VCF object}
\usage{
  makePrior(COSMIC, regions, pi.mut = 0.1, pi0 = 1e-04)
}
\arguments{
  \item{COSMIC}{A VCF object from COSMIC VCF export.}

  \item{regions}{A GRanges object with the regions (gene)
  of interest.}

  \item{pi.mut}{Probability that a gene is mutated}

  \item{pi0}{Background probability of a locus being
  mutated. Default 1e-4, corresponding to an expected value
  of 1 SNV per 1e4 bases.}
}
\value{
  A vector of prior values with length given by the length
  of the regions GRanges object.
}
\description{
  This function computes the prior probability of detecting
  a true variant from a variation data base. It assumes a
  VCF file with a CNT slot for the count of a given base
  substitution. Such a VCF file can be downloaded at
  ftp://ngs.sanger.ac.uk/production/cosmic/CosmicCodingMuts_v64_02042013_noLimit.vcf.gz.
  The prior probability is simply defined as pi.mut *
  CNT[i]/sum(CNT). On sites with no count, a background
  probability of pi0 is used.
}
\note{
  Experimental code, subject to changes
}
\examples{
## Make prior (not run)
#COSMIC <- readVcf("PATHTO/CosmicCodingMuts_v64_02042013_noLimit.vcf.gz", genome="GChr37")
#prior <- makePrior(COSMIC[info(COSMIC)$GENE=="TP53"], regions=GRanges(17, IRanges(7571720,7578811)))
#plot(prior[,1], type="h")
}
\author{
  mg14
}

