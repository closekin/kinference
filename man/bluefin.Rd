\name{bluefin}
\alias{bluefin}
\alias{data}
\docType{data}
\title{Southern Bluefin Tuna data}
\description{
An anonymised \code{snpgeno} dataset of Southern Bluefin Tuna (SBT) genotypes for 1038 individuals at 1510 loci.

This number of loci would be too small for finding HSPs if we only had biallelic loci with no nulls. However, SBT genotyping is "6-way"; that is, true nulls are present (in fact very common), and the genotyping process attempts to distinguish single-nulls from homozygotes based on read-depth (which is pretty high, typically 50-100 per copy). It does so pretty well but not perfectly, so there is an estimated error rate between single-nulls and homozygotes which is allowed for in the kin-finding calculations. The separation of single-nulls from homozygotes adds considerable statistical power, though also a lot of statistical pain. Calling of genotypes was done by in-house CSIRO software.

Note that this genotyping system is "legacy"; CSIRO is unlikely to use it in future projects, since economics (and statistical simplicity) now favours more loci at lower read-depth, and technological changes have reduced the prevalence of true nulls. However, it is still being used successfully for ongoing SBT kin-finding and CKMR.
}
\usage{bluefin
}
\format{An object of class \code{gbasics::snpgeno()}
}
\keyword{data}
