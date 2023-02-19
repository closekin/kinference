\name{est_ALF_ABCO}
\alias{est_ALF_ABCO}
\title{Estimate allele frequencies, including nulls}
\description{Uses "ABCO" genotypes, ie up to 3 scorable alleles plus possible nulls, eg from \code{geno_deambig_ABC}. NALF (null allele frequency) is estimated from HWE deviations, so this requires a decent sample size. But, it doesn't require the elaborate 6-way genotyping and the massive read-depths needed for that.
}
\usage{est_ALF_ABCO(lociar, geno_amb = lociar@geno_amb)
}
\arguments{\item{  lociar}{a \code{loc.ar} object with the @geno_amb attribute.}
}
\value{Returns the input, adding a 4-column matrix \code{pambig} to the "locinfo" attribute, plus attributes \code{gobs} and \code{gpred} showing observed and expected counts of each genotype per locus.
}
\seealso{geno_dembig_ABC
}
\keyword{misc}
