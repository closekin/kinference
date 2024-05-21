\name{est_ALF_ABCO}
\alias{est_ALF_ABCO}
\title{Estimate allele frequencies, including nulls}
\description{Performs estimation of ALFs, given a snpgeno or loc.ar object with \code{geno_amb} attribute. Uses "ABCO" genotypes, ie up to 3 scorable alleles plus possible nulls, eg from \code{geno_deambig_ABC}. NALF (null allele frequency) is estimated from HWE deviations, so this requires a decent sample size. But, it doesn't require the elaborate 6-way genotyping and the massive read-depths needed for that. Largely superseded by est_ALF_ABO_quick, except in cases where C alleles are used.
}
\usage{est_ALF_ABCO(lociar, geno_amb = lociar@geno_amb)
}
\arguments{\item{  lociar}{a \code{loc.ar} object with the @geno_amb attribute.}
}
\value{Returns the input, adding a 4-column matrix \code{pambig} to the "locinfo" attribute, plus attributes \code{gobs} and \code{gpred} showing observed and expected counts of each genotype per locus.
}
\seealso{\code{est_ALF_ABO_quick, } re_est_ALF', and \code{ est_ALF_6way}'.
}
\examples{
data( bluefin)
library( atease)
library( mvbutils)
bluefin@locinfo$pbonzer <- NULL ## remove pre-existing ALFs
bluefin@locinfo$snerr <- NULL ## remove pre-existing snerr
## bluefin <- est_ALF_ABCO( bluefin) ## error :(
## geno_amb is a version of the snpg where all single-nulls
## are re-coded as ambiguous single-null / homozygote genotypes
gbasics::define_genotypes()
bluefin@diplos ## check what genotypes exist in the dataset
geno_amb <- bluefin
geno_amb@diplos <- genotypes_ambig ## note genotypes_ambig has
## possibility for C alleles, but bluefin doesn't start with these.
## Turn each genotype into its single-null-ambiguous equivalent.
geno_amb[ bluefin == "AA"] <- "AAO"
geno_amb[ bluefin == "AO"] <- "AAO"
geno_amb[ bluefin == "BB"] <- "BBO"
geno_amb[ bluefin == "BO"] <- "BBO"
geno_amb[ bluefin == "AB"] <- "AB"
geno_amb[ bluefin == "OO"] <- "OO"
bluefin@geno_amb <- geno_amb
bluefin <- est_ALF_ABCO( bluefin)
head( bluefin@locinfo$pbonzer)
}
\keyword{misc}
