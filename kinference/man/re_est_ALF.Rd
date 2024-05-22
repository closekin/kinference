\name{re_est_ALF}
\alias{re_est_ALF}
\title{Re-estimate allele frequencies after read-in with load_whopper}
\description{A thin wrapper around est_ALF_ABCO. Performs estimation of allele frequencies after creating the required \code{geno_amb} attribute. Largely superseded by est_ALF_ABO_quick.
}
\usage{re_est_ALF(snpg)
}
\arguments{\item{  snpg}{an snpg object}
}
\value{Returns the input, adding a 4-column matrix \code{pambig} to the "locinfo" attribute, plus attributes \code{gobs} and \code{gpred} showing observed and expected counts of each genotype per locus.
}
\seealso{\code{est_ALF_ABO_quick, } re_est_ALF', and \code{ est_ALF_6way}'.
}
\examples{
data( bluefin)
library( atease)
bluefin <- bluefin[,1:100] ## just 100 loci, otherwise too slow for R CMD check
bluefin@locinfo$pbonzer <- NULL ## remove pre-existing ALFs
bluefin@locinfo$snerr <- NULL ## remove pre-existing snerr
bluefin <- re_est_ALF( bluefin)
}
\keyword{misc}
