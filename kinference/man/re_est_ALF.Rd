\name{re_est_ALF}
\alias{re_est_ALF}
\title{Re-estimate allele frequencies after read-in with load_whopper}
\description{A thin wrapper around est_ALF_ABCO. Performs estimation of
  allele frequencies after creating the required 'geno_amb'
  attribute. Largely superseded by est_ALF_ABO_quick.
}
\usage{re_est_ALF(snpg)
}
\arguments{\item{  snpg}{an snpg object}
}
\value{Returns the input, adding a 4-column matrix \code{pambig} to the
  "locinfo" attribute, plus attributes \code{gobs} and \code{gpred}
  showing observed and expected counts of each genotype per locus.
}
\seealso{\code{\link{est_ALF_ABO_quick}, \code{ \link{ re_est_ALF}}, and \code{ \link{ est_ALF_6way}}}.
}
\examples{
data( bluefin)
bluefin@locinfo$pbonzer <- NULL ## remove pre-existing ALFs
bluefin@locinfo$snerr <- NULL ## remove pre-existing snerr
bluefin <- re_est_ALF( bluefin)
}
\keyword{misc}
