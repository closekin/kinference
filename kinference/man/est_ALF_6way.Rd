\name{est_ALF_6way}
\alias{est_ALF_6way}
\title{Estimate ALFs from 6-way genotypes and snerr}
\description{Performs 6-way re-estimation of ALFs, given 6-way
  genotypes, starting 4-way estimates of ALFs (from est_ALF_ABO_quick),
  and snerr. Used as a second-pass estimate after ALFs have already been
  calculated based on 4-way genotypes.
}
\usage{est_ALF_6way(snpg, control = list())
}
\arguments{
\item{  snpg}{a \code{snpgeno} object with six-way genotypes (i.e.,
  \code{ snpg@diplos == genotypes6} ), with \code{snerr} and
  \code{pbonzer} included
}
\item{  control}{as per \code{nlminb}}
}
\seealso{
  \code{ \link{ est_ALF_ABCO}}, \code{ \link{ re_est_ALF}}, and \code{
    \link{ est_ALF_ABO_quick}}
}
\examples{
data( bluefin)
head( bluefin$locinfo$snerr) ## has to exist for 6-way genotypes
bluefin@locinfo$pbonzer <- NULL ## remove pre-existing allele freq
## estimates!
bluefin <- est_ALF_ABO_quick( bluefin)
bluefin <- est_ALF_6way( bluefin)
head( bluefin@locinfo$pbonzer)
}
\keyword{misc}
