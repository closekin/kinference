\name{est_ALF_6way}
\alias{est_ALF_6way}
\title{Estimate ALFs from 6-way genotypes and snerr}
\description{Performs 6-way re-estimation of ALFs, given 6-way genotypes, starting 4-way estimates of ALFs (from est_ALF_ABO_quick), and estimates of "snerr" (single-to-null error rates for apparent homozgyotes, per locus). Used as a second-pass estimate after ALFs have already been estimated roughly from 4-way genotypes (hence the phrase "re-estimation"). Used in CSIRO pipelines, e.g. for SBTuna, where the genotyping and estimation of snerr has already been done by routines in the \pkg{genocalldart} package. If this makes no sense to you, then just stay away from 6-way genotyping!
}
\usage{est_ALF_6way(snpg, control = list())
}
\arguments{
\item{  snpg}{a \code{snpgeno} object with 6-way genotypes (i.e., \code{diplos(snpg)} matches \code{get_genotype_enconding()$genotypes6}), with \code{snerr} and \code{pbonzer} included}
\item{  control}{as per \code{nlminb}}
}
\seealso{\code{\link{est_ALF_ABCO}}, \code{\link{re_est_ALF}}, and \code{\link{est_ALF_ABO_quick}}
}
\examples{
head( bluefin$locinfo$snerr) ## has to exist for 6-way genotypes
bluefin$locinfo$pbonzer <- NULL ## remove pre-existing allele freq
## estimates!
bluefin <- est_ALF_ABO_quick( bluefin)
head( bluefin$locinfo$pbonzer)
bluefin <- est_ALF_6way( bluefin)
head( bluefin$locinfo$pbonzer)
}
\keyword{misc}
