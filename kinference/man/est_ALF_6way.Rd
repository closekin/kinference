\name{est_ALF_6way}
\alias{est_ALF_6way}
\title{Estimate ALFs from 6-way genotypes and snerr}
\description{Performs "straight" estimation of ALFs, given 6-way genotypes and snerr. Won't allow for changes in C-allele frequency from one population to the next. In principle, should just use \code{genocalldart::choose_geno6_thresholds} but fix the count-related thresholds and re-estimate ALFs. (What??)
}
\usage{est_ALF_6way(snpg, control = list())
}
\arguments{
\item{  snpg}{a \code{snpgeno} object with \code{snerr} included}
\item{  control}{as per \code{nlminb}}
}
\keyword{misc}
