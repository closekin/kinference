\name{check6and4}
\alias{check6and4}
\title{Locus QC check}
\description{Checks 6-way and 4-way genotype frequencies against HWE expectations, and generates plots of observed / expected frequencies. Recently moved into kinference from genocalldart.
}
\usage{
check6and4(
  geno6,
  thresh_pchisq_6and4,
  return_what = c("just_pvals", "all"),
  extra_title = ""
)
}
\arguments{
\item{  geno6}{a \code{snpgeno} object with 4-way and 6-way genocalls}
\item{  thresh_pchisq_6and4}{thresholds for \code{bad} and \code{really bad} p-values}
\item{  return_what}{one of \code{just_pvals} or \code{all}; see value}
\item{  extra_title}{a character string to be added to the bottom-right corner of all plots. Best if < 25 characters.}
}
\value{Creates per-locus vectors \code{pval6} and \code{pval4} for 6-way and 4-way genotypes respectively. If \code{return_what="just_pvals"}, these are returned in a list; if \code{return_what="all"}, they are added as columns to \code{geno6$locinfo}.
}
\keyword{misc}
