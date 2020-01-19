\name{juicyQC}
\alias{juicyQC}
\title{Combine two QC steps: experimental}
\description{\code{juicyQC} is a bivariate combination of \code{\link{hetzminoo_fancy}} and \code{\link{ilglk_geno}}, which \emph{may} be better at identifying "suspect" samples. It may or may not end up in the final kinference toolchain.
}
\usage{juicyQC(snpg, focusees, boring = 1 \%upto\% nrow(snpg))
}
\arguments{
\item{  snpg}{a param}
\item{  focusees}{a param}
\item{  boring}{a param, default \code{1 \%upto\% nrow(snpg)}. See mvbutils for \code{\%upto\%}.}
}
\keyword{misc}
