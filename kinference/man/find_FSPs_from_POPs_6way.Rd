\name{find_FSPs_from_POPs_6way}
\alias{find_FSPs_from_POPs_6way}
\title{Split FSPs from POPS; experimental, obsolete, and private!}
\description{
This version is not finished and probably abandoned; it will be superceded by newer continuum-likelihood method. The reason for its existence, is that (to my surprise) the 4way-based optimal-weight splitter for FSP/POP doesn't work particularly well.

Tries to split Full-Sibling Pairs from Parent-Offspring Pairs, using 6-way genotypes. Takes advantage of the fact that POP-detection will also pick up FSPs (as both relationship classes share at least one allele at every locus), so FSP-detection can be run only within the subset of animals identified as potential POPs, hence low computational demand.
}
\usage{find_FSPs_from_POPs_6way(snpg, candiPOPs, SDwt_POP = 0.5)
}
\arguments{
\item{  snpg}{a \code{snpgeno} object}
\item{  candiPOPs}{a 2-column matrix of of row-numbers in \code{snpg}, for pairs known to be either POPs or FSPs (eg from \code{find_POPs}).}
\item{  SDwt_POP}{scalar for determining "prior weighting" of POPs and FSPs. 0.5 (the default) is probably OK.}
}
\seealso{simcheck_FSP_POP
}
\keyword{misc}
