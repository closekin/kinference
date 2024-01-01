\name{drop_dups_pairwise_equiv}
\alias{drop_dups_pairwise_equiv}
\title{Grouping pairwise duplicates}
\description{Constructs equivalence classes to show from pairwise equivalences, and returns the "surplus" elements; if you then drop those elements, only one element from each eq-class will be retained.
}
\usage{drop_dups_pairwise_equiv(ij, want_groups = FALSE)
}
\arguments{
\item{  ij}{2-column matrix or data.frame; possibly row numbers in a dataset, or strings (now that \code{find_HSPs} etc can optionally return "row ID" strings)}
\item{  want_groups}{if \code{TRUE}, also return the equivalence-classes themselves, as attribute \code{groups}.}
}
\details{Input should be row numbers in a \code{snpgeno} objects of duplicates, as a two-column data.frame or matrix with each row being a pair of duplicates, or the output from \code{\link{find_duplicates}} (a 3-col matrix). Identifies \code{groups} of equivalent observations (e.g., if i and j are duplicates, and j and k are duplicates, then i, j, and k are all equivalent). Outputs a vector of the row numbers for all-but-one of each group.
}
\value{Surplus elements in \code{ij}, perhaps plus attributes \code{groups} if \code{want_groups=TRUE}. You can look at that to figure out which elements are being retained (one "representative" from each equiv class). If \code{ij} has no rows, an empty integer vector is returned, without any attributes.
}
\examples{
pairs <- matrix( c(
294, 289,
328, 294,
904, 857,
905, 904),
    ncol=2, byrow=TRUE)
drop_dups_pairwise_equiv( pairs, TRUE)
#[1] 289 328 857 905
#attr(,"groups")
#attr(,"groups")$`5`
#[1] 294 328 289
#
#attr(,"groups")$`6`
#[1] 904 905 857
}
\keyword{misc}
