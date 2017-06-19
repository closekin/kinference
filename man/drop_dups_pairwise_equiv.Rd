\name{drop_dups_pairwise_equiv}
\alias{drop_dups_pairwise_equiv}
\title{Equivalence class sifter}
\description{
Constructs equivalence classes from pairwise equivalences, and returns the "surplus" elements; if you then drop those elements, only one element from each eq-class will be retained. Requires 2-col matrix showing equivalent pairs.

Code is taken from Numerical Recipes so I should rewrite it perhaps (original algorithm is by Knuth).
}
\usage{drop_dups_pairwise_equiv(ij, want_groups = FALSE)
}
\arguments{
\item{ ij}{2-column matrix or data.frame; probably "row numbers" in a dataset, though might work with character strings too}
\item{ want_groups}{if TRUE, also return the equivalence-classes themselves, as attribute \code{groups}.}
}
\value{Surplus elements in \code{ij}, perhaps plus attributes \code{groups} if \code{want_groups=TRUE}. You can look at that to figure out which elements are being retained (one "representative" from each equiv class).
}
\examples{
pairs <- matrix( c(
294, 289,
328, 294,
904, 857,
905, 904),
    ncol=2, byrow=T)
drop_dups_pairwise_equiv( pairs, T)
#[1] 289 328 857 905
#attr(,"groups")
#attr(,"groups")$`5`
#[1] 294 328 289
#
#attr(,"groups")$`6`
#[1] 904 905 857
#
}
\keyword{
}
