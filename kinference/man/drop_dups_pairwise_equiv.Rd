\name{drop_dups_pairwise_equiv}
\alias{drop_dups_pairwise_equiv}
\title{Grouping duplicate samples}
\description{
\code{\link{find_duplicates}} only does pairwise comparisons. However, tissue from the same animal may turn up in multiple samples, so that one sample may turn up in many duplicate-pairs, and the pairs are linked. This function constructs equivalence classes--- each corresponding notionally to one \emph{animal}--- showing which samples belong in each class. It can either return the entire set of classes, or it can pick just one sample from each class and then return the "surplus" duplicate samples; if you then drop those elements, only one sample from each animal will be retained.

MVB adds: it is some years since I checked this code, but I think it keeps merging classes whenever a sample in one class is flagged as a duplicate of a sample in a different class. Thus, it will be sensitive to false-positive duplicates (though less so to false-negative ones). It's up to you to make sure that the input really contains true duplicates!
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
\seealso{chain.pairwise
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
