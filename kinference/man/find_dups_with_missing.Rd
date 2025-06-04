\name{find_dups_with_missing}
\alias{find_dups_with_missing}
\title{Duplicate-finding with some missing genotypes}
\description{
This function is not for CKMR datasets (where missingness is not allowed) but rather for "Gene-Tagging" (individual mark-recapture using genotypes as tags). For GT, the genotyping method has to be cheap rather than high-quality, so that some genotypes just end up missing (ie the genotyping pipeline has decided they are unscorable). Rather than "imputing" those missing genotypes (which you'd have to do for finding kin-pairs), for the specific case of finding duplicate samples it is probably better to just compare across the loci that \emph{are} both called in each pair of samples.

\code{find_dups_with_missing} looks at all pairwise comparisons and "accepts" any where there are not too many inconsistencies (which can arise within true duplicates, from allelic dropout when read-depth is low). This means that a very-poor-quality-DNA sample might actually match with lots of others, simply because it has so few genotypes called! So, subset your data beforehand to remove Bad Eggs.
}
\usage{
find_dups_with_missing(
  snpg,
  subset1 = 1 \%upto\% nrow(snpg),
  subset2 = subset1,
  max_diff_ppn,
  limit = 10000,
  ij_numeric= is.null( rowid_field( snpg))
)
}
\arguments{
\item{  snpg}{a \code{snpgeno} object}
\item{  subset1, subset2}{numeric vectors of which samples to use (not logical, not negative). Defaults to all of them. Iff \code{subset1} and \code{subset2} are identical, only half the comparisons are done (i.e., not \emph{i} with \emph{j} then \emph{j} with \emph{i}). Some sanity checks are made.}
\item{  max_diff_ppn}{What \emph{proportion} of non-missing (ie scored-in-both) loci to treat as the threshold for duplicity?}
\item{  limit}{if you hit this many "duplicates", it will stop, to avoid blowing out memory. It means you set \code{max_diff_ppn} too high. For consistency, we should probably have called this \code{keep_n} as per other \code{find_XXX} functions.}
\item{  ij_numeric}{if FALSE, use the rowid field (see \code{with_rowid_field}) to label the pair-members, rather than their row numbers.}
}
\value{A dataframe with columns \code{ppn_diff}, then \code{i} and \code{j} which show pairs of samples that are within \code{max_diff_ppn}. (Note that \code{i} and \code{j} refer to rows in \code{snpg} itself, not to whatever was passed in the subsets. This is what any normal person would expect, but long ago the software worked differently...)
}
\keyword{misc}
