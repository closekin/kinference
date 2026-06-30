\name{chain_pairwise}
\alias{chain_pairwise}
\alias{get_chain}
\title{Sib-groups within HSPs}
\description{
For checking veracity of \emph{potential} half-sibs or other kin-pairs. \code{chain_pairwise} organizes them into chains within which each sample can be linked to another by a succession of direct pairwise links. The general idea is that real HSPs will be in clusters of 2 or 3; a spurious sample with a "lucky" genotype that wants to be everybody's mate will appear in a big but incomplete chain of mostly false-positives, where the direct pairwise links between the other chain-members are weak. You'd only run \code{chain_pairwise} for pairs with a PLOD (or whatever statistic is being used) within a particular suspect range, so each chain may have false-negatives (i.e. missing direct links), but the general idea should be clear.

\code{chain_pairwise} and \code{get_chain} are diagnostic tools for when \code{\link{kinference}} goes wrong. For a comprehensive approach to complicated links between \emph{genuine} sib-groups (full and half), e.g. within larval samples such as for Atlantic Bluefin Tuna, see the \pkg{sibgrouper} package (if not on R-universe, then contact the \code{\link{kinference}} authors).
}
\usage{
chain_pairwise( thing)
get_chain( thing, seed)
}
\arguments{
\item{  thing}{output from \code{find_HSPs} or \code{find_POPs} etc, or some subset thereof}
\item{  seed}{one sample ID, interpreted as a row-number in \code{thing}. To do:also allow names, via \code{info} attr.}
}
\details{\code{get_chain} finds the chain for one specific sample.
}
\value{\code{chain_pairwise} returns a list of matrices, each for one chain; the rows and columns of each matrix are the samples in that chain. A "+" in the matrix indicates that those two samples have a direct pairwise link (i.e., they appear together in one row of \code{thing}); a "." means not. The rows and columns of each matrix are sorted so that the linkiest samples are on the bottom and right. \code{get_chain} returns the row-subset of \code{thing} that is chained to \code{seed}.
}
\keyword{misc}
