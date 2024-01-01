\name{est_ALF_ABO_quick}
\alias{est_ALF_ABO_quick}
\title{Estimate allele frequencies with nulls; much faster than 'est_ALF_ABCO' (qv).}
\description{Fast estimation of A (major), B (minor), and O (null) frequencies for a set of loci. Uses the EM algorithm plus Aitken acceleration; this means the whole calculation can be vectorized across loci, which more than compensates for any EM inefficiency (and Aitken helps \emph{a lot}).
}
\usage{
est_ALF_ABO_quick(x = NULL, AB, AAO, BBO, OO,
    tol = 0.0000001, EMtol = 0.001,  quietly = FALSE)
}
\arguments{
\item{  x}{a \code{snpgeno} object, or NULL to use the next 4 args explicitly.}
\item{  AB, AAO, BBO, OO}{vectors (over all loci) of counts of these genotypes. Can't mix with non-null \code{x}.}
\item{  tol}{final convergence tolerance (in Aitken steps)}
\item{  EMtol}{tolerance within the EM steps; after this is achieved, try an Aitken step}
\item{  quietly}{if TRUE, then at the end print information on the number of iterations required}
}
\details{I should add some more...
}
\value{If \code{x} is supplied, then its \code{locinfo} attribute will be augmented with the \code{pbonzer} (allele frequency) matrix required by most functions in the \pkg{kinference} package. Note that \code{pbonzer} has 4 columns always, so here the 3rd column ("C") is set to zero. If \code{x} is not supplied, then a 3-column matrix matrix is returned. Rowsums of the matrix are always 1 in either case.
}
\seealso{\code{\link{est_ALF_ABCO}}
}
\examples{
# Should have one
}
\keyword{misc}
