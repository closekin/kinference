\name{est_ALF_ABO_quick}
\alias{est_ALF_ABO_quick}
\alias{est_ALF_nonulls}
\title{Estimate allele frequencies with nulls}
\description{
\code{est_ALF_ABO_quick} is the recommended way to estimate allele frequencies from called genotypes, \emph{provided that} you are using biallelic SNPs wtih 4-way genotyping and you believe that null alleles are a real, repeatable thing in your data. If you have biallelic SNPs but you don't believe you have repeatable, heritable null alleles, then use \code{est_ALF_nonulls} instead. If you are using genuine 6-way genotyping with nulls, see \bold{Beyond 4 way genotypes}.

So... \code{est_ALF_ABO_quick} is for fast maximum-likelihood estimation of A (nominally major), B (nominally minor), and O (bona fide null) allele frequencies for a set of loci. It uses the EM algorithm plus Aitken acceleration; this means the whole calculation can be vectorized across loci, which (compared to direct maximization of the log-likelihood) more than compensates for the notorious inefficiency of EM (and also Aitken helps \emph{a lot}). See the MS or the code for more details. Missing data not allowed.

\code{est_ALF_nonulls} is very simple; it just looks at the ratio of A to B alleles across all samples at each locus, and sets the null frequency to 0. This \emph{can} cope with missing data, and any sample-loci recorded as "OO" (double null) is treated as such. Thus, you could do \code{est_ALF_nonulls} followed by a random-imputation step to fill in the missings, if you really can't unmissingize them (eg by using a more confident genotype-calling algorithm to your raw data). See \bold{Examples}.
}
\usage{
est_ALF_ABO_quick(x = NULL, AB, AAO, BBO, OO,
    tol = 0.0000001, EMtol = 0.001,  quietly = FALSE,
    MAX_AITKEN= 40, return_unconverged= FALSE)
est_ALF_nonulls( x=NULL, AB, AA, BB,
  pbonzer_format= FALSE)
}
\arguments{
\item{  x}{a \code{snpgeno} object, or NULL to use the next 4 or 3 args explicitly. \code{diplos(x)} should either be \code{genotypes6} or \code{genotypes4_ambig}. For \code{est_ALF_nonulls}, the latter is treated as if AAO always truly means AA, and BBO means BB; OO is treated as missing.}
\item{  AB, AAO, BBO, OO, AA, BB}{vectors (per locus) of counts of these genotypes. Can't mix with non-null \code{x}.}
\item{  tol}{final convergence tolerance (in Aitken steps)}
\item{  EMtol}{tolerance within the EM steps; after this is achieved, try an Aitken step}
\item{  quietly}{if TRUE, then at the end print information on the number of iterations required}
\item{  MAX_AITKEN}{Maximum number of Aitken-accelerations to allow. Most loci converge within 10. We have seen a few where over 30 Aitkens are needed; however, the results after just 10 were still pretty good. You can try pushing it higher, but do check your weird loci for extreme weirdness (see next argument).}
\item{  return_unconverged}{if TRUE and some loci still haven't converged after \code{MAX_AITKEN} iterations, then just return the indices of those loci. By default (ie if this is FALSE), you'll just get the estimates with a warning, and that's probably fine.}
\item{  pbonzer_format}{if \code{est_ALF_nonulls} is called directly on \code{AB} etc rather than on \code{x}, then \code{pbonzer_format} determines the format of the returned allele frequency estimates:TRUE means you get a 4-column matrix suitable for \code{kinference} (to go in \code{$locinfo$pbonzer}) and FALSE means you just get a vector of the major (A) allele frequencies.}
}
\section{Beyond 4 way genotypes}{
\code{est_ALF_ABO_quick} does actually accept 6-way genotypes, but it would be a bit weird to use it, because the first thing it then does is to re-encode the genotypes as 4-way, thus sacrificing some statistical information. It starts by merging the single-nulls with true homozygotes, and then tries to "unmerge" them statistically! It would be better to use \code{\link{est_ALF_6way}}, though you might need to first run \code{\link{est_ALF_ABCO}} or \code{est_ALF_ABO_quick} to get starting values.

\code{\link{est_ALF_ABCO}} (qv) is a much slower version that can handle triallelic SNPs (ie potentially with a C allele). Like \code{est_ALF_ABO}, it does not distinguish between single-nulls and homozygotes. It should give the same results as \code{est_ALF_ABO} for loci without a C allele. Syntax is a bit different because of its role in the legacy 6-way pipeline, so you might prefer to use \code{\link{re_est_ALF}} (qv) which does the same thing but with different (easier?) formatting of the input.

\code{\link{est_ALF_6way}} (qv) uses 6-way genotypes, where single nulls are called separately from homozygotes, but potentially with error. It requires error estimates in \code{$locinfo$snerr}. Unless you are dealing with a legacy 6-way dataset, you don't want to go near this!
}
\value{If \code{x} is supplied, then its \code{locinfo} attribute will be augmented with the \code{pbonzer} (allele frequency) matrix required by most \code{\link{kinference}} functions. Note that \code{pbonzer} has 4 columns always, so here the 3rd column ("C") is set to zero. If \code{x} is not supplied, then a 3-column matrix matrix is returned. Rowsums of the matrix are always 1 in either case.
}
\examples{
#dropbears$locinfo$pbonzer <- NULL ## no population allele frequency estimates!
#dropbears <- est_ALF_ABO_quick( dropbears)
#head( dropbears$locinfo$pbonzer)
#dropbears <- est_ALF_nonulls( dropbears)
#head( dropbears$locinfo$pbonzer)
## Randomly make some values missing...
#nb <- nrow( dropbears)
#nl <- ncol( dropbears)
#nmiss <- round( 0.1 * nb * nl)
#missij <- cbind( rsample( nmiss, 1:nb, replace=TRUE),
#    rsample( nmiss, 1:nl, replace=TRUE))
#missbears <- dropbears
#missbears[ missij] <- 'OO'
#missbears <- est_ALF_nonulls( missbears)
#plot( missbears$locinfo$pbonzer[,1], dropbears$locinfo$pbonzer[,1]) # very similar
#abline( 0, 1)
#imputor <- function( x)\{
#  # I don't *recommend* this, but you *could* randomly impute "missing" (OO) genotypes like this
#  # if you are sure there's no nulls
#  # I haven't tested the code!
#  # Caveat emptor... and read it carefully to work out what it's
#  # (hopefully)
#  # doing...
#  pA <- x$locinfo$pbonzer[,1]
#  pB <- 1-pA
#  pAB <- 2*pA*pB
#  pAA <- sqr( pA)
#  pBB <- sqr( pB)
#  misso <- which( x=='OO', arr.ind=TRUE)
#  r <- runif( nrow( misso))
#  x[ misso] <- 'BB' # default
#  x[ misso[ r < (pAB+pAA)[ misso[,2]],] <- 'AA'
#  x[ misso[ r < pAB[ misso[,2]],] <- 'AB'
#return( x)
#\}
}
\keyword{misc}
