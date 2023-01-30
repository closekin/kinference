\name{split_FSPs_from_HSPs}
\alias{split_FSPs_from_HSPs}
\alias{split_FSPs_from_POPs}
\alias{split_FSPs_from_candiHSPs}
\title{Find full-sib pairs among parent-offspring pairs or among half-sib pairs}
\description{For pairs already picked as likely parent-offspring pairs (POPs), i.e., those with a weighted pseudo-exclusion (WPSEX) statistic less than some threshold, they might be full sibling pairs (FSPs). This function checks potential POPs with very low WPSEX values for their potential to be FSPs.
}
\usage{
split_FSPs_from_HSPs( snpg, candiHSPs)
split_FSPs_from_POPs( snpg, candiPOPs, gerr=NULL, use_obsolete_version=TRUE)
}
\arguments{
\item{  snpg}{a \code{snpgeno} object}
\item{  candiHSPs, candiPOPs}{candidate kin-pairs- normally, a dataframe with rows being pairs and columns \emph{i} and \emph{j} (and possibly others) e.g. from find_POPs() or find_HSPs(). Can also be a 2-column matrix (each row again one pair).}
\item{  gerr}{genotyping error rate:0.01 would mean 1\%. There's no default so it's up to you; e.g., looking at replicate samples is a good way to estimate it. Results of \code{split_FSP_from_POPs} should not be very sensitive to the value, however.}
\item{  use_obsolete_version}{the original code for POPs-vs-FSPs was based on a different non-likelihood-based statistic; see \bold{Obsolete note}. It turned out to have low statistical power, but did not require specifying a \code{gerr}. For replicability purposes, you can still run it by setting this parameter to \code{TRUE}.}
}
\details{
The idea of \code{split_FSPs_from_POPs}--- though this is not the only possible workflow--- is that pairs which are \emph{either} POPs \emph{or} FSPs should stand out very clearly from everything else, via \code{\link{find_POPs}}. Then the job is to pick between those possibilities. The workflow is supposed to be:

\itemize{
\item nail POPs/FSPs first with \code{\link{find_POPs}}
\item pick between them with \code{split_FSPs_from_POPs} (update: this doesn't work very well yet... use age info if you can)
\item look for HSPs and filter out already-known POPs and FSPs
}

Hence the other function, \code{split_FSPs_from_HSPs}, is theoretically unnecessary in that you have already run \code{\link{find_POPs}} and \code{split_FSPs_from_POPs} so you should know which of your "HSPs" are really something else. But, nevertheless it's handy to have.

However, an equally reasonable workflow might be:

\itemize{
\item nail HSPs and everything stronger with \code{\link{find_HSPs}}
\item split HSPs from POPs/FSPs with \code{split_FSPs_from_HSPs}
\item \code{split_FSPs_from_POPs}
}

All \code{split_} functions return expected values under different possible kin-types (not variances, since these cannot be predicted for all kin-types).
\subsection{Obsolete note}{
The statistic for \code{split_FSPs_from_POPs}- _which isn't as powerful as I'd hoped; I'm going to redo it- is based on the weighted sum of the number of exactly-matching 4-way genotypes, with weights chosen to have high power for this particular discrimination. Weighting is optimized for the unlikely scenario that POPs and FSPs are equally likely a priori, but in practice the weights are not sensitive to this. The test is deliberately crude and robust- e.g. it avoids exclusion-based checks- on the assumption that you have enough loci to pick HSPs, so the more-related kin-types should be slam-dunks. \bold{But} it doesn't seem powerful enough. More worked needed...

\code{split_FSPs_from_HSPs} again uses 4-way genotypes only (to avoid having to worry about errors) but in a properly optimal PLOD designed for FSP/HSP discrimination- its expectation is positive for FSPs and negative for HSPs. Theoretical means for those are returned as attributes (variances cannot be predicted). Haven't added means for POPs or UPs since you're not "supposed" to have those in the mix by the time you run \code{split_FSPs_from_HSPs}, but maybe I should fix that at some point.
}
}
\examples{
# pops_or_fsps <- find_POPs( mysnpg, ...)
## do histograms etc to find likely ones
# discro <- split_FSPs_from_POPs( mysngp, pops_or_fsps \%where\% (wpsex < 0.042), gerr=0.01)
# hist( discro, nc=30, col='grey')
# abline( v=discro@E_POP, col='red')
# text( discro@E_POP, par( 'usr')[4], 'POP', col='red', pos=1) # below
# abline( v=discro@E_FSP, col='lightblue')
# text( discro@E_FSP, par( 'usr')[4], 'FSP', col='lightblue', pos=1) # below
# abline( v=discro@E_HSP, col='pink', lty='dashed')
# text( discro@E_HSP, par( 'usr')[4], 'HSP', col='pink', pos=1) # below
###
# h_or_f <- find_HSPs( mysnpg, ...)
## do histograms etc to find likely sib pairs
# discro2 <- split_FSPs_from_HSPs(  mysngp, h_or_f \%where\% (PLOD > 55))
# hist( discro2$PLOD_FH, nc=20, col='grey') # HSPs "should" be < 0, FSPs > 0
# abline( v=discro2@E_HSP, col='orange')
# text( discro2@E_HSP, par( 'usr')[4], 'POP', col='orange', pos=1) # below
# abline( v=discro2@E_FSP, col='lightblue')
# text( discro2@E_FSP, par( 'usr')[4], 'FSP', col='lightblue', pos=1) # below
}
\keyword{misc}
