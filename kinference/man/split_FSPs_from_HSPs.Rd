\name{split_FSPs_from_HSPs}
\alias{split_FSPs_from_HSPs}
\alias{split_FSPs_from_POPs}
\alias{split_HSPs_from_HTPs}
\title{Discriminate between kinships of known close-kin}
\description{For pairs already picked as likely close-kin via one of the \code{find_XXX} functions, but whose exact kinship is uncertain (e.g. they might clearly be either POPs or FSPs, but it's not obvious which), these \code{split_XXX_from_YYY} functions apply a more powerful likelihood-based pairwise test statistic to each pair, to help decide what it is. All these functions use 4-way genotypes, i.e. not relying on 6-way genotyping. They possibly should be adapted to cope with 3-way (ie not trusting double-nulls) but currently they aren't (so they do trust double-nulls).
}
\usage{
split_FSPs_from_HSPs(snpg, candipairs)
split_FSPs_from_POPs( snpg, candiPOPs, gerr, use_obsolete_version=FALSE)
split_HSPs_from_HTPs(snpg, candipairs)
}
\arguments{
\item{  snpg}{a \code{snpgeno} object}
\item{  candipairs}{normally, a dataframe with rows being pairs and columns \emph{i} and \emph{j} (and possibly others) e.g. from \code{find_POPs} or \code{find_HSPs}. Can also be a 2-column matrix (each row again one pair).}
\item{  gerr}{genotyping error rate, where 0.01 would mean \code{1\%}. It's there to make the POP case robust; you have to choose it, but the precise value should not matter. See \bold{Details}.}
\item{  use_obsolete_version}{the original code for POPs-vs-FSPs was based on a different non-likelihood-based statistic; see \emph{Obsolete note}. It turned out to have low statistical power, but did not require specifying a \code{gerr}. For replicability purposes, you can still run it by setting this parameter to \code{TRUE}.}
}
\details{
The idea of \code{split_FSPs_from_POPs}- though this is not the only possible workflow- is that pairs which are \emph{either} POPs \emph{or} FSPs should stand out very clearly from everything else, via \code{\link{find_POPs}}. Then the job is to pick between those possibilities. The workflow is supposed to be:

 nail POPs/FSPs first with \code{\link{find_POPs}}`

\describe{
\item{pick between them with \code{split_FSPs_from_POPs} (update}{this doesn't work very well yet... use age info if you can)}
}

 look for HSPs (and potentially some HTPs) and filter out already-known POPs and FSPs

 filter out HTPs from the remaining set of HSPs with \code{split_HSPs_from_HTPs} and/or \code{\link{autopick_threshold}}

However, an equally reasonable workflow might be:

 nail HSPs and everything stronger (and potentially some HTPs) with \code{\link{find_HSPs}}

 split HSPs/HTPs from POPs/FSPs with \code{split_FSPs_from_HSPs}

 filter out HTPs from the remaining set of HSPs with \code{split_HSPs_from_HTPs} and/or \code{\link{autopick_threshold}}

 use \code{split_FSPs_from_POPs} to split FSPs from POPs

All \code{split_} functions return expected values under different possible kin-types (not variances, since these cannot be predicted for all kin-types).

The \code{gerr} parameter in \code{split_FSPs_from_POPs} is there to alleviate the problem that a single locus displaying apparent Mendelian exclusion is in theory reason enough to prove that a pair is \emph{not} a POP (if a likelihood-based criterion is used). But, of course, we can have genotyping errors (and, in a reasonably big dataset, mutations). Allowing for a small amount of error gives the method much mor flexibility, without paying a high price in statistical efficiency (provided \code{gerr} is small). The technical interpretation is that, if a genotyping error occurs at a locus, then the true genotype is replaced by a randomly-drawn genotype from the marginal distro of genotypes at that locus. Real genotyping errors don't work like that, but it is mathematically convenient and achieves the desired effect of robustifying the FSP-vs-POP statistic. The value to use is up to you; you can experiment; if you really want to estimate it, then look at replicate genotypes.
\subsection{Obsolete note}{
The statistic for \code{split_FSPs_from_POPs} - which isn't as powerful as I'd hoped; I'm going to redo it - is based on the weighted sum of the number of exactly-matching 4-way genotypes, with weights chosen to have high power for this particular discrimination. Weighting is optimized for the unlikely scenario that POPs and FSPs are equally likely a priori, but in practice the weights are not sensitive to this. The test is deliberately crude and robust- e.g. it avoids exclusion-based checks- on the assumption that you have enough loci to pick HSPs, so the more-related kin-types should be slam-dunks. \emph{But} it doesn't seem powerful enough. More worked needed...

\code{split_FSPs_from_HSPs} and \code{split_HSPs_from_HTPs} use 4-way genotypes only (to avoid having to worry about errors) but in a properly optimal PLOD designed for FSP/HSP or HSP/HTP discrimination- its expectation is positive for the higher-order kinship and negative for the lower-order kinship. Theoretical means for the two kinships named in the function's name are returned as attributes (variances cannot be predicted). Haven't added means for POPs or UPs since you're not "supposed" to have those in the mix by the time you run \code{split_*_from_*} functions, but maybe I should fix that at some point.
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
