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
data( bluefin)
## stripped-down data-cleaning for example - see the
## vignette for approach for real data!
pvals <- check6and4( bluefin, thresh_pchisq_6and4 = c( 0.001, 0.0001))
bluefin_1 <- bluefin[ , pvals$pval4 > 0.01]
ilglks <- ilglk_geno( bluefin_1)
bluefin_2 <- bluefin_1[ ilglks > -1050,]
bluefin_3 <- est_ALF_ABO_quick( bluefin_2)
bluefin_4 <- bluefin_3[ , bluefin_3@locinfo$pbonzer[,"B"] > 0.02]

bluefin_5 <- kin_power( bluefin_4, k = 0.5)
dups <- find_duplicates( bluefin_5, max_diff_loci = 20)
bluefin_6 <- bluefin_5[ -c(drop_dups_pairwise_equiv( dups)) ,]

## via find_HSPs for isolated split of 'FSPsOrPOPs'
hsps <- find_HSPs( bluefin_6, keep_thresh = 0)

histoPLOD( hsps, log = FALSE, lb = 0, fullsib_cut = 70, HSP_distro_show
= TRUE)
## note the gap between PLOD = 60 and PLOD = 75, separating likely HSPs
## from likely POPs and FSPs
FSPsOrPOPs <- hsps[ hsps$PLOD > 70,]
HSPsOrHTPs <- hsps[ hsps$PLOD < 70,] ## implicit >0 from keep_thresh in
## find_HSPs call

## via find_POPs_lglk for 'laser-focus on POPs'
maybePOPs <- find_POPs_lglk( bluefin_6, gerr = 0.01, keep_thresh = -30)
hist(maybePOPs$PLOD, breaks = 20)
abline( v = maybePOPs@mean_POP, col = kinPalette()["POP"], lwd = 2)
abline( v = maybePOPs@mean_UP, col = kinPalette()["UP"], lwd = 2)
## note the 6 pairs clustered around the expected mean PLOD for
## POPs (all > 60). They are probably POPs, but may include FSPs.
maybePOPs <- maybePOPs[ maybePOPs$PLOD > 60,]

## split_FSPs_from_POPs
### using FSPsOrPOPs
splitFSPsOrPOPs <- split_FSPs_from_POPs(bluefin_6, candiPOPs = FSPsOrPOPs, gerr = 0.01)
hist(splitFSPsOrPOPs$PLOD)
## add line of expected PLOD if FSP
abline( v = splitFSPsOrPOPs@E_FSP, lwd = 2, col = kinPalette()["FSP"])
## add line of expected PLOD if POP
abline( v = splitFSPsOrPOPs@E_POP, lwd = 2, col = kinPalette()["POP"])
## The expected PLOD for true POPs is negative, and for true FSPs is positive

### using maybePOPs
splitmaybePOPs <- split_FSPs_from_POPs(bluefin_6, candiPOPs = maybePOPs, gerr = 0.01)
hist(splitmaybePOPs$PLOD)
## add line of expected PLOD if FSP
abline( v = splitmaybePOPs@E_FSP, lwd = 2, col = kinPalette()["FSP"])
## add line of expected PLOD if POP
abline( v = splitmaybePOPs@E_POP, lwd = 2, col = kinPalette()["POP"])
## The expected PLOD for true POPs is negative, and for true FSPs is positive

## split_HSPs_from_HTPs
splitHSPsOrHTPs <- split_HSPs_from_HTPs( bluefin_6, HSPsOrHTPs)
hist(splitHSPsOrHTPs$PLOD, breaks = 15)
## add line of expected PLOD if HSP
abline( v = splitHSPsOrHTPs@E_HSP, lwd = 2, col = kinPalette()["HSP"])
## add line of expected PLOD if HTP
abline( v = splitHSPsOrHTPs@E_HTP, lwd = 2, col = kinPalette()["HTP"])
## the expected PLOD for true HSPs is positive, and for true HTPs is negative
}
\keyword{misc}
