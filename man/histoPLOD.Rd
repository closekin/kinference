\name{histoPLOD}
\alias{histoPLOD}
\title{Histogram PLODs}
\description{
\code{histoPLOD} shows a histogram of PLODs across pairwise comparisons, where PLOD has been calculated by \code{\link{find_HSPs}} (qv), \code{\link{find_POPs_lglk}} (qv), or something similar. The "histogram" can either be on a log-scale or the more familiar unlogged scale; the problem with the latter is that there are often so many UPs (we have hundreds of millions for tuna) that the interesting kin-bumps become squashed to the point of invisibility. The plots should show various bumps for the different kinships, overlapping to some extent; the theoretical locations of the bump-centres are shown by vertical lines for the most common kinships.

The typical workflow would be to first call with \code{log=TRUE}, to check that the UP bump looks good and that CK bumps (if any) are centred in the right places. For the UP bump only, the width can also be predicted theoretically, and so the entire predicted curve is shown (and it had better be a good match to the empirical distribution, otherwise there is a QC problem). Then you can focus in on interesting bits and actual pairs, by setting \code{log=FALSE} and zooming with the \code{lb} and \code{ub} parameters.

For \code{log=FALSE}, the lower bound \code{lb} should be set to exclude almost all of the UP bump, which will otherwise swamp the signal from the close-kin pairs. Specifically for HSPs or other 2KPs, it is also possible to show the entire "expected" distribution by setting \code{HSP_distro_show=TRUE}. However, unlike the UP bump, its variance and vertical scaling have to be calculated empirically. \code{histoPLOD} does that from PLODs between that mean and some upper limit \code{fullsib_cut}, which must be chosen manually. Hopefully the FSP/HSP gap is clear, so it won't matter much what you choose within that.

The \code{kinference-vignette} has examples.
}
\usage{
histoPLOD( PLODs, log = c(FALSE, TRUE)[0],
    mean_show = c( 'POP', 'UP',
        if( !isFALSE( PLODs@trustplot_FSP_HSP_means)) c( 'HSP', 'FSP')),
    UP_distro_show = c(SPA = TRUE, Normal = FALSE),
    HSP_distro_show = !log && !is.null( PLODs$mean_HSP),
    lb,
    ub = NULL,
    fullsib_cut = NULL,
    bin = 5,
    main = deparse1(substitute(PLODs), width.cutoff = 50),
    ...)
}
\arguments{
\item{  PLODs}{dataframe from \code{find_HSPs} or conceivably a future \code{find_kin3} etc}
\item{  log}{TRUE or FALSE to call \code{PLOD_loghisto} or \code{HSP_histo} respectively}
\item{  mean_show}{for which kinships should expected PLODs values be shown? There's no harm in showing all of them, but some kinships won't make sense in particular applications, so you can turn them off with this argument. (Note that the default won't show FSP or HSP "means" for \code{find_POPs_lglk}, since they are calculated wrongly at present.)}
\item{  main}{graph title, defaulting to the name of \code{PLODs} argument}
\item{  ...}{passed to the plotting routine, which is \code{plot} for \code{log=TRUE} and \code{hist} for \code{log=FALSE}.}
\item{  UP_distro_show}{Only for \code{log=TRUE}, this controls which approximations to show for the theoretical PLODs distribution of UPs. In practice, they look damn similar! The remaining args only apply when \code{log=FALSE}:}
\item{  HSP_distro_show}{Only for \code{log=FALSE}:show the theoretical PLODs distro for HSPs, based on \emph{empirical} variance and number of "definite" 2KPs.}
\item{  lb, ub}{only for \code{log=FALSE}. \code{lb} is a \emph{mandatory} cutoff for which PLODs to include. \code{ub} is optional, but can be used to visually exclude very high PLODs:eg to exclude POPs when the real interest lies in HSPs.}
\item{  fullsib_cut}{only if \code{log=FALSE} and \code{HSP_distro_show=TRUE}, then use this to determine which PLODs to include when calculating empirical variance of "HSP PLODs".}
\item{  bin}{bin width for histogram. Only for \code{log=FALSE} since most PLODs (which will be for UPs) are already binned during \code{find_HSPs} (qv).}
}
\keyword{misc}
