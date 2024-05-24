\name{histoPLOD}
\alias{histoPLOD}
\title{Histogram PLODs}
\description{This is a wrapper for the two old (and strangely-named, and too-HSP-specific) functions \code{\link{HSP_histo}} (qv) and \code{\link{PLOD_loghisto}} (qv), which are likely to become hidden inside the namespace so that \code{histoPLOD} will be the only way to call them. For now, see their help- while you still can. Eventually their help will move into this docu instead.
}
\usage{
histoPLOD(PLODs, log = c(FALSE, TRUE)[0],
    mean_show = c( 'POP', 'UP', if( !isFALSE( PLODs@trustplot_FSP_HSP_means)) c( 'HSP', 'FSP')),
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
\item{  mean_show}{for which kinships should expected PLODs values be shown? (Default won't show FSP or HSP "means" for \code{find_POPs_lglk}, since they are calculated wrongly for now.)}
\item{  main}{graph title, defaulting to the name of \code{PLODs} argument}
\item{  ...}{passed to the plotting routines}
\item{  UP_distro_show}{Only for \code{log=TRUE}, this controls which approximations to show for the theoretical PLODs distribution of UPs. In practice, they look damn similar! The remaining args only apply when \code{log=FALSE}:}
\item{  HSP_distro_show}{?show theoretical PLODs distro for HSPs, based on \emph{empirical}variance?}
\item{  lb, ub}{\code{lb} is \emph{mandatory} cutoff for which PLODs to include}
\item{  fullsib_cut}{only if \code{HSP_distro_show=TRUE}, then use this to determine which PLODs to include when calculating empirical variance of "HSP PLODs".}
\item{  bin}{bin width for histogram}
}
\seealso{\code{\link{PLOD_loghisto}}, \code{\link{HSP_histo}}
}
\keyword{misc}
