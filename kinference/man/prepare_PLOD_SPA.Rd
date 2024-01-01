\name{prepare_PLOD_SPA}
\alias{prepare_PLOD_SPA}
\title{Prepare for kin-finding}
\description{\code{prepare_PLOD_SPA} is something you used to have to run before using some kin-finding/QC tools, to set up your \code{snpgeno} object for fancy maths woooo (saddlepoint approximations). There are no meaningful options, you just have to run this. It can be \emph{slightly} slow which is why it was a separate step. However, nowadays I don't think you need to run it at all, because it's built into \code{\link{kin_power}} (qv).
}
\usage{prepare_PLOD_SPA(geno6, n_pts_SPA_renorm= 201, sd_half_range= 10)
}
\arguments{
\item{  geno6}{a \code{snpgeno} object that has been thru \code{kin_power}}
\item{  n_pts_SPA_renorm}{how accurate to make the approximation. Default should be fine.}
\item{  sd_half_range}{How many SD's into the tails to push the approximation. The default of 10 is massively far. Normally this is fine, but if you get an error (probably from an NA cropping up during extreme calculations), then trying making it smaller.}
}
\value{Another \code{sngeno} object with an environment \code{Kenv}, which contains functions (with their own preloaded data) allowing null distributions (eg PLODs for true UPs) to be calculated. Various sanity checks are incorporated to try to stop you from stuffing up with out-of-synch loci etc later.
}
\keyword{misc}
