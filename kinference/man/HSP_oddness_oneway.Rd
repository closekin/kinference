\name{HSP_oddness_oneway}
\alias{HSP_oddness_oneway}
\alias{HSP_oddness_twoway}
\title{Oddness metrics over an HSP histo}
\description{Plots an absolute-frequency histogram for the output of \code{\link{find_HSPs}} (qv), with additional coloured regions showing individuals with odd-looking CLOD scores (see \code{check_FPosity()} ), ilglk stat (see \code{\link{ilglk_geno}}), or hetz stat (see \code{\link{hetzminoo_fancy}}). You should choose a PLOD range for plotting that excludes the UP bump, since UPs can swamp the signal from HSPs and HTPs. \code{HSP_oddness_oneway} shows cases where at least one member of the pair has an unusual score in the oddness metric, whereas \code{HSP_oddness_twoway} shows cases where both members do.
}
\usage{
HSP_oddness_oneway(hsps, snpg, lb, ub = max(hsps$PLOD) + 10, bin = 5,
  CLOD_prop = 0.001, ilglk_prop = 0.001, hetz_prop = 0.001, ...)
HSP_oddness_twoway(hsps, snpg, lb, ub = max(hsps$PLOD) + 10, bin = 5,
  CLOD_prop = 0.001, ilglk_prop = 0.001, hetz_prop = 0.001, ...)
}
\arguments{
\item{  hsps}{the output of a call to \code{find_HSPs()}}
\item{  snpg}{the \code{snpgeno} or \code{SPAgeno} object from which \code{hsps} was built}
\item{  lb}{PLOD lower bound for plot extent. Should exclude the UP bump}
\item{  ub}{PLOD upper bound for plot extent. Defaults to the maximum PLOD score plus a little padding}
\item{  bin}{hist bin width. Default 5. lb, ub, and bin together define \code{breaks}, so you can't pass \code{breaks} via \code{...}}
\item{  CLOD_prop}{the quantile of CLOD below which animals are highlighted. Default 0.001}
\item{  ilglk_prop}{the quantile of ilglk stat below which animals are highlighted. Default 0.001}
\item{  hetz_prop}{the quantile of hetz stat below which animals are highlighted. Default 0.001}
\item{  ...}{additional pars, passed to \code{hist()}}
}
\seealso{PLOD_oddness_oneway
}
\keyword{misc}
