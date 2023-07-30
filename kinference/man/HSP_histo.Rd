\name{HSP_histo}
\alias{HSP_histo}
\title{PLOD histogram}
\description{Plots an absolute-frequency histogram for the output of \code{find_kinpairs}, with the lower bound set by the user. Lower bounds should be set to exclude (as much as possible) the UP bump, as this will otherwise swamp the signal from the HSP bump. Users must manually set a lower bound for full-sibling PLODs (\code{fullsib_cut}) on order to exclude full-siblings from the variance estimate for HSP PLODs.
}
\usage{
HSP_histo(
  kinpairs,
  lb,
  ub = max(kinpairs$PLOD) + bin,
  fullsib_cut,
  bin = 5,
  HSPmean = TRUE,
  HSPdist = TRUE,
  POPmean = TRUE,
  FSPmean = TRUE,
  main = "",
  ...
)
}
\arguments{
\item{  kinpairs}{the output of a call to \code{find_kinpairs}}
\item{  lb}{PLOD lower bound for plot extent. Should exclude the UP bump}
\item{  ub}{PLOD upper bound for plot extent. Defaults to maximum PLOD score plus a little padding}
\item{  fullsib_cut}{PLOD score above which there are only full-sibs}
\item{  bin}{hist bin width. Default 5. lb, ub, and bin together define \code{breaks}, so you can't pass \code{breaks} via \code{...}}
\item{  HSPmean}{plot the mean PLOD for kinpairs? Default TRUE}
\item{  HSPdist}{plot the distribution of PLOD for kinpairs? Default TRUE}
\item{  POPmean}{plot the mean PLOD for POPs? Default TRUE}
\item{  FSPmean}{plot the mean PLOD for FSPs? Default TRUE}
\item{  main}{graph title, passed straight to \code{hist()}}
\item{  ...}{additional pars, passed to \code{hist()}}
}
\seealso{PLOD_loghisto
}
\keyword{misc}
