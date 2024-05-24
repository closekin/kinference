\name{FSP_POP_histo}
\alias{FSP_POP_histo}
\title{Plot for splitting FSPs from POPs}
\description{Plots an absolute-frequency histogram for the output of \code{find_FSPs_from_POPs_v2()}.
}
\usage{FSP_POP_histo(fsps2, bin = 100, FSPmean = TRUE, POPmean = TRUE, main = "", ...)
}
\arguments{
\item{  fsps2}{the output of a call to \code{find_FSPs_from_POPs_v2()}}
\item{  bin}{hist bin width. Used to define \code{breaks} (along with \code{xlim}, if given), so you can't manually pass in \code{breaks}.}
\item{  FSPmean}{plot the expected mean FP stat for FSPs? Default TRUE}
\item{  POPmean}{plot the expected mean FP stat for POPs? Default TRUE}
\item{  main}{ plot title. Passed to \code{hist.default}.}
\item{  ...}{additional pars, passed to \code{hist} (?should we force naming them in a list?)}
}
\seealso{PLOD_loghisto
}
\keyword{misc}
