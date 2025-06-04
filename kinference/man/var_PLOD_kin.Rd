\name{var_PLOD_kin}
\alias{var_PLOD_kin}
\title{Predict variance of PLOD for HTPs and HCPs}
\description{
\code{var_PLOD_kin} infers the extent of genetic linkage based on the empirical (i.e. observed) variance of PLODs among (presumed) true HSPs, by comparing it to the (smaller) variance that would theoretically apply if all loci were inherited independently. It then predicts the empirical variance for HTPs and/or HC1Ps, based on one of two extreme assumptions about the detailed nature of recombination. The truth will likely be between those extremes. \code{var_PLOD_kin} is called internally by \code{\link{autopick_threshold}} (qv), and few users will want to call it directly.

Although the biggest linkage effect comes from the main "order" of kinship (2KP, 3KP, etc, as controlled by the \code{n_meio} argument), this is one situation that is also affected somewhat by the "subtype" of kinship (HSP vs GGP vs FTP, etc). For CKMR purposes, the commonest type of kin of given order are typically those born closest in time, so for a given order of kinship, the \code{var_PLOD_kin} always assumes that the relevant kin-pairs will consist of the subtype with \emph{single} shared ancestor and minimax number of generations since shared ancestor. This means HSPs for \code{n_meio=2} (FTPs have 2 shared ancestors; GGPs entail 2 generations of gap, whereas HSPs have only 1), HTPs for \code{n_meio=3}, HC1Ps for \code{n_meio=4}, etc.

\code{var_PLOD_kin} can also be used for "simulation", e.g. to explore the effects of different types of linkage, or of different numbers of markers; see the \code{linfo} and \code{C_equiv} arguments, plus the code itself, for hints.
\subsection{Nature of linkage}{
The two extreme alternatives considered by \code{var_PLOD_kin} each use one parameter to capture "linkage", but the parameter has quite different meanings. The two versions are:

\itemize{
\item separate chromosomes, each containing the same number of loci, with \emph{no} recombination. The linkage parameter is the number of such chromosomes.
\item one single chromosome containing all the loci, equally spaced, with random recombination at multiple places. The linkaage parameter is the probability of recombination between adjacent loci.
}
}
}
\usage{
var_PLOD_kin(
  linfo,
  emp_V_HSP = V_noX(C_equiv, 2),
  n_meio,
  debug = FALSE,
  C_equiv = NULL
)
}
\arguments{
\item{  linfo}{either a \code{snpgeno} object, or its \code{$locinfo} attribute (or a fake one). Respectively, either \code{linfo$locinfo} or \code{linfo} itself should be a dataframe containing columns \code{e0}, \code{e1}, \code{v0}, \code{v1}, \code{count}. Each row is one "type" of locus, i.e., with roughly the same values of e/v 0/1, and \code{count} says how many such loci there are (for simulation purposes; with real data, \code{count==1}). \code{ e/v 0/1} are means and variances of the per-locus LOD (note no P) when the locus is or isn't co-inherited. See \emph{Details}}
\item{  emp_V_HSP}{empirical variance of PLOD among definite HSPs, normally from \code{autopick_threshold} (qv). You're supposed to be running this on real data, so that \code{emp_V_HSP} is an actual number; however, for testing purposes, you can set up an artificial version via \code{C_equiv} below.}
\item{  n_meio}{Target number of meioses:2 for HSPs, 3 for HTPs, 4 for HCPs, etc. This is by far the main driver of variance, but technically the only one; see \emph{Subtypes of kin}.}
\item{  debug}{Logical flag, not for regular folks to use! Defaults to FALSE. The \pkg{debug} package needs to be loaded beforehand.}
\item{  C_equiv}{for artificial test, with \code{emp_V_HSP} set to the no-crossover variance from \code{C_equiv} chromos (need not be integer). Ignored if \code{emp_V_HSP} is set.}
}
\section{Gruesome details}{The "per-locus LOD" (whose properties are stored in the columns \code{e0}, \code{e1}, \code{v0}, \code{v1} in \code{linfo}) is created by calling \code{\link{kin_power}} (qv). The normal use-case would be from \code{kin_power(...,k=0.5)}, so that the (P)LOD pertains to 2KP/UP comparisons. However, if you called \code{kin_power(...,k=0.25)} then the (P)LOD would be designed for 3KP/UP comparisons, and so on. In fact, you could even hand-tweak the calculations to contain LODs for 3KP/2KP comparisons, which \emph{might} in principle improve the resolution (but you'd have to fiddle manually; you could actually do it based on two calls to \code{\link{kin_power}}, one with \code{k=0.5} and one with \code{k=0.25}, and manipulating the results).
}
\value{Matrix with two rows \code{V0} and \code{Vx}, and one column for each element of \code{n_meio} (which is always augmented to include 2), named "M2" etc. The two rows pertain respectively to the no-crossover multiple-chromosome scenario, and the single-chromosome multiple-crossover scenario. The matrix also has an attribute \code{info}, which is a numeric vector of elements named \code{V_UP}, \code{V_HSP}, \code{C_hat}, and \code{rho_hat} (note that \code{V_HSP} should duplicate the first column of the matrix) \code{C_hat} is estimated equivalent number of chromosomes for the no-crossover scenario, and \code{rho_hat} is per-locus crossover rate for the all-crossover scenario.
}
\examples{
# COMPLETELY MADE-UP e/v values! Nothing to do with genetics :)
var_PLOD_kin( data.frame( count=45, ev01= I( cbind( e0=-1, e1=2, v0=0.03, v1=0.02))), C_equiv=22, n_meio=3:4)
#        M2    M3    M4
#  V0 208.2 156.6  91.9
#  Vx 208.2 186.2 101.8
#  attr(,"info")
#      V_UP    V_HSP    C_hat  rho_hat
#    1.3500 208.2273  22.0000   0.2616
}
\keyword{misc}
