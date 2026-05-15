\name{autopick_threshold}
\alias{autopick_threshold}
\title{PLOD threshold for HSPs}
\description{
This function proposes a PLOD threshold for excluding almost all 3KPs and 4KPs from true HSPs, and computes the associated false-negative probability (i.e., that a true HSP will have a PLOD below that threshold).

The modus operandi is to fit a mixture distribution to observed PLODs within some limited range that is expected to contain only 2KPs, 3KPs, and \emph{perhaps} a few 4KPs (and specifically excluding 1KPs). It is implicitly assumed that there \emph{is} a clear 2KP bump, and that it peaks close to the theoretical mean; if not, there is no point proceeding until you find more 2KPs (and perhaps tackle any QC issues). The threshold is then "chosen" (or proposed; the actual choice is really up to you) so that the expected number of false-positives from 3KPs plus 4KPs (i.e., with PLODs above the threshold) matches whatever you decide. A histogram with expected values is plotted (unless you tell it not to).

The mixture distribution treats the expected PLOD values for true 2KPs, 3KPs, and (if used) 4KPs as known (since those expectations can be calculated just from allele frequencies). However, because of linkage, the corresponding variances cannot be predicted a priori. Instead, the variance of the 2KP bump is estimated empirically from those PLODs which are above the 2KP mean, on the assumption that there will be few if any 3KPs with PLODs that large (hence the need for an upper limit on the included PLODs). Then, \code{\link{var_PLOD_kin}} (qv) is used to infer the extent of linkage and to predict what the 3KP and 4KP variances might be, under two extreme assumptions about the detailed nature of linkage. The two assumptions lead to upper and lower bounds on the variances of 3KPs and 4KPs. Next, a Normal mixture distribution is fitted by ML, to estimate the height and variance of the 3KP bump, subject to those bounds; it can optionally allow for the presence of 4KPs as well. Once those parameters have been estimated, \code{autopick_threshold} calculates the threshold PLOD above which \code{FPtol_pairs} of 3KPs (and perhaps 4KPs) would be expected. Finally, it calculates the proportion of 2KP PLODs that would be expected to fall below that threshold.

Since it is better to make the threshold too high rather than too low, you can use the \code{selecto} option to choose how the 3KP variance is estimated, subject to its bounds. The options are (i) strict ML estimation, to give the best fit to the overall PLODs, or (ii) choosing whichever 3KP variance leads to the \emph{highest} threshold (but still fitting the 3KP bump height by ML conditional on that variance).

The calculations in \code{\link{var_PLOD_kin}} are in fact geared specifically to HSPs, HTPs, and HC1Ps, as opposed to the other subtypes of 2/3/4KPs. Accordingly, it is a good idea to restrict the pairs to those which, if the are indeed 2KPs, are a priori likely to be HSP rather than GGP or HTP. (If they are not 2KPs, then it probably matters less exactly what subtype they are.) Having said that, the sky will \emph{probably} not fall in if some GGPs are also included; in some circumstances it might even be worth doing so deliberately, in order to increase the 2KP sample size. But the general principle of kin-finding applies: it is best to concentrate on the comparisons that yield the outcomes you're interested in, rather than including lots of extraneous comparisons and then trying to filter the results.

Despite the name \code{auto...}, \emph{you} still have to supply sensible values for a couple of parameters, based on looking at your data and understanding what you are trying to do. So it's not \emph{completely} automated process, and it never will be! Choosing a threshold is \emph{not} an "optimization process" with explicit bias/variance tradeoffs; rather, it's about ensuring that you have adequate "engineering tolerance" in the next stage of CKMR. The False-Negative Probability will, to a great extent, compensate for the choice of threshold (i.e. removing any bias in the fitted CKMR model) \emph{unless} you set the threshold too low, and thus end up with some 3KPs in your set of "definite 2KPs".
\subsection{Fitrange and use4th}{Assuming that you are using a 2KP-oriented PLOD, then the range of PLODs you fit to should extend from somewhere above 0 (which is very close to the expected PLOD for 3KPs), up to the RHS of the 2KP bump, but not so large as to include any FSPs and POPs. You can use the \code{fitrange_PLOD} parameter to control this, instead of subsetting the \code{kin} argument manually. If you set the lower limit high enough, then you don't need to worry about 4KPs intruding (their contribution would be negligible), so you can get away with fitting a 2-component mixture (simpler, less to go wrong...) by setting \code{use4th=FALSE}. But if you push the lower range closer to 0 (which does give you a larger sample size for fitting), then you might need to set \code{use4th=TRUE}. The (substantial) downside of doing that, is that there are often more PLODs close to 0 than near the 2KP mean, so the mixture-fit (which has to make more assumptions when also using 4KPs, and those assumptions may not be perfect) will concentrate its efforts on getting a good fit near 0, rather than near the 2KP mean which is what really matters. It is worth experimenting.
}
}
\usage{
autopick_threshold(
  x,
  kin,
  fitrange_PLOD,
  FPtol_pairs,
  use4th,
  selecto = c("ML", "paranoid"),
  NVAR = 10,
  plot_bins = NULL,
  shading_density = 10,
  want_all_results = FALSE,
  ...
)
}
\arguments{
\item{  x}{a \code{snpgeno} or its \code{locinfo} attribute. Must already have been prepared by running \code{kin_power}.}
\item{  kin}{a dataframe of "close-ish" kin-pairs and their PLODs, presumably from running \code{find_HSPs}; must have a column "PLOD".}
\item{  fitrange_PLOD}{two numbers, specifying the range of PLODs from \code{kin} to use in \emph{fitting} (though all are \emph{plotted}, by default)}
\item{  FPtol_pairs}{how many expected False-Positive 3KPs (and 4KPs, if \code{use4th=TRUE}) should the threshold exclude?}
\item{  use4th}{whether to allow for 4KPs when fitting.}
\item{  selecto}{whether to choose the threshold based on the best mixture fit ("ML"), or the most conservative ("paranoid").}
\item{  NVAR}{how many variances to try, between the limits set by \code{var_PLOD_kin}}
\item{  plot_bins}{bin-width for histogram plotting. Default NULL means no plot.}
\item{  shading_density}{By default, all PLODs in \code{kin} will be included in the histogram, even though only a subset are used in fitting. The histogram bars \emph{not} used in fitting (i.e., below \code{fitrange_PLOD[1]}) will be lightened in colour, according to this parameter. Results are graphics-device-dependent, so you may need to experiment away from the default; larger numbers usually mean lighter shading. Setting \code{shading_density=NA} should result in a light transparent rectangle covering the entire LHS of the graph, which you might prefer. You can also set \code{xlim} as usual, to remove those left-hand bars altogether.}
\item{  want_all_results}{if TRUE, return dataframe(s) containing results for each variance explored. This lets you examine "sensitivity".}
\item{  ...}{other parameters passed to \code{hist}, eg \code{xlim}, \code{ylim}, \code{col}. Many others will be ignored, and some will cause problems.}
}
\value{The proposed threshold, with lots of attributes. You can use those to calculate false-neg probabilities for \emph{other} possible thresholds, as per \bold{Examples}; you \emph{don't} have to accept the one that is proposed here! Threshold choice is \emph{up to you} (and badly-chosen thresholds are not the fault of \code{\link{kinference}})!
}
\seealso{\code{\link{kin_power}}, \code{\link{var_PLOD_kin}}
}
\examples{
library( atease) # for convenient 'thresh@info' below
dropbears1 <- kin_power( dropbears, k = 0.5)
hsps <- find_HSPs( dropbears1, keep_thresh = -10)
histoPLOD( hsps, log=TRUE) ## observe HSP - POP gap centred ~ PLOD = 120.
## Use that air-gap to set fitrange_PLOD.
thresh <- autopick_threshold( x= dropbears1, kin= hsps,
  fitrange_PLOD= c(0, 120), FPtol_pairs= 1, use4th= TRUE, plot_bins= 5)
thresh ## the 2nd order / 3rd order threshold value
thresh@info["Pr_FNeg"] ## ...
## ... the estimated 2KP false-neg rate, given that threshold
}
\keyword{misc}
