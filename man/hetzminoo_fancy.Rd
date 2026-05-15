\name{hetzminoo_fancy}
\alias{hetzminoo_fancy}
\title{QC checks on sample heterozygosity}
\description{
This test looks for samples with anomalous numbers of heterozygotes and/or double-nulls, which can result from (i) degraded DNA or (ii) sample contamination. Useful both for finding outlier samples, and for checking whether the loci are collectively working as they should (as is assumed by all the calculations in \code{\link{kinference}}). The histogram should coincide nicely with its predicted line.

The basis for the "hetzminoo" statistic is explained in the accompanying MS. Note that the \code{target} argument provides two variants, "rich" and "poor", designed to test for contamination and degradation respectively. With "rich", you should look for outlying samples on the \emph{RHS} of the histogram (too many heterozygotes): with "poor", on the \emph{LHS} (too few). The choice of target affects the weightings across loci, as explained in the MS. In practice, there is often little visual difference, and a bad sample looks bad bad in both. Nevertheless, you \emph{should} run both variants. See \code{kinference-vignette} (qv) for examples.
}
\usage{
hetzminoo_fancy(
  snpg,
  target = c("rich", "poor"),
  hist_pars = list(),
  showPlot = TRUE
)
}
\arguments{
\item{  snpg}{a \code{snpgeno} object, with allele frequencies already estimated and \code{kin_power} (qv) already run.}
\item{  target}{which potential problem to focus on.}
\item{  hist_pars}{list of parameters to pass to \code{hist}. If you are very sneaky, you can pass in an \code{expression} to be evaluated inline instead (ie fly-hacking). No, there's no example showing that!}
\item{  showPlot}{show the plot? Default TRUE}
}
\value{A vector of "hetzminoo" scores. You can then use it to subset your data by removing samples with unpleasant values.
}
\keyword{misc}
