\name{check6and4}
\alias{check6and4}
\title{Locus QC check}
\description{
Compares 6-way (if available) and 4-way genotype counts to HWE expectations, based on already-estimated allele frequencies. Plots a histogram of p-values across all loci, and (more importantly) plots observed / expected counts for each genotype by locus.

An overall p-value is calculated for each locus, based on chisq(1) of the G-statistic; it is not obvious how many DoF should be used. Anyway, the p-value itself is just indicative; few if any CKMR datasets actually "fit properly" because sample sizes are so large that even trivial misfits lead to significant departures from HWE null hypothesis, even though the final outcomes of kin-finding (ie PLODdograms) can look perfectly good. This is usually obviuos from the initial histogram of p-values, which in theory should be uniform for the "good" loci; in practice, it usually has heavy left-skew. Thus, the p-values are not to be taken literally, but just as a relative guide; the worst loci have the smallest p-values, and you can use the p-value as a criterion for dropping some loci. Since QC is iterative (samples, loci, samples, loci, ...) you can always revisit the decision later.

The p-value threshold (\code{thresh_pchisq_6and4}) simply determines the colour used to plot each locus. You can manually changing the threshold and re-run until you have a visually satisfactory pattern of orange vs green. Green means good (or good enough); orange means bad, ie below the lower threshold. (You are actually allowed to supply two numbers for the threshold, in which case loci with in-between values will display pink; that is probably a design flaw, because it's confusing to have more than one threshold). Once you have identified a threshold that \emph{looks good}, then you can use the \code{pval4} or \code{pval6} return-value to keep or discard loci, depending whether they are above or below that threshold.

The plots never look absolutely perfect, and there is no absolute criterion for "how bad is too bad". So, judgement and experience are required. But remember that keep-or-drop decisions aren't final; the whole QC process can be somewhat iterative, and the ultimate test is the PLODdogram(s) at the end. If it looks bad, ie bumps in the wrong places, then you may not have been restrictive enough (ie your threshold might be too high); whereas if the bumps are in the right places but are too wide for kin-finding, then your threshold might be too low.

With 6-way genotyping (eg SBTuna), calculations are done separately for the 6-way and 4-way versions. It is quite possible for a locus to look bad in the 6-way version, but good in the 4-way version; if so, don't just throw it out entirely, but try setting \code{useN=4} for that locus, or \code{useN=3} if null frequency is dangerously low.

For examples, see the \code{kinference-vignette}.
}
\usage{
check6and4(
  geno6,
  thresh_pchisq_6and4,
  return_what= c("just_pvals", "all"),
  extra_title= "",
  show6= FALSE
)
}
\arguments{
\item{  geno6}{a \code{snpgeno} object with 4-way (or optionally 6-way) genotypes}
\item{  thresh_pchisq_6and4}{a pair of thresholds for "bad" and "really bad" p-values. These determine the color in which each locus appears in all subplots.}
\item{  return_what}{one of \code{just_pvals} or \code{all}; see value}
\item{  extra_title}{a character string to be added to the bottom-right corner of all plots. Best if < 25 characters.}
\item{  show6}{show the plots for 6-way goodness-of-fit? Defaults to TRUE. If diplos is anything other than genotypes6, should be FALSE.}
}
\value{Creates per-locus vectors \code{pval6} and \code{pval4} for 6-way and 4-way genotypes respectively. If \code{return_what="just_pvals"}, these are returned in a list; if \code{return_what="all"}, they are added as columns to \code{geno6$locinfo}.
}
\keyword{misc}
