\name{ilglk_geno}
\alias{ilglk_geno}
\title{Check individual multilocus genotypes for typicality}
\description{
\code{ilglk_geno} computes the per-sample log-likelihood across its entire genotype, i.e. sum log Pr g(i,l); and compares the distribution across individuals to the theoretical distribution given allele frequencies. Some mismatch is normal (and can arise just from noise in allele-frequency estimates), but substantial mismatch is bad. You get to define "substantial". \code{ilglk_geno} can also detect outlier individuals, usually with lglks that are much too low rather than too high; it's not obvious what could generate genomes that are "too typical" at the individual level.

You can use \code{locator(1)} to click the histogram to figure out where to adjust the \code{xlim/ylim} values to change the range of the data to inspect more closely. In other words, you can then re-run the function with its \code{...hist_par} argument set accordingly.
}
\usage{ilglk_geno(snpg, hist_pars = list(), showPlot = TRUE)
}
\arguments{
\item{  snpg}{a \code{snpgeno}. Genotype encoding (see \code{\link{diplos}}) must be one of 4-way, 3-way, or 6-way, as determined by \code{diplos} and the \code{$locinfo$useN} field.}
\item{  hist_pars}{\code{list()} passed to \code{hist} for controlling histogram, e.g. \code{hist_pars=list(xlim=c(-12000, -6000))}, or use \code{FALSE} to not plot.}
\item{  showPlot}{show the histogram? Defaults to TRUE, but overrideen by \code{hist_pars=FALSE}.}
}
\value{Vector of log-likelihood for each individual; also usually (but optionally), a histogram of log-likelihood values across individuals.
}
\examples{
## get rid of really bad loci
pvals <- check6and4( bluefin, thresh_pchisq_6and4 = c( 0.001, 0.0001))
bluefin_1 <- bluefin[ , pvals$pval4 > 0.01] # drastic QC!
## check for samples that are not like the others
ilglks <- ilglk_geno( bluefin_1)
## looks like anything with a lglk < -1030 is definitely abnormal
bluefin_2 <- bluefin_1[ ilglks > -1030,]
ilglks <- ilglk_geno( bluefin_2) ## much better, but not perfect -
## see other cleaning steps in the vignette
}
\keyword{misc}
