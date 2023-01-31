\name{ilglk_geno}
\alias{ilglk_geno}
\title{Check individual multilocus genotypes for typicality}
\description{\code{ilglk_geno} computes log-likelihood of entire 4-way (not 6-way) genotype of each individual, i.e., sum log Pr[ g(i,l)]; and compares the distribution across individuals to theoretical distro given allele frequencies. Significant mismatch is bad. Can also detect outlier individuals, usually with lglks that are much too low (ie rather than too high- I'm not sure what could generate "too typical a genome" at the individual level).
}
\usage{ilglk_geno(snpg, hist_pars = list(), showPlot = TRUE)
}
\arguments{
\item{  snpg}{a \code{snpgeno} (6-way genotype)}
\item{  hist_pars}{\code{list()} passed to \code{hist} for controlling histogram, eg \code{hist_pars=list(xlim=c(-12000, -6000))}, or use \code{FALSE} to not plot.}
\item{  showPlot}{show the histogram? Defaults to TRUE, but overrideen by \code{hist_pars=FALSE}.}
}
\details{
You can use \code{locator(1)} to click the histogram to figure out where to adjust the \code{xlim}/'ylim' values to change the range of the data to inspect more closely- ie you then re-run the function with its \code{...hist_par} argument set accordingly.

Currently, the SPA calcs are a wee bit slow because of heavy use of \code{vecless} which in version 1.0 is sluggish. The lglks themselves are computed in C and are blisteringly fast. If the SPA line (expected distro) doesn't appear, let us know- needs fixing! There might e.g. be too many loci, so that the calculation is falling over.

Haven't added any formal uh-oh criteria yet; that could be done via the SPA, as in \code{dump_badhetz_fish}. But, reading off from the graph is probably fine...
}
\value{Vector of log-likelihood for each individual; also usually (but optionally), a histogram of log-likelihood values across individuals.
}
\examples{
## Need an example!
}
\keyword{misc}
