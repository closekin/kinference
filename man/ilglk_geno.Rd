\name{ilglk_geno}
\alias{ilglk_geno}
\title{Check individual genotypes for aggregate typicality}
\description{
Compute lglk of entire 4way genotype of each individual, ie sum log Pr[ g(i,l)]; and compare the distro of lglks across individuals with its predicted shape given ALFs. Significant mismatch is bad. Can also detect outliers. Up to you what criterion to use for that. NB running \code{locator(2)} lets you click the histogram to figure out where you'd like to cut.

Currently, the SPA calcs are a wee bit slow because of heavy use of \code{vecless} which in version 1.0 is sluggish. The lglks themselves are computed in C and are blisteringly fast.
}
\usage{ilglk_geno(snpg, indiv_lglk_hist_pars = list())
}
\arguments{
\item{ snpg}{a \code{snpgeno} (6way genotype)}
\item{ indiv_lglk_hist_pars}{list like in \code{dump_badhetz_fish}, for controlling histogram}
}
\value{Vector of lglk for each individual. I haven't added any formal uh-oh criteria yet; that could be done via the SPA, as in \code{dump_badhetz_fish}. But, reading off from the graph is probably fine...
}
\keyword{misc}
