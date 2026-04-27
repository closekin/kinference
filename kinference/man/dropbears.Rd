\name{dropbears}
\alias{dropbears}
\docType{data}
\title{Red-rumped Dropbear data}
\description{
An anonymised \code{snpgeno} dataset of Red-rumped Dropbear (\emph{Thylarctos plummetus}, ssp. \emph{haemorhous}) genotypes for 480 individuals at 2000 loci. Owing difficulties retaining field staff for dropbear research, these data represent only a single sampling year.

Genotypes are "4-way", ie single-nulls and homozygotes are not distinguished, and double-nulls are potentially allowed for, although there aren't any in this dataset.

Allele frequencies were estimated using \code{\link{est_ALF_ABO_quick}} (qv) which \emph{does} estimate a null-allele frequency. It's arguable whether that was the best choice, given that nulls seems so rare; we might have used \code{\link{est_ALF_nonulls}} (qv) instead. See \code{kinference-vignette} for more.

\preformatted{
  diplos( dropbears) # 4-way encoding
  dropbears[ 1:5, 1:3] # subset of covariates and genotypes
  sum( dropbears=='OO') # no double-nulls...
  sum( dropbears=='AB') # ... but lots of heterozygotes!
}
}
\usage{dropbears
}
\format{An object of class \code{gbasics::snpgeno()}
}
\author{Shane M Baylis <email: shane.baylis@csiro.au>}
\keyword{data}
