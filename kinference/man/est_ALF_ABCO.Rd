\name{est_ALF_ABCO}
\alias{est_ALF_ABCO}
\alias{re_est_ALF}
\title{Estimate allele frequencies olde-style, including nulls}
\description{
Most users should avoid these! \code{re_est_ALF} calculates maximum-likelhood estimates of minor, null, and 3rd-allele freqs from a \code{snpgeno} object. Its workhorse is \code{est_ALF_ABCO} which has a weirder syntax unless you are using an old CSIRO "pipeline", and is documented here for that reason only.

\code{re_est_ALF} will accept genotypes 4-way genotypes including "AAO" and "BBO", 6-way genotypes (but there are better options in that case; see below), or triallelic "ABCO" genotypes. The latter corresponds to \code{genotypes_ambig} as seen in the code of \code{define_genotypes}; they allow for an optional 3rd allele and nulls, but do not distinguish between single-nulls and homozygotes. Historically, for CSIRO users, "ABCO"-type genotypes are produced by \code{genocalldart::geno_deambig_ABC}.

Null-allele frequency has to be estimated from HWE deviations, so good estimates require a decent sample size.

The original use-case for these functions was datasets where bona fide 3rd alleles are common; even though they are not used in any \code{\link{kinference}} step (because they get recoded to nulls, i.e. "neither A nor B" en route), it's useful to have them around for ALF estimation. There's a bit of history behind this, which I shan't go into here.

If there are no 3rd alleles, then \code{re_est_ALF} is not the best choice. If you really do have 6-way genotypes (i.e. differentiating single-nulls from homozygotes, at least approximately) then you \emph{could} use \code{re_est_ALF} but the problems are:

\itemize{
\item \code{re_est_ALF} does not use the extra statistical information on single-null vs homozygote that is available with 6-way genotypes, whereas \code{\link{est_ALF_6way}} does;
\item and in the absence of that information, it's much slower than \code{\link{est_ALF_ABO_quick}} would be, so why not use that instead?!
}
}
\usage{
re_est_ALF( snpg)
est_ALF_ABCO( lociar, geno_amb = attr( lociar, 'geno_amb'))
}
\arguments{
\item{  snpg}{a \code{snpgeno} object.}
\item{  lociar}{a \code{loc.ar} or \code{snpgeno} object, normally with a \code{geno_amb} attribute from 4-way genotyping (see next).}
\item{  geno_amb}{a set of 4-way genotypes. In CSIRO's 6-way genotyping pipeline (in package \pkg{genocalldart}, not for general use), 4-way genotypes get called first and stored in a \code{geno_amb} attribute, before making initial allele frequency estimates with \code{est_ALF_ABCO}. those are then used as starting values in 6-way genotype-calling and allele-frequency estimation, and the \code{geno_amb} attribute gets discarded..}
}
\value{\code{re_est_ALF} returns the input, adding a 4-column matrix \code{pbonzer} to the \code{$locinfo} attribute, plus attributes \code{gobs} and \code{gpred} showing observed and expected counts of each genotype per locus. \code{pbonzer} is what you want for subsequent calculations. \code{est_ALF_ABCO} is similar but creates a new attribute \code{$locinfo$pambig}, instead of creating/modifying \code{pbonzer}.
}
\seealso{\code{\link{est_ALF_ABO_quick}} for normal users, \code{\link{est_ALF_6way}} for special 6-way people.
}
\examples{
head( bluefin$pbonzer) # C alleles sometimes present
bluefin$locinfo$pbonzer <- NULL ## remove pre-existing ALFs
bluefin$locinfo$snerr <- NULL ## remove pre-existing snerr
bluefin <- re_est_ALF( bluefin)
head( bluefin$locinfo$pbonzer) # C-alleles have gone to 0
}
\keyword{misc}
