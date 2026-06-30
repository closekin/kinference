\name{kin_power}
\alias{kin_power}
\title{Locus selection for kin-finding}
\description{
\code{kin_power} is an essential preliminary step for various QC and PLOD calculations on real datasets; it calculates summary quantities for each locus based on allele frequency estimates. For real data, the results per se aren't of direct interest to most users. However, as its name suggests, \code{kin_power} can also be used to predict how well a proposed set of loci would work for finding 2KPs. (It can also do the same for 3KPs or other specified weaker kin, although the answer will always be "nope" for those cases anyway.)

\code{kin_power} returns its input \code{snpgeno} object after adding extra columns to the \code{locinfo} attribute, related to the per-locus mean and variance of LOD (presumably an HSP/UP PLOD, though not inevitably) for different true kinships. It respects the per-locus decision about how precisely to genotype (\code{useN=6/4/3}).
}
\usage{
kin_power( lociar, want_LOD_table=TRUE, k,
  hack_LOD= NULL, sd_half_range= 10)
}
\arguments{
\item{  lociar}{\code{snpgeno} objects with the necessary ingredients}
\item{  want_LOD_table}{can't think why you'd set this to FALSE}
\item{  k}{target average kinship for LOD; 0.5 for HSPs, 0.25 for HTPs, etc.}
\item{  hack_LOD}{Don't mess around with this; it's for internal black magic}
\item{  sd_half_range}{Normally leave this alone, but in case the final \code{prepare_PLOD_SPA} (qv) step gives an error (rare but possible), try reducing it below the default.}
}
\value{
A \code{snpgeno} object with augmented columns in \code{$locinfo}. Those are:
\item{  E_UP, V_UP}{mean & variance for UPs}
\item{  E_HSP, E_POP, E_FSP}{as you would expect}
\item{  Ediff}{\code{E_HSP-E_POP} ie the "absolute" power of that locus}
\item{  sdiff}{\code{(E_HSP-E_POP)/sqrt(V_UP)} which is arguably better than \code{Ediff} for ranking loci It also attaches \code{LOD}, \code{PUP}, and \code{ev01} elements (each a matrix) to \code{$locinfo}. They have been made dull (see \code{mvbutils::make_dull}) to improve your viewing experience, but they work fine for all normal purposes (and you can always \code{unclass} them to remove the S3 class \code{dull}).}
}
\examples{
## Next will fail cozza missing pre-calculated objects necessary for kin-finding:
try( hsps <- find_HSPs( dropbears, keep_thresh = 0))
dropbears_1 <- kin_power( dropbears, k = 0.5)
## works now
hsps <- find_HSPs( dropbears_1, keep_thresh = 0)
}
\keyword{misc}
