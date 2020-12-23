\name{var_PLOD_kin}
\alias{var_PLOD_kin}
\title{Predict variance of PLOD for HCPs and HTPs}
\description{Aim is to work out how much your putative half-sibling pairs (HSPs) might be contaminated by half-cousin pairs (HCPs) or half-thiatic pairs (HTPs). HSP-selection is presumably based on the pairwise PLODs for HSP:UP, taking all pairs where that PLOD exceeds some threshold. Given the allele freqs, the mean PLOD is predictable when the truth is UP, HCP, HTP, or HSP. The variance is only predictable for UPs, because linkage makes loci non-independent for kin. However, an empirical variance can be estimated for HSPs based on the observed PLODs above some super-high threshold, e.g., the mean PLOD when truth=HSP (this preliminary variance step could be iterated using different super-high thresholds). Based on the empirical variance for HSPs and the analytical variance for UPs, we basically know how much linkage there might be, so we can predict the PLOD variances for the other kin-pair types. The wrinkle is that it also depends somewhat on the finer-scale organization of the genome, i.e. whether it's lots of chromosomes with no crossover, or a few chromosomes with lots of crossover. \code{var_PLOD_kin} therefore calculates two versions, one assuming the genome is entirely made up of equal-sized chromosome with zero crossover, and the other assuming the genome is a single chromosomes with crossover according to a random process.
}
\usage{
var_PLOD_kin(linfo, emp.V.HSP = V.noX(C.equiv, 2), kin.true = c("HCP",
  "HTP", "HSP"), debug = FALSE, C.equiv = NULL)
}
\arguments{
\item{  linfo}{locus info, a \code{data.frame} with columns \code{e0}, \code{e1}, \code{v0}, \code{v1}, \code{count}. Each row is one "type" of locus, i.e., with roughly the same values of e/v 0/1, and \code{count} says how many such loci there are. e/v 0/1 are means and variances of the per-locus LOD (note no P) when the locus is or isn't co-inherited.}
\item{  emp.V.HSP}{empirical variance of PLOD for deffo HSPs. You're supposed to be running this on real data, so that \code{emp.V.HSP} is an actual number; however, for testing purposes, you can set up an artificial version via \code{C.equiv} below.}
\item{  kin.true}{Var[PLOD|kin.true]. The default, \code{"HCP"}, corresponds to 4 meioses, \code{"HTP"} to 3, and \code{"HSP"} to 2; the latter is only for debugging, since it should reproduce the original empirical variance!}
\item{  debug}{Logical flag. Defaults to FALSE.}
\item{  C.equiv}{for artificial test, with \code{emp.V.HSP} set to the no-crossover variance from \code{C.equiv} chromos (need not be integer). Ignored if \code{emp.V.HSP} is set.}
}
\details{
It's assumed that lots of loci are being used, so that the mix of loci on each "chromo", or the splatter of loci along the single "megachromo", always matches the overall population, on law-of-large-numbers grounds.

Stuff like uncertainly in allele frequencies, and in the PLOD variance for HSPs, needs to be accounted for externally, by repeatedly drawing from the posteriors and re-calculating the PLODs and re-running this function.

If the variance estimates show really good separation between the kin-pair types, then one could refine the "preliminary variance" step by reducing the super-high threshold (and assuming a truncated-Normal distribution). This might be worthwhile if the preliminary variance otherwise has to be based on a very small number of no-brainer HSPs. The "logical conclusion" of That Kind Of Thing is some kind of MLE involving estimating the population of different types of kin, and we really don't want to go there for now (since that should include the population dynamics shebang). In other words, we'd end up linking the genetic kin-finding model to the population dynamics model, which makes life statistically harder. And god knows it's hard enough. Anyway, if we were taking that approach, it might well be better to avoid PLODs altogether and instead go for inferences about the actual ppn of co-inherited loci, from which estimates-of-co-inherited-variance and inferences about kin-ppns can be made.
}
\value{Vector with names \code{V0}, \code{Vx}, \code{C.hat}, \code{rho.hat}, \code{n.meio.<XXX>}. First two are variances under the no- and all-crossover scenarios; \code{C.hat} is estimated equivalent number of chromosomes, and \code{rho.hat} is per-locus crossover rate, under the same scenarios. \code{n.meio.<XXX>} (where \code{"XXX"} is set to \code{kin.true}) shows how many meioses are involved; yes, a perverse way to return that piece of info.
}
\examples{
\dontrun{
# H1CPs, for "simulation" equiv to 22 no-Xover chromos
var_PLOD_kin( data.frame( count=45, e0=-1, e1=2, v0=0.03, v1=0.02), C.equiv=22)
#        VUP       V.HSP          V0          Vx       C.hat     rho.hat n.meio.HCP
#     1.3500    208.2273     91.9010    101.8270     22.0000      0.2616      4.0000
# HTPs
var_PLOD_kin( data.frame( count=45, e0=-1, e1=2, v0=0.03, v1=0.02), C.equiv=22, kin='HTP')
#       VUP      V.HSP         V0         Vx      C.hat    rho.hat n.meio.HTP
#    1.3500   208.2273   156.5642   186.1779    22.0000     0.2616     3.0000
# Vars would go up with fewer equiv chromos
var_PLOD_kin( data.frame( count=45, e0=-1, e1=2, v0=0.03, v1=0.02), C.equiv=5, kin='HTP')
#       VUP      V.HSP         V0         Vx      C.hat    rho.hat n.meio.HTP
#   1.35000  912.37500  684.67500  786.47854    5.00000    0.04951    3.00000
}
}
\keyword{misc}
