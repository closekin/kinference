#' Predict variance of PLOD for HCPs and HTPs
#'
#' Aim is to work out how much your putative half-sibling pairs (HSPs) might be contaminated by half-cousin pairs (HCPs) or half-thiatic pairs (HTPs). HSP-selection is presumably based on the pairwise PLODs for HSP:UP, taking all pairs where that PLOD exceeds some threshold. Given the allele freqs, the mean PLOD is predictable when the truth is UP, HCP, HTP, or HSP. The variance is only predictable for UPs, because linkage makes loci non-independent for kin. However, an empirical variance can be estimated for HSPs based on the observed PLODs above some super-high threshold, e.g., the mean PLOD when truth=HSP (this preliminary variance step could be iterated using different super-high thresholds). Based on the empirical variance for HSPs and the analytical variance for UPs, we basically know how much linkage there might be, so we can predict the PLOD variances for the other kin-pair types. The wrinkle is that it also depends somewhat on the finer-scale organization of the genome, i.e. whether it's lots of chromosomes with no crossover, or a few chromosomes with lots of crossover. \code{var.PLOD.kin} therefore calculates two versions, one assuming the genome is entirely made up of equal-sized chromosome with zero crossover, and the other assuming the genome is a single chromosomes with crossover according to a random process.
#'
#' It's assumed that lots of loci are being used, so that the mix of loci on each "chromo", or the splatter of loci along the single "megachromo", always matches the overall population, on law-of-large-numbers grounds.
#'
#' Stuff like uncertainly in allele frequencies, and in the PLOD variance for HSPs, needs to be accounted for externally, by repeatedly drawing from the posteriors and re-calculating the PLODs and re-running this function.
#'
#' If the variance estimates show really good separation between the kin-pair types, then one could refine the "preliminary variance" step by reducing the super-high threshold (and assuming a truncated-Normal distribution). This might be worthwhile if the preliminary variance otherwise has to be based on a very small number of no-brainer HSPs. The "logical conclusion" of That Kind Of Thing is some kind of MLE involving estimating the population of different types of kin, and we really don't want to go there for now (since that should include the population dynamics shebang). In other words, we'd end up linking the genetic kin-finding model to the population dynamics model, which makes life statistically harder. And god knows it's hard enough. Anyway, if we were taking that approach, it might well be better to avoid PLODs altogether and instead go for inferences about the actual ppn of co-inherited loci, from which estimates-of-co-inherited-variance and inferences about kin-ppns can be made.
#'
#' @param linfo locus info, a \code{data.frame} with columns \code{e0}, \code{e1}, \code{v0}, \code{v1}, \code{count}. Each row is one "type" of locus, i.e., with roughly the same values of e/v 0/1, and \code{count} says how many such loci there are. e/v 0/1 are means and variances of the per-locus LOD (note no P) when the locus is or isn't co-inherited.
#' @param emp.V.HSP empirical variance of PLOD for deffo HSPs. You're supposed to be running this on real data, so that \code{emp.V.HSP} is an actual number; however, for testing purposes, you can set up an artificial version via \code{C.equiv} below.
#' @param kin.true Var[PLOD|kin.true]. The default, \code{"HCP"}, corresponds to 4 meioses, \code{"HTP"} to 3, and \code{"HSP"} to 2; the latter is only for debugging, since it should reproduce the original empirical variance!
#' @param C.equiv for artificial test, with \code{emp.V.HSP} set to the no-crossover variance from \code{C.equiv} chromos (need not be integer). Ignored if \code{emp.V.HSP} is set.
#'
#' @return Vector with names \code{V0}, \code{Vx}, \code{C.hat}, \code{rho.hat}, \code{n.meio.<XXX>}. First two are variances under the no- and all-crossover scenarios; \code{C.hat} is estimated equivalent number of chromosomes, and \code{rho.hat} is per-locus crossover rate, under the same scenarios. \code{n.meio.<XXX>} (where \code{"XXX"} is set to \code{kin.true}) shows how many meioses are involved; yes, a perverse way to return that piece of info.
#'
#' @examples
#' \dontrun{
#' # H1CPs, for "simulation" equiv to 22 no-Xover chromos
#' var.PLOD.kin( data.frame( count=45, e0=-1, e1=2, v0=0.03, v1=0.02), C.equiv=22)
#' #        VUP       V.HSP          V0          Vx       C.hat     rho.hat n.meio.HCP
#' #     1.3500    208.2273     91.9010    101.8270     22.0000      0.2616      4.0000
#'
#' # HTPs
#' var.PLOD.kin( data.frame( count=45, e0=-1, e1=2, v0=0.03, v1=0.02), C.equiv=22, kin='HTP')
#' #       VUP      V.HSP         V0         Vx      C.hat    rho.hat n.meio.HTP
#' #    1.3500   208.2273   156.5642   186.1779    22.0000     0.2616     3.0000
#'
#' # Vars would go up with fewer equiv chromos
#' var.PLOD.kin( data.frame( count=45, e0=-1, e1=2, v0=0.03, v1=0.02), C.equiv=5, kin='HTP')
#' #       VUP      V.HSP         V0         Vx      C.hat    rho.hat n.meio.HTP
#' #   1.35000  912.37500  684.67500  786.47854    5.00000    0.04951    3.00000
#' }
#' @export
"var.PLOD.kin" <-
function( linfo, emp.V.HSP=V.noX( C.equiv, 2),
         kin.true=c( 'HCP', 'HTP', 'HSP'), debug=FALSE, C.equiv=NULL) {

  kin.true <- match.arg( kin.true)
  n.meio <- c( HCP=4, HTP=3, HSP=2)[ kin.true]

# linfo should be a DF with cols e0, e1, v0, v1, count
  names( linfo) <- names( linfo) %&% '.l'
  extract.named( linfo)
  L <- sum( count.l)
  pi <- count.l / L

  e0 <- pi %*% e0.l
  e1 <- pi %*% e1.l
  v0 <- pi %*% v0.l
  v1 <- pi %*% v1.l

  # noX = no crossover, i.e. all in separate equal chromos
  V.noX <- function( C, meioses) {
      p <- 2 ^ (1-meioses) # Marginal prob of coinheritance = 1/2 for HSPs, 1/8 for HCPs
      EofV <- L * v0 + L * p *(v1-v0)
      VofE <- sqr( L * (e1-e0)) * p * (1-p) / C
    return( EofV + VofE)
  }

  # MoM for HSPs:
  if( !is.null( C.equiv)) {
    C.hat <- C.equiv <- min( L, C.equiv)
  } else {
    C.hat <- if( V.noX( 1, 2) < emp.V.HSP) 1 else
        if( V.noX( L, 2) > emp.V.HSP) L else
        find.root( V.noX, target=emp.V.HSP, start=1, step=1, fdirection='decreasing',
            min.x=1, max.x=L, meioses=2)
  }
  V0 <- V.noX( C.hat, meioses=n.meio)

  # Entirely Xover on one single chromo
  # Should these be done separately by locus type, then averaged??
  ell <- 1 %upto% (L-1)
  e2.1 <- pi %*% (sqr( e1.l) + v1.l)
  e2.0 <- pi %*% (sqr( e0.l) + v0.l)
  e1.0 <- e0 # pi %*% e0.l
  e1.1 <- e1

  phi <- function( r, mul, s) 0.5 * (1 + s*exp( -mul*r*ell))
  PSstar <- P0star <- array( 0, c( L-1, 2, 2, 2))
  si2 <- slice.index( PSstar, 2)
  si3 <- slice.index( PSstar, 3)
  si4 <- slice.index( PSstar, 4)

  V.allX <- function( rho, meioses) {
    # PSstar[ L,i,j,k] is Pr[ composite state @ L = 1 | starting states @ 0 are i, j, k]
    # j & k series irrel for HSP
    # slice.index() == 1 or 2
    # so 3 - 2*s.i. == -1 or +1
    # and 2-s.i. == 1 or 0

    # This code could be much more efficient; pre-compute P0star, and much reduce the phi-calcs
    PSstar[] <<- phi( rho, 4, 3 - 2*si2) *
        (if( meioses>2) phi( rho, 2, 3 - 2*si3) else 1) *
        (if( meioses>3) phi( rho, 2, 3 - 2*si4) else 1)

    P0star[] <<-  (2-si2) *
        (if( meioses>2) 2-si3 else 1) *
        (if( meioses>3) 2-si4 else 1)

    p11.8 <- apply( PSstar * P0star, 1, sum)
    p01.8 <- apply( PSstar * (1-P0star), 1, sum)
    p00.8 <- apply( (1-PSstar) * (1-P0star), 1, sum)

    p11 <- p11.8 / 8
    p01.10 <- p01.8 / 4 # both ways
    p00 <- p00.8 / 8

    p <- 2 ^ (1-meioses) # Marginal prob of coinheritance = 1/2 for HSPs, 1/8 for HCPs

    EL <- L * ( e1.1*p + e1.0*(1-p))
    EL2 <- L * ( e2.1*p + e2.0*(1-p)) +
        2 * (L-ell) %*% ( sqr( e1.1)*p11 + e1.1*e1.0*p01.10 +  sqr( e1.0)*p00)
    VL <- EL2 - sqr( EL)
  }

  if( debug) {
    mtrace( V.allX) # surely wanna
    0 # help debugging...
  }

  # If no Odis (e.g. with very few loci!) then no point in going further
  rho.hat <- if( C.hat > L-1) 100 else
      find.root( V.allX, target=emp.V.HSP, start=1/L, step=0.2/L,
          fdirection='decreasing', min.x=0, meioses=2)
  Vx <- V.allX( rho.hat, meioses=n.meio)

return( unlist( returnList( VUP=L*v0, V.HSP=emp.V.HSP, V0, Vx, C.hat, rho.hat, n.meio)))
}
