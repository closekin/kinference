
#' @export
"add_pairprob_error" <-
function( nlocal=sys.parent()) mlocal({
## Needs pp_true and snerr

  g_err <- named( genotypes_C) ? 0 # recorded because it makes perr easy to define
  perr <- 0*nchar( g_err) ? 1 # named

  g_1 <- substring( genotypes_C, 1, 1)
  g_2 <- substring( genotypes_C, 2, 2)
  AA_BB <- g_1==g_2 & (g_1==A | g_1==B)
  g_err[ AA_BB] <- g_1[ AA_BB] %&% 'O'
  perr[ AA_BB] <- snerr[ sprintf( '%1$s%1$s2%1$sO', g_1[ AA_BB])] ? 0 # XX2XO

  AO_BO <- (g_1==A | g_1==B) & g_2==O
  g_err[ AO_BO] <- g_1[ AO_BO] %&% g_1[ AO_BO]
  perr[ AO_BO] <- snerr[ sprintf( '%1$sO2%1$s%1$s', g_1[ AO_BO])] ? 0 # XO2XX

  # Original version: not sure it's selfrecordable
#  make_pp_err <- function( g1, g2) {
#      g1_err <- g_err[ g1]
#      g2_err <- g_err[ g2]
#      pp_true[ cbind( g1, g2) ] * (1-perr[ g1]) * (1-perr[ g2]) +
#      pp_true[ cbind( g1_err, g2) ] * perr[ g1_err] * (1-perr[ g2]) +
#      pp_true[ cbind( g1, g2_err) ] * (1-perr[ g1]) * perr[ g2_err] +
#      pp_true[ cbind( g1_err, g2_err) ] * perr[ g1_err] * perr[ g2_err]
#    }
#
  # pp_err[] <- outer( named( genotypes_C), named( genotypes_C), make_pp_err)

  # Need 2-stage setup here so selfrec works: "pure" assignment can't be overloaded in R
  pp_err <- matrix( -1, length( genotypes_C), length( genotypes_C), dimnames=rep( list( genotypes_C), 2)) ? 1

  gg <- list() ? 0 # this is a trick to allow with(gg,...) in playback, since the "real" gg isn't kept

  gg <- expand.grid( g1=genotypes_C, g2=genotypes_C, stringsAsFactors=FALSE)
  gg <- within( gg, {
     g1_err <- g_err[ g1]
     g2_err <- g_err[ g2]
  })
  extract.named( gg)

  pp_err[ cbind( g1, g2)] <-
      pp_true[ cbind( g1, g2) ] * (1-perr[ g1]) * (1-perr[ g2]) +
      pp_true[ cbind( g1_err, g2) ] * perr[ g1_err] * (1-perr[ g2]) +
      pp_true[ cbind( g1, g2_err) ] * (1-perr[ g1]) * perr[ g2_err] +
      pp_true[ cbind( g1_err, g2_err) ] * perr[ g1_err] * perr[ g2_err]        ? 0


  # Merge C with O. Has to be mlocal() to avoid <<- which buggers playback. Things before 'nlocal' are args; things after are temporaries
  add_up <- function( i, result, ..., nlocal=sys.parent(), dimor, other)  mlocal({
      result <- match( result, genotypes_C)
      dimor <- slice.index( pp_err, i)

      for( other in FOR( list( ...), match( ., genotypes_C))) {
        pp_err[ dimor==result] <- pp_err[ dimor==result] + pp_err[ dimor==other] ? 0
        pp_err[ dimor==other] <- 0 # avoid double-use ? 0
      }
    })

  for( i in 1:2) {
    add_up( i, OO, CC, CO)
    add_up( i, AO, AC)
    add_up( i, BO, BC)
  }

# IE:
# pp_err2[,OO] <- pp_err2[,CC] + pp_err2[,CO] + pp_err2[,OO]
# pp_err2[,CC] <- pp_err2[,CO] <- 0
# pp_err2[AO,] <- pp_err2[AC,] + pp_err2[AO,]
# pp_err2[AC,] <- 0

  pp6_err <- pp_err[ genotypes6, genotypes6] ? 0
})


# @export
#' @importFrom vecless set_recording
#' @importFrom mvbutils %&%
"calc_g6probs_IBD0_scalar" <-
function( P, snerr, record=FALSE) {
## SCALAR-ONLY VERSION... this is hard enough!
## Though can be called with 1-row matrix args, eg with( x@locinfo[1,], calc_g6probs_IBD1( pbonzer, snerr))

  # snerr = P( misclassifying true XX as XO, and vice versa)

  set_recording( cq( P, snerr, pp_true, pr2, pp_err, perr, pp6_err), record)
  define_genotypes()           ? 0

  P <- drop( P) # for scalar version
  P <- P              ? 0
stopifnot( my.all.equal( names( P), names( ABCO)))

  snerr <- drop( snerr) # for scalar version

  snerr <- snerr            ? 0
  pp_true <- matrix( 0, length( genotypes_C), length( genotypes_C),
                    dimnames=list( genotypes_C, genotypes_C))         ? 1


  g_1 <- substring( genotypes_C, 1, 1)
  g_2 <- substring( genotypes_C, 2, 2)
  pr2 <-  nchar( named( genotypes_C))    ? 1 # named

  is_het <- g_1 != g_2
  pr2[] <- P[ g_1] * P[ g_2]             ? 0
  pr2[ is_het] <- 2 * pr2[ is_het]       ? 0

  pp_true <- matrix( 0, length( genotypes_C), length( genotypes_C),
                    dimnames=rep( list( genotypes_C), 2))        ? 1

  extract.named( expand.grid( gp1=genotypes_C, gp2=genotypes_C,
                             stringsAsFactors=FALSE))
  pp_true[ cbind( gp1, gp2)] <- pr2[ gp1] * pr2[ gp2]                       ? 0

  # NB that for this UP case, XX/XO errors shouldn't change the overall probs because the cutoffs are chosen to do exactly that!
  add_pairprob_error()

return( pp6_err)
}


# @export
#' @importFrom vecless set_recording
"calc_g6probs_IBD1_scalar" <-
function( P, snerr, record=FALSE) {
## SCALAR-ONLY VERSION... this is hard enough!
## Though can be called with 1-row matrix args, eg with( x@locinfo[1,], calc_g6probs_IBD1( pbonzer, snerr))

  set_recording( cq( P, snerr, pp_true, pp3, pp_err, perr, pp6_err), record)
  define_genotypes()           ? 0

  P <- drop( P) # for scalar version
  P <- P               ? 0
stopifnot( my.all.equal( names( P), names( ABCO)))

  snerr <- drop( snerr) # for scalar version
  snerr <- snerr        ? 0
  pp_true <- matrix( 0, length( genotypes_C), length( genotypes_C),
                    dimnames=list( genotypes_C, genotypes_C))          ? 1

  # Note: AB/AB can share either A or B... so have to accumulate

  # always write as S,U1 and S,U2 so probs are unambig
  # then accumulate into sorted versions since *several* S,U1,U2 could contribute to one XY/ZW
  # ... or *none*, of course, if there is no shared allele in XY/ZW

  # 3-col DF of (shared, 1st unshared, 2nd unshared) alleles in the pair
  su1u2 <- as.matrix( expand.grid( ABCO, ABCO, ABCO))
  implied <- su1u2[ ,c(1,2,1,3)]
  swap12 <- implied[,2] < implied[,1]
  implied[ swap12, 1:2] <- implied[ swap12, 2:1]
  swap34 <- implied[,4] < implied[,3]
  implied[ swap34, 3:4] <- implied[ swap34, 4:3]
  cimplied <- cbind( implied[,1] %&% implied[,2], implied[,3] %&% implied[,4])

  pp3 <- P[ su1u2[,1]] * P[ su1u2[,2]] * P[ su1u2[,3]]             ? 0

  # Each row can contribute to at most 2 implied XY/ZW combinations
  # eg AB/AB <- ABB or BAA
  # AB/CD <- nothing
  # AB/AC <- ABC

  implstr <- cimplied[,1] %&% cimplied[,2]
  has_impl2 <- which( duplicated( implstr))

  first_impl <- match( implstr[ -has_impl2], implstr)
  second_impl <- match( implstr[has_impl2], implstr, 0)
  cimpl2 <- cimplied[ second_impl,]

  pp_true[ cimplied[ first_impl,]] <- pp3[ -has_impl2]            ? 0
  pp_true[ cimpl2] <- pp_true[ cimpl2] + pp3[ has_impl2]          ? 0

  # Allow for XX <-> XO errors--- hopefully the only ones! (Watch out for scaffoldy version)
  # Assumes errors INDEPENDENT even tho there could be heritability (eg weak grabbing of mutated primer)
  # CC <-> CO is ignored in this version, since C gets scored as O

  add_pairprob_error()

return( pp6_err)
}


