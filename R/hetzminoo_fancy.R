#' Heterozygotes minus "OO" checking
#' @export
# This test looks at whether the allele frequencies in a given fish seem right, or if there are discrepancies due to (i) degraded DNA or (ii) sample contamination.
#' @param snpg object of type \code{\link[gbasics]{snpgeno}}
"hetzminoo_fancy" <-
function( snpg, target=c( 'rich', 'poor'), hist_pars=list(), multhresh=1) {
###################
  define_genotypes()
  extract.named( snpg@locinfo[ cq( use6, PUP4, pbonzer)])
  p0 <- pbonzer[,'O'] + pbonzer[,'C']
  pA <- pbonzer[,'A']
  pB <- pbonzer[,'B']

  v <- 2*pA*pB + sqr( p0) - sqr( 2*pA*pB-sqr( p0))
  target <- match.arg( target)
  edash <- if( target=='rich') {
      (2*pA*p0+sqr(pA)) * (1-sqr(1-pB)) + (2*pB*p0+sqr( pB))*(1-sqr(1-pB)) + sqr( p0) * (1-sqr( p0))
    } else {
      edash <- 2*(pA*pB + pA*p0 + pB*p0) # minus sign
    }

  ww <- edash / v
  ww <- c( ww / sum( ww) ) # else get 1-col matrix
stopifnot( all( ww>0))

  use_loci <- which( ww > 0) # all of them, for now
  temp_snpg <- snpg[ , use_loci]
  recode4to6temp <- function( x) { x[ x=='AO'] <- AA; x[ x=='BO'] <- BB; x}
  temp_snpg <- recode4to6temp( temp_snpg) # (AA,AO) -> AA; (BB,BO) -> BB

  delta <- (temp_snpg==AB) - (temp_snpg==OO)
  whmo[ f]:= delta[ f, l] %[l]% ww[ l]

  # Null distro: P[S==1] = pAB; P[S==-1] = pOO; P[S==0] = 1-pAB - pOO
  # E[ exp( tt*S)] = (e(tt)+1) * pAB + 1 + (e(-tt) -1)* pOO

  pAB <- 2*pA*pB
  pOO <- sqr( p0)
  four_pab_poo <- 4*pAB*pOO
  compaboo <- 1 - pAB - pOO

  # now setup the functions for the SPA
  # use the above as shortcuts

  K <- function( tt) {
      # here a *1 inside exp omitted
      etwab[ l, j] := pAB[l] * exp( tt[j] * ww[ l])
      # here a *-1 inside exp omitted
      etwoo[ l, j] := pOO[l] * exp( -tt[j] * ww[ l])
      # compaboo computed above, P(not AB & not OO)*exp(0)
      KK[ j]:= SUM_ %[l]% log( compaboo[l] + etwab[ l, j] + etwoo[ l, j])
    return( c( KK)) # without the c(), you get a scalar xtensor, and trouble...
    }

  # derivative of K
  dK <- function( tt) {
      etwab[ l, j] := pAB[l] * exp( tt[j] * ww[ l])
      etwoo[ l, j] := pOO[l] * exp( -tt[j] * ww[ l])
      denom[ l, j] := compaboo[l] + etwab[l,j] + etwoo[l,j]
      dKK[ j] := ww[l] %[l]% ((etwab[l,j]-etwoo[l,j])/denom[l,j])
    return( c( dKK))
    }

  # 2nd derivative of K
  ddK <- function( tt) {
      etwab[ l, j] := pAB[l] * exp( tt[j] * ww[ l])
      etwoo[ l, j] := pOO[l] * exp( -tt[j] * ww[ l])
      denom[ l, j] := compaboo[l] + etwab[l,j] + etwoo[l,j]
      ddKK[ j] := sqr( ww[l]) %[l]% ( (compaboo[l] * (etwab[ l, j] + etwoo[l,j]) + four_pab_poo[ l]) / sqr( denom[l, j]) )
    return( c( ddKK))
    }

  K <- compile_vecless( K(0))
  dK <- compile_vecless( dK(0))
  ddK <- compile_vecless( ddK(0))

  dens_SPA <- renorm_SPA( K, dK, ddK, 'func')

  # optional graphics and/or user-specified outputs
  switch( mode( hist_pars),
    list = {
        hist_pars <- add_list_defaults( hist_pars,
            main=sprintf( '%s: multhresh=%5.2f', target, multhresh),
            xlim= range( whmo), # so cutoff lines show
            xlab='', nclass=50)
        lv <- do.call( 'hist', c( list( x=whmo), hist_pars))
        with( lv, lines( mids, diff( breaks) * dens_SPA( mids) * sum( counts), col='green'))
        # abline( v=ncuts, col='red')
      },
    expression = eval( badhetz_hist_pars),
    NULL = NULL
  )

return( c( whmo))
}
