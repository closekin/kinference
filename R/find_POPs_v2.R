# see find_POPs for documentation
#' @export
#' @importFrom gbasics sqr
#' @importFrom atease @
#' @importFrom vecless := compile_vecless
#' @importFrom stats runif
#' @importFrom mvbutils cq %upto% %that.are.in% my.all.equal extract.named %without.name%
"find_POPs_v2" <-
function( snpg, subset1=1 %upto% nrow( snpg), subset2=subset1,
    alpha,
    # pOC_max,
    one_in_X_eta,
    rough_n_pairs_to_keep= NA,
    eta= NULL,
    keep_thresh= NULL,
    nq,
    quick=TRUE) {
###################
  # Sanity...
stopifnot( is.numeric( subset1) && is.numeric( subset2))
stopifnot( all( !duplicated( subset1)) && all( !duplicated( subset2)))
stopifnot( my.all.equal( subset1, subset2) || !length( intersect( subset1, subset2)))

  # Decide based #apparent exclusions of AA/BB form, using 4way genos, though
  # ... it's really AAO/BBO so not a true exclu but
  # ... close among pop-loci
  # Sticking with 4way genos so that genotyping errors are low

  # Instead of Nexclu, uses a wted sum of "exclus" in 4way genos to max expected diff between POP and UP
  # This version doesn't allow for geno errors, but does realize that AO/BO could happen
  # Doesn't bother with OO-AB

  define_genotypes()
  extract.named( snpg@locinfo[ cq( use6, PUP4, pbonzer)])
  p0 <- pbonzer[,'O'] + pbonzer[,'C']
  pA <- pbonzer[,'A']
  pB <- pbonzer[,'B']

  # "Exclusion" whenever AAO & BBO, but this *could* be AO/BO
  # I'm calling them "ex" for now anyway, hence pex_BLAH
  pex <- cbind(
      POP= 2 * p0 * pA * pB,
      UP= 2 * (2*pA*p0 + pA*pA) * (2*pB*p0 + pB*pB)
    )

  # Want to keep the POP and UP means as far apart as possible on the scale of SDs...
  # ... but, which SD? Make it alpha * SD[UP] + (1-alpha) * SD[POP]

  # Optimal wt would depend on p0 and to some extent on pA
  # wt should be 1 if p0==0 and 0 if p0==1
  delta <- pex[,'UP'] - pex[,'POP'] # mathematically I think this *can't* be -ve
  SD <- sqrt( pex * (1-pex))
  SD_combo <- alpha * SD[,'UP'] + (1-alpha) * SD[,'POP'] # %*% c( alpha, 1-alpha)
  V_combo <- sqr( SD_combo)
  ww <- delta / V_combo # considerable algebra appears to show this is optimal
  ww <- c( ww / sum( ww) ) # else get 1-col matrix
stopifnot( all( ww>0))

  pop_loci <- which( ww > 0) # all of them, for now
  # pop_loci <- which( pbonzer[,'O'] + pbonzer[,'C'] < pOC_max)

  temp_snpg <- snpg[ , pop_loci]
  recode4to6temp <- function( x) { x[ x=='AO'] <- AA; x[ x=='BO'] <- BB; x}
  temp_snpg <- recode4to6temp( temp_snpg) # (AA,AO) -> AA; (BB,BO) -> BB

  # Remove extranea
  attributes( temp_snpg) <- attributes( temp_snpg)[ 'dim']
  temp_snpg <- t( temp_snpg)

  # Distro of #excl loci for UPs
  opphetz <- c( 'AAO/BBO', 'BBO/AAO') %that.are.in% colnames( PUP4)
  pex_up <- PUP4[ pop_loci, opphetz]
  rr <- pex_up / (1-pex_up)
  log1m_pexup <- log1p( -pex_up)

  K <- function( tt) {
      # vecless2 notation could 1-step this via SUM
      KK[ l, j]:= log1m_pexup[ l] + log1p( rr[ l] * exp( tt[ j] * ww[ l]))
      K[ j] := KK[ +., j]
    return( c( K)) # without the c(), you get a scalar xtensor, and trouble...
    }

  dK <- function( tt) {
      retw[l,j] := rr[ l] * exp( tt[ j] * ww[ l])
      dKK[ l, j] := ww[ l] * retw[ l, j] / (1+retw[ l, j]) # guessing this is more accurate tahn 1-1/1+x
      dK[ j] := dKK[ +., j]
    return( c( dK))
    }

  ddK <- function( tt) {
     retw[l,j] := rr[ l] * exp( tt[ j] * ww[ l])
      ddKK[ l, j] := sqr( ww[ l]) * retw[ l, j] / sqr( 1+retw[ l, j]) # guessing this is more accurate tahn 1-1/1+x
      ddK[ j] := ddKK[ +., j]
   return( c( ddK))
   }

  n_loci <- ncol( snpg)
  n_sim_check <- 1000
  Ktest <- function( tt) {
    x <- matrix( runif( n_sim_check * n_loci) < pex_up, n_loci, n_sim_check)
    ewx <- exp( x*ww*tt)
    colSums( log( ewx))
  }

  if( quick) {
    K <- compile_vecless( K(0))
    dK <- compile_vecless( dK(0))
    ddK <- compile_vecless( ddK(0))
  }

  symmo <- my.all.equal( subset1, subset2)
  inv_CDF <- renorm_SPA_cumul( K, dK, ddK)$inv_CDF

  set_thresholds( keeping='lo') # cf 'hi' for HSPs

  # Prepare for diagnostics of #excl
  # prolly not needed but HSP C code already does it
  qq <- (2:nq-1)/nq
  pciles <- inv_CDF( qq) # inv_CDF_SPA2 may struggle with LOWER tail...

  # Trying special-cases here to minimize copying
  if( symmo){
    if( !my.all.equal( subset1, 1 %upto% ncol( temp_snpg))) {
      temp_snpg <- temp_snpg[, subset1]
    }
    result <- POP_wt_paircomps_lots(
        geno1= temp_snpg,
        geno2= temp_snpg,
        w= ww,
        symmo= TRUE,
        eta= eta,
        max_keep_wpsex= keep_thresh,
        bins= pciles,
        AAO= match( 'AA', snpg@diplos), # NB NB: AO has been recoded to AA
        BBO= match( 'BB', snpg@diplos)
      )
  } else { # different subsets
    result <- POP_wt_paircomps_lots(
        geno1= temp_snpg[ ,subset1],
        geno2= temp_snpg[ ,subset2],
        w= ww,
        symmo= FALSE,
        eta= eta,
        max_keep_wpsex= keep_thresh,
        bins= pciles,
        AAO= match( 'AA', snpg@diplos), # NB NB: AO has been recoded to AA
        BBO= match( 'BB', snpg@diplos)
      )
  }

  # bigs is a misnomer for POPs; smalls is more like it (since we want few exclusions)
  result$bigs <- with( result, data.frame( wpsex=big_wpsex, i=big_i, j=big_j))
  # Check nABOO, only for interesting pairs
  snpg_i <- snpg[ subset1[ result$bigs$i], pop_loci]
  snpg_j <- snpg[ subset2[ result$bigs$j], pop_loci]
  isABOO <- ((snpg_i==OO) & (snpg_j==AB)) + ((snpg_i==AB) & (snpg_j==OO))
  result$bigs$nABOO <- rowSums( isABOO)

  result <- result %without.name% cq( big_wpsex, big_i, big_j)
  result <- within( result, {
    bins <- pciles
    eta <- eta
    keep_thresh <- keep_thresh
    n_loci <- length( pop_loci)
    mean_theory <- dK( 0)
    var_theory <- ddK( 0)
  })
  result$call <- sys.call() # iffy within within

return( result)
}

