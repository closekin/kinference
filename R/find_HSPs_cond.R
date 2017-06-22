## TODO
# - unfinished
# - should this exist?
#' @importFrom atease @ @<-
#' @importFrom mvbutils cq %without.name% returnList
#' @importFrom gbasics sqr
"find_HSPs_cond" <- function(snpg, subset1=1 %upto% nrow(snpg), subset2=subset1,
    one_in_X_eta,
    rough_n_pairs_to_keep,
    eta= NULL,
    keep_thresh= NULL,
    nq= 50,
    bins= NULL) {
## snpg should have been thru 'prepare_PLOD_SPA' so it has @PPS
stopifnot( 'Kenv' %in% names( attributes( snpg)))

  # Sanity...
stopifnot( is.numeric( subset1) && is.numeric( subset2))
stopifnot( all( !duplicated( subset1)) && all( !duplicated( subset2)))
stopifnot( my.all.equal( subset1, subset2) || !length( intersect( subset1, subset2)))

  og <- options( vecless.print=FALSE)
  on.exit( options( og))

  # Here I'm using L-R tail approx SPA for CDF
  # ... although Kenv$inv_CDF is likely more accurate for "moderate" tails but I don't quite trust it in the extremes
  # ... actually they are pretty similar
  # ... Possibly, Kenv$inv_CDF should check if arg exceeds the range it was fitted to, and if so call
  # ... inv_CDF_SPA2() instead
  # ... but the range used in fitting is very goddamn wide (say +/- 10 SD) !

  define_genotypes()

  for( iwhat in cq( K, dK, ddK, inv_CDF)) {
    assign( iwhat, snpg@Kenv[[ iwhat]])
  }
  set_thresholds( keeping='hi')

  # For 4way loci, temporarily treat XO as XX...
  # ... have already adjusted the LOD entries so that new_LOD6( XX/..) <- LOD4( XXO/..)
  # ... use the LOD that's in Kenv, where SPA is calculated

  make_CLOD <- function( LOD, PUP) {
    # Could've/should've be done at the start in hsp_power, but here will do
    Pg[ l, gi] := sqrt( PUP[ l, gi, gi])
    PHSP[l,gj,gi] := exp( LOD[ l, gj, gi]) * Pg[ l, gj] * Pg[ l, gi]
    Pg2_g1_H[ l, gj, gi] := PHSP[ l, gj, gi] / Pg[ l, gi]

    e_CLOD[ l, gi] := LOD[ l, gi, gj] %[gj]% Pg[ l, gj]  # since gj indept gi
    e2_CLOD[ l, gi] := sqr( LOD[ l, gi, gj]) %[gj]% Pg[ l, gj]
    e_CLOD_HSP[ l, gi] := LOD[ l, gi, gj] %[gj]% Pg2_g1_H[ l, gj, gi]
  returnList( e_CLOD, e2_CLOD, e_CLOD_HSP)
  }

  extract.named( snpg@locinfo[ cq( use6, LOD6, LOD4)])
  use4 <- !use6
  temp_snpg <- snpg
  recode4to6temp <- function( x) { x[ x=='AO'] <- AA; x[ x=='BO'] <- BB; x}
  temp_snpg[ , use4] <- recode4to6temp( snpg[, use4]) # (AA,AO) -> AA; (BB,BO) -> BB
  temp_LOD <- LOD # from Kenv; already done in prepare_PLOD_SPA

  # Remove extranea
  attributes( temp_snpg) <- attributes( temp_snpg)[ 'dim']
  temp_snpg <- t( temp_snpg)

  if( is.null( bins)) {
    qq <- (2:nq-1)/nq
    bins <- inv_CDF( qq)
  }
  binprobs <- CDF( bins)

  mean_theory <- dK( 0)
  var_theory <- ddK( 0)

  # Trying special-cases here to minimize copying
  if( symmo) {
    if( !my.all.equal( subset1, 1 %upto% ncol( temp_snpg))) {
      temp_snpg <- temp_snpg[, subset1]
    }

#    result <- HSP_cond_paircomps_lots(
#      vec_LOD= LOD,
#      geno1= temp_snpg,
#      geno2= temp_snpg,
#      e_CLOD= e_CLOD,
#      e2_CLOD= e2_CLOD,
#      e_CLOD_HSP= e_CLOD_HSP,
#      e_typical_PLOD= mean_theory,
#      v_typical_PLOD= var_theory,
#      symmo= TRUE,
#      eta= eta,
#      min_keep_PLOD= keep_thresh,
#      bins= bins)
#  } else { # different subsets
#stop( "Fix the non-symm code, bozo...")
#    result <- HSP_cond_paircomps_lots( this+will+fail,
#        pair_geno= temp_LOD@mg,
#        LOD= t( temp_LOD),
#        geno1= temp_snpg[ , subset1],
#        geno2= temp_snpg[ , subset2],
#        symmo= FALSE,
#        eta= eta,
#        min_keep_PLOD= keep_thresh,
#        bins= bins
#      )
  }

  result$bigs <- with( result, data.frame( PLOD=big_PLOD, i=big_i, j=big_j))
  result <- result %without.name% cq( big_PLOD, big_i, big_j)
  result$bins <- bins
  result$binprobs <- binprobs
  result$eta <- eta
  result$keep_thresh <- keep_thresh
  result$call <- sys.call()

return( result)
}


