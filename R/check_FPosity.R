## TODO
# - CDF calls commented out
# - to export?
#' @importFrom mvbutils cq %without.name% %&%
#' @importFrom gbasics sqr
#' @importFrom atease @ @<-
# @export
"check_FPosity" <- function( snpg, nsim=0){
## snpg should have been thru 'prepare_PLOD_SPA' so it has @PPS
stopifnot( 'Kenv' %in% names( attributes( snpg)))

  og <- options( vecless.print=FALSE)
  on.exit( options( og))

  define_genotypes()
  for( iwhat in cq( LOD, PUP, PUPLOD, PUPLOD2)) {
    assign( 'O' %&% iwhat, snpg@Kenv[[ iwhat]])
  }
  mg <- OLOD@mg

  use6 <- snpg@locinfo$use6
  use4 <- !use6
  temp_snpg <- snpg
  recode4to6temp <- function( x) { x[ x=='AO'] <- AA; x[ x=='BO'] <- BB; x}
  temp_snpg[ , use4] <- recode4to6temp( snpg[, use4]) # (AA,AO) -> AA; (BB,BO) -> BB

  # For 4way loci, temporarily treat XO as XX...
  # ... have already adjusted the LOD entries so that new_LOD6( XX/..) <- LOD4( XXO/..)
  # ... use the LOD that's in Kenv, where SPA is calculated

  # Genofreqs could've/should've be done at the start in hsp_power, but here will do

  # LOD and PUP are stored in compacted 2D form to save space... need to fix that
  # Can't quite do this with vecless!
  OLOD[ is.na( OLOD)] <- 0 # set to NA for 4way loci
  n_loci <- nrow( OPUP)
  NPUP <- PUP <- LOD <- array( 0, c( n_loci, 6, 6))
  for( ig in 1:6) {
    gjseq <- mg[ , ig]
    XXi[ l, gj] := OPUP[ l, gj=gjseq] # Shouldn't work with new vecless syntax... but does !?
    PUP[ l, {ig}, gj] := XXi[ l, gj]
    NPUP[,ig,] <- OPUP[ , mg[,ig]]

    XXi[ l, gj] := OLOD[ l, gj=gjseq]
    LOD[ l, {ig}, gj] := XXi[ l, gj]
  }

  Pg[ l, gi] := sqrt( PUP[ l, gi, gi])
  # PHSP[l,gj,gi] := exp( LOD[ l, gj, gi]) * Pg[ l, gj] * Pg[ l, gi]
  # Pg2_g1_H[ l, gj, gi] := PHSP[ l, gj, gi] / Pg[ l, gi]

  e_CLOD[ l, gi] := LOD[ l, gi, gj] %[gj]% Pg[ l, gj]  # since gj indept gi
  e2_CLOD[ l, gi] := sqr( LOD[ l, gi, gj]) %[gj]% Pg[ l, gj]
  # e_CLOD_HSP[ l, gi] := LOD[ l, gi, gj] %[gj]% Pg2_g1_H[ l, gj, gi]
  # returnList( e_CLOD, e2_CLOD, e_CLOD_HSP)

  # Int version needed for vecless lookups
  geno <- as.integer( c( snpg))
  dim( geno) <- dim( snpg)

  # Ugly way to do lookups
  # my_e_CLOD[ i] := SUM_ %[l]% e_CLOD[ l, geno[ i, l]]
  # my_e2_CLOD[ i]:= SUM_ %[l]% e2_CLOD[ l, geno[ i, l]]
  g1seq <- seq_along( genotypes6)
  my_e_CLOD[ l, i] := e_CLOD[ l, g1] %[g1]% (geno[i,l] == g1seq[ g1])
  my_e2_CLOD[ l, i]:= e2_CLOD[l,g1] %[g1]% (geno[i,l] == g1seq[ g1])
  my_v_CLOD[ l, i]:= my_e2_CLOD[ l, i] - sqr( my_e_CLOD[ l, i])
  my_rat_CLOD[ i]:= (SUM_ %[l]% my_e_CLOD[ l, i] ) / sqrt( SUM_ %[l]% my_v_CLOD[ l, i])

  if( nsim) {
    gsim <- matrix( 0L, nsim, n_loci)
    for( il in seq_len( n_loci)) {
      gsim[,il] <- rsample( nsim, 1:6, prob=Pg[il,], replace=TRUE)
    }

    sim_e_CLOD[ l, i] := e_CLOD[ l, g1] %[g1]% (gsim[i,l] == g1seq[ g1])
    sim_e2_CLOD[ l, i]:= e2_CLOD[l,g1] %[g1]% (gsim[i,l] == g1seq[ g1])
    sim_v_CLOD[ l, i]:= sim_e2_CLOD[ l, i] - sqr( sim_e_CLOD[ l, i])
    sim_rat_CLOD[ i]:= (SUM_ %[l]% sim_e_CLOD[ l, i] ) / sqrt( SUM_ %[l]% sim_v_CLOD[ l, i])
  }


  # Vectorized individual KGFs, for each sample (columns) and numerous t-values (row)
  K <- function( tt) {
    # need to lookupize this until vecless 2.0 is out...
    ETT[ it, i, l, g] := exp( tt[ it] * LOD[ l, geno[ i, l], g])
    log_S[ it, i, l] := log( Pg[ l, g] %[g]% ETT[ it, i, l, g])
    K[ it, i] := log_S[ it, i, +.]
  return( K)
  }

  dK <- function( tt) {
    ETT[ it, l, g12] := exp( tt[ it] * LODOK[ l, g12])
    S[ it, l] := PUP[ l, g12] %[g12]% ETT[ it, l, g12]
    SL[ it, l] := PUPLOD[ l, g12] %[g12]% ETT[ it, l, g12]
    rowSums( SL/S)
  }

  ddK <- function( tt) {
    ETT[ it, l, g12] := exp( tt[ it] * LODOK[ l, g12])
    S[ it, l] := PUP[ l, g12] %[g12]% ETT[ it, l, g12]
    SL[ it, l] := PUPLOD[ l, g12] %[g12]% ETT[ it, l, g12]
    SLL[ it, l] := PUPLOD2[ l, g12] %[g12]% ETT[ it, l, g12]
    rowSums( (SLL/S-gbasics::sqr( SL/S)))
  }



stop()

  # Remove extranea
  attributes( temp_snpg) <- attributes( temp_snpg)[ 'dim']
  temp_snpg <- t( temp_snpg)

  if( is.null( bins)) {
    qq <- (2:nq-1)/nq
#CDF COMMENTED OUT    bins <- inv_CDF( qq)
  }
#CDF COMMENTED OUT  binprobs <- CDF( bins)

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

  result <- with( result, data.frame( PLOD=big_PLOD, i=big_i, j=big_j))
  result <- result %without.name% cq( big_PLOD, big_i, big_j)

  result@bins <- bins
  result@binprobs <- binprobs
  result@eta <- eta
  result@keep_thresh <- keep_thresh
  result@call <- sys.call()

return( result)
}
