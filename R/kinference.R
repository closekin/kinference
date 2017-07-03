
# first two lines here are for Rcpp to work
#' @useDynLib kinference
#' @importFrom Rcpp evalCpp
#' @importFrom mvbutils %without.name% ?
"add_list_defaults" <- function( l, ...) {
###### Add defaults to list 'l' if not already in 'l'
  defaults <- list(...)
  l <- c( l, defaults %without.name% names( l))
return( l)
}

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

  result$bigs <- with( result, data.frame( PLOD=big_PLOD, i=big_i, j=big_j))
  result <- result %without.name% cq( big_PLOD, big_i, big_j)
  result$bins <- bins
  result$binprobs <- binprobs
  result$eta <- eta
  result$keep_thresh <- keep_thresh
  result$call <- sys.call()

return( result)
}




#' @importFrom atease @ @<-
#' @importFrom mvbutils %is.a%
"define" <- function( expr, expand_dim=FALSE) {
  expr <- substitute( expr)
  stopifnot( ((expr %is.a% '<-') || ((expr %is.a% 'call') && (expr[[1]]==as.name( '<<-')) ) ) && (expr[[2]] %is.a% 'name'))
  obj <- expr[[2]] # name/symbol
  name <- as.character( obj)

  eval.parent( expr)
  if( expand_dim && length( prefixdims)) { # known from environment
    prefixdims <- prefixdims # since these live in the parent, and substitute() won't find them
    dimorlen <- function( x) { if( is.null( d <- dim( x))) d <- length( x); d}
    eval.parent( substitute(
        obj <- structure( rep( obj, prod( prefixdims)), dim=c( prefixdims, dimorlen( obj)))
      ))
  }

  e <- environment( sys.function())
  eval.parent( substitute( {
    dimnames( obj) <- NULL # they just slow things down
    oldClass( obj) <- 'playback'
    obj@where <- e
    obj@whoami <- name
  }))
eval.parent( obj)
}


#' @importFrom mvbutils cq extract.named named mlocal
# @export
"define_genotypes" <- function( nlocal=sys.parent()) mlocal({
  ABCO <- named( cq( A, B, C, O))
  extract.named( ABCO) # A, B, C, and O

  genotypes <- cq( OO, AO, BO, AB, AA, BB, AAO, BBO, AC, BC, CO, CC, CCO)
  genotypes_ambig <- cq( OO, AB, AC, BC, AAO, BBO, CCO)
  genotypes4_ambig <- cq( OO, AB, AAO, BBO)
  genotypes6 <- cq( AA, AB, AO, BB, BO, OO)
  genotypes_C <- cq( AA, AB, AC, AO, BB, BC, BO, CC, CO, OO)

  for( ig in genotypes) {
    # assign( ig, structure( as.raw( match( ig, genotypes)), class='ABOSNP'))
    assign( ig, structure( ig, class='noquote')) # for nicer printing
  }
})



#' @importFrom mvbutils %is.not.a% %where%
#' @export
"get_chain" <- function( thing, seed) {
  if( thing %is.not.a% 'data.frame') {
    thing <- thing$bigs
  }
  extract.named( thing)

  oset <- integer()
  set <- seed
  while( length( set) != length( oset)) {
    newj <- j[ i %in% set]
    newi <- i[ j %in% set]
    oset <- set
    set <- unique( c( set, newi, newj))
  }

  thing %where% (i %in% set | j %in% set)
}

#' @importFrom stats pnorm dnorm qnorm
#' @importFrom gbasics logit inv.logit
#' @export
"inv_CDF_SPA2" <- function( p, K, dK, ddK, tol=formals( ridder)$tol) {

######## Invert L-R SPA approx to CDF on "s-scale"
######## Avoids "double iteration" of nonlinearity
######## Should work for vector p (the target) but only if your K etc do

  isqrt_2pi <- 1/sqrt(2*pi)
  x <- sqrt_ddK_s <- Leg_trans <- u <- w <- 0*p # vectorized

  CSPA <- function( s) {
      x <<- dK(s)
      Leg_trans <<- s*x - K(s)
      w <<- sign(s) * sqrt( 2*Leg_trans)
      sqrt_ddK_s <<- sqrt( ddK( s))
      u <<- s * sqrt_ddK_s
    return( pnorm( w) + dnorm( w) * (1/w - 1/u) - p_target)
    }

  p_target <- 0
  isK2 <- 1/sqrt(ddK(0))
  seps <- 0.001 * isK2 # fraction of 1 SD; (tol/2) tends to give numeric errors...
  p0 <- (CSPA( seps) + CSPA( -seps)) / 2 # CSPA(0)==NA--- avoid!
  is_lower <- p0 > p
  # K3 <- 3 * (sqr(u)-sqr(w)) / seps^3 # yeh not bad FWIW

  hi <- lo <- 0*p

  # Shifted start...
  q0 <- qnorm( inv.logit( logit( p) - logit( p0)))
  bingo <- q0==0.5 # this case won't converge
  mean_x <- dK( 0)

  # L-R does not work at x==E[X]..!
  if( any( bingo)) { # bingo!
    # This was in the original code, and needed vectorizing...
    # I think p[bingo] corresponds exactly (??) to mean(X)
    if( all( bingo)) {
return( 0*x + mean_x)
    }

    # change p for those cases to something that will converge and carry on
    # sub back correct x (ie mean) on exit
    p[ bingo] <- p[ !bingo][1]
    q0[ bingo] <- q0[ !bingo][1]
  }

  p_target <- p # so CSPA gives 0 at solution

  # Only risk I can see, is that if true s ~= 0, start may be on wrong side... and CSPA calcs go wrong.
  s0 <- 0.9 * q0 *isK2 # "zeroth-order" approx, times 0.9 for guess at lower bound
  while( any(
      bad <- !is.finite(
        C0 <- CSPA( s0)))) {
    s0[ bad] <- s0[ bad]/2
  }

  toobig <- C0 > 0
  multor <- 2 * xor( toobig, s0>0) - 1 # -1 or +1
  bracketed <- s0 != s0
  step <- 1.2 ^ multor
  snext <- s0

  repeat{
    snext[ !bracketed] <- s0[ !bracketed] * step[ !bracketed]
    bracketed[ !bracketed] <- xor( toobig, CSPA( snext) > 0)[ !bracketed]
    if( all( bracketed))
  break
    s0[ !bracketed] <- snext[ !bracketed]
  }

  hi <- pmax( s0, snext)
  lo <- pmin( s0, snext)

  s <- ridder( CSPA, lo, hi, tol=tol, skip_bounds=TRUE) # root finder
  CSPA( s)

  x[ bingo] <- mean_x # any that hit first time
return( x)
}


#' @importFrom atease @ @<-
#' @importFrom mvbutils do.on
#' @export
"pick_FSPs_from_HSPs" <- function( snpg, HSPs) {
  # For pairs already picked as HSPs, ie PLOD(HSP,UP) > eta: they might be FSPs. H
  # How would a FSP / HSP comparison look?
  # Don't need full pairwise screening for FSPs (do post hoc on a few hundred HSPs), hence all in R.

  # HSPs should be M*2 matrix of rows in snpg that are HSPs or FSPs

  # Transform to 4way genotypes
  # based on code in find_duplicates
  # careful, since "factor level" of AB and OO is different in 4way vs 6way
  sibg <- just_sibg <- snpg[ c( HSPs),]

  sibg@diplos <- genotypes4_ambig
  sibg[ just_sibg==AO] <- AAO
  sibg[ just_sibg==AA] <- AAO
  sibg[ just_sibg==BO] <- BBO
  sibg[ just_sibg==BB] <- BBO
  sibg[ just_sibg==OO] <- OO # need to do OO & AB too, since codes are different in 4way vs 6way
  sibg[ just_sibg==AB] <- AB


  extract.named( snpg@locinfo[ cq( PUP4, LOD4)])
  PHSP4 <- exp( LOD4) / PUP4 # Pr[gg|HSP] <- 0.5 * PUP4 + 0.5 * Pr[gg|kappa=1]
  P_k1 <- 2*PHSP4 - PUP4
  P_k1[ P_k1 < 0] <- 0 # rounding error

  samo <- do.on( strsplit( colnames( LOD4), '/'), .[1]==.[2])
  PFSP4 <- ...
  # Number of ID locis should be

  nsibs <- nrow( HSP)

  n_loci_id <- sibg[ 1:nsib,] == sibg[ nsib + 1:nsib,]
  }


#' @importFrom atease @ @<-
#' @importFrom gbasics make_genopairer sqr
"predict_hsp_util" <-
function( pIBD0, pIBD1, want_LOD_table=FALSE, k=0.5) {
  # This version ignores the possibility of errors involving AB or OO...
  # ... which should be pretty rare

  define_genotypes()
  nl <- nrow( pIBD1)
  Phsp <- pIBD1 * k + pIBD0 * (1-k)
  Pup <- pIBD0

  LOD <- log( Phsp / Pup)
  LOD[ Pup==0] <- 0 # if Pup=0 then p*log(p) = 0; only happens when r=0
  if( want_LOD_table) {
    # LOD is 3D: nloci * ng1 * ng2
    # gpLOD is 2D: nloci * n_genopairs
    # Need only certain "columns" of 2D-fied LOD
    mg <- make_genopairer( dimnames( pIBD0)[[2]])
    ngp <- max( mg)
    wanted <- match( 1:ngp, mg)

    gpLOD <- gpPUP <- matrix( 0, nl, ngp)
    LOD_as_2D <- matrix( LOD, nl, prod( dim( pIBD0)[-1]))
    gpLOD <- LOD_as_2D[ , wanted]
    PUP_as_2D <- matrix( Pup, nl, prod( dim( pIBD0)[-1]))
    gpPUP <- PUP_as_2D[ , wanted]

    # Off-diagnonals appear twice, and prob should be doubled...
    omg <- mg
    omg[ wanted] <- 0
    double_wanted <- 1:ngp %in% omg
    # wrong for some reason:    double_wanted <- wanted %in% mg[ duplicated( c( mg))]
    gpPUP[ ,double_wanted] <- gpPUP[,double_wanted] * 2

    dimnames( gpLOD) <- dimnames( gpPUP) <- list( dimnames( pIBD0)[[1]], mg@what)
    gpLOD@mg <- mg # why not
  }

  # EPLOD is sum( LOD * Pup) but we want to keep it by locus for now
  E.HSP[l] := LOD[l,i,j] %[i,j]% Phsp[l,i,j]
  E.UP[l] := LOD[l,i,j] %[i,j]% Pup[l,i,j]
  E2.UP[l] := (LOD*LOD)[l,i,j] %[i,j]% Pup[l,i,j]
  V.UP <- E2.UP - sqr( E.UP)
  Ediff <- E.HSP - E.UP

  # Standardized difference ie locus power: not so useful post hoc, but possibly interesting for 6 vs 4 comps
  sdiff <- (E.HSP - E.UP) / sqrt( V.UP)

#  Ediff <- unclass( Ediff)
#  V.UP <- unclass( V.UP)
#  sdiff <- unclass( sdiff)
#

  retval <- data.frame( Ediff, V.UP, sdiff)
  if( want_LOD_table) {
    retval@LOD <- gpLOD
    retval@PUP <- gpPUP
  }
return( retval)
}


#' @importFrom gbasics ridder integ
#' @importFrom stats splinefun
# @export
"renorm_SPA" <- function(K, dK, ddK, return_what=c( 'func', 'mulfuncby'),
                         tol=formals( ridder)$tol
  # , ... ; should really allow extra args to K & co, and build them in...
){
  isqrt_2pi <- 1/sqrt( 2*pi)
  absmax <- 10 / sqrt( ddK( 0)) # x between +/- 10SD of mean

  # mc <- as.list( match.call( expand.dots=FALSE)$...)
  K <- Vectorize( K) # names( formals( K)) %except% names( list( ...)) kinda thing
  dK <- Vectorize( dK)
  ddK <- Vectorize( ddK)

  itotto <- 1
  sfunc <- function( s, ddK_s=ddK(s)) {
      x <- dK( s)
      itotto * exp( K( s) - s * x) * sqrt( ddK_s) * isqrt_2pi
    }

  itotto <- 1 / integ( sfunc( x), -absmax, absmax)

  return_what <- match.arg( return_what)
  if( return_what=='mulfuncby')
return( itotto)

  # Otherwise we need the x-ready version, with (vectorized) root-finding
  xfunc <- function( x) {
    # Lower & upper bounds for s

    # Aim for exact, miss, then try to get equal dist the other side
    iddK0 <- 1/ddK( 0)
    s1 <- (x-dK(0)) * iddK0

    dK1 <- dK( s1)
    # Paranoia: might be perfect!
    bingo <- dK1==x
    if( any( bingo)) { # otherwise, if K()==Vectorize(...), it fucks up and turns it all into a list FFS
      s1[ bingo] <- 1.99 * s1[ bingo]
      dK1[ bingo] <- dK( s1[ bingo])
    }

    ddK1 <- ddK( s1)
    s2 <- s1 - 2*(dK1-x) / ddK1
    dK2 <- dK( s2)
    while( any( same_sign <- (dK1-x)*(dK2-x)>1) ) { # unlikely; try a bit further
      s2 <- s2 - same_sign * (dK1-x)/ddK1
      dK2[ same_sign] <- dK( s2[ same_sign])
    }

    # 'ridder' wants a vector func with no args
    dK_min_x <- function( s) dK(s)-x
    s <- ridder( dK_min_x, pmin( s1, s2), pmax( s1, s2), tol=tol) # root finder
    ddK_s <- ddK( s)

    # Undo the reordering
  return( sfunc( s, ddK_s) / ddK_s)
  }

return( xfunc)
}

#' @importFrom mvbutils returnList
# @export
"renorm_SPA_cumul" <- function( K, dK, ddK, sd_half_range=10, n_pts=2001) {
  x <- 0
  SPA_s_dxds <- function( s) {
    x <<- dK( s)
    # SPA is (K(x) - s*x) / sqrt( 2*pi*ddK( s))
    # but dx/ds = ddK( s)
    exp( K(s) - s*x) * sqrt( ddK( s) / (2*pi))
  }

  sd <- sqrt( ddK( 0))
  # To norm, we'd integ X over say +/- 10 sd
  # but s ~= (x-mu) / (sd)^2 hence s-range below

  spoints <- seq( -sd_half_range, sd_half_range, length=n_pts) / sd
  pdf_s <- SPA_s_dxds( spoints) * diff( spoints[1:2]) # diff() gets the integral about right
  C <- 1/sum( pdf_s) # should be close to 1
  cdf <- cumsum( pdf_s) * C # nicely renormalized

  inv_CDF_bod <- splinefun( cdf, x,  method='hyman')
  inv_CDF <- function( p) inv_CDF_bod( p) # sensibler arg name, and no deriv arg
  CDF <- splinefun( x, cdf, method='hyman')
returnList( CDF, inv_CDF)
}

#' @importFrom mvbutils cq mlocal %is.a%
"set_thresholds" <-
function( keeping, nlocal=sys.parent()) mlocal({
stopifnot( keeping %in% cq( hi, lo))

  symmo <- my.all.equal( subset1, subset2)
  if( is.null( eta) || is.null( keep_thresh)) {
    probinverts <- numeric()
    if( is.null( eta)) {
      probinverts <- 1/one_in_X_eta
    }
    if( is.null( keep_thresh)) {
      probinverts <- c( probinverts,
          rough_n_pairs_to_keep / (length( subset1) * length( subset2) / (1+symmo)))
    }

    if( keeping == 'hi') {
      probinverts <- 1-probinverts
    }

    XX <- try( inv_CDF_SPA2( probinverts, K, dK, ddK))
    if( XX %is.a% 'try-error') {
  warning( "Couldn't set thresholds via Lugannini-Rice SPA (will use alternative); ' %&%
      'probably too extreme for this distro too.")
      # Use renormed sum-of-pdf:
      XX <- inv_CDF( probinverts)
    }

    if( is.null( eta)) {
      eta <- XX[1]
      XX <- XX[-1]
    }
    if( is.null( keep_thresh)) {
      keep_thresh <- ( if( keeping == 'lo') max else min)( XX[1], eta)
    }
  }
})

## TODO
# - should this exist?
#' @importFrom atease @ @<-
#' @importFrom gbasics rsample
#' @importFrom stats var
#' @importFrom mvbutils scatn
"simtest_Kstuff" <- function( ck, n, nq=20) {
  # ck needs locinfo$LOD
  extract.named( ck@Kenv) # 4ways pretending to be 6ways
  mg <- LOD@mg
  n_loci <- nrow( PUP)

  # Direct simulation of genotypes is tricky because of errors
  # ... though for completeness SHOULD really try that here
  # Instead, use PUP table

  g12_code <- matrix( 0L, n, n_loci)
  for( il in 1:n_loci) {
    this_g12 <- rsample( n, colnames( PUP), prob=PUP[il,], replace=TRUE)
    g1 <- substring( this_g12, 1, 2)
    g2 <- substring( this_g12, 4, 5)
    g12_code[ ,il] <- mg[ cbind( g1, g2)]
  }

  # NYI in vecless:
  # LOD_obs[ l, i] := LOD[ l, g12[ l, i]]
  # PLOD[ i] := LOD_obs[ +., i]

  LOD_obs <- matrix( 0, n, n_loci)
  LOD_obs[] <- LOD[ cbind( rep( 1:n_loci, each=n), c( g12_code)) ]
  PLOD <- rowSums( LOD_obs)

  scatn( 'Emp mean %6.4f SPA %6.4f', mean( PLOD), dK( 0))
  scatn( 'Emp var %6.4f SPA %6.4f', var( PLOD), ddK( 0))

  # %les:
  qq <- (2:nq-1)/nq
  pciles <- inv_CDF( qq)

  bin <- 1+findInterval( PLOD, pciles)
  counts <- tabulate( bin, nbins=nq)
  scatn( 'Expect %5.1f PLODs in each of %i percentile-bins; got this:', n/nq, nq)
  print( counts)

return( invisible( PLOD))
}


