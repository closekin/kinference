# Calculate some log-odds

calculate_IBD <- function(lociar){

  define_genotypes()
  li <- lociar@locinfo
  li1 <- li[1,]

  temp0 <- with( li1, calc_g6probs_IBD0_scalar( pbonzer, snerr, record=TRUE))
  cg6p0 <- make_playback( calc_g6probs_IBD0_scalar, temp0)

  temp1 <- with( li1, calc_g6probs_IBD1_scalar( pbonzer, snerr, record=TRUE))
  cg6p1 <- make_playback( calc_g6probs_IBD1_scalar, temp1)

  pIBD0 <- with( li, cg6p0( pbonzer, snerr))
  pIBD1 <- with( li, cg6p1( pbonzer, snerr))

  return(list(pIBD0 = pIBD0,
              pIBD1 = pIBD1))
}



calculate_LOD_HSP <- function(lociar, k=0.5){

  LODs <- calculate_IBD(lociar)

  pIBD0 <- LODs$pIBD0
  pIBD1 <- LODs$pIBD1

  nl <- nrow( pIBD1)
  Phsp <- pIBD1 * k + pIBD0 * (1-k)
  Pup <- pIBD0

  LOD <- log( Phsp / Pup)
  LOD[ Pup==0] <- 0 # if Pup=0 then p*log(p) = 0; only happens when r=0

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

  # EPLOD is sum( LOD * Pup) but we want to keep it by locus for now
  # expected value of the PLOD for HSP (mean of distn)
  E.HSP[l] := LOD[l,i,j] %[i,j]% Phsp[l,i,j]
  # expected value of the PLOD for UP (mean of distn)
  E.UP[l] := LOD[l,i,j] %[i,j]% Pup[l,i,j]
  E2.UP[l] := (LOD*LOD)[l,i,j] %[i,j]% Pup[l,i,j]
  V.UP <- E2.UP - sqr( E.UP)
  Ediff <- E.HSP - E.UP

  # Standardized difference ie locus power: not so useful post hoc,
  #  but possibly interesting for 6 vs 4 comps
  sdiff <- (E.HSP - E.UP) / sqrt( V.UP)

  retval <- data.frame( Ediff, V.UP, sdiff)
  retval@LOD <- gpLOD
  retval@PUP <- gpPUP

return( retval)

}
