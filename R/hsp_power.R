#' @importFrom atease @ @<-
#' @importFrom vecless make_playback
# @export
"hsp_power" <- function( lociar,
    want_LOD_table, # T/F
    k # 0.5 for HSPs
){
############
  define_genotypes()
  li <- lociar@locinfo
  li1 <- li[1,]

  temp0 <- with( li1, calc_g6probs_IBD0_scalar( pbonzer, snerr, record=TRUE))
  cg6p0 <- make_playback( calc_g6probs_IBD0_scalar, temp0)

  temp1 <- with( li1, calc_g6probs_IBD1_scalar( pbonzer, snerr, record=TRUE))
  cg6p1 <- make_playback( calc_g6probs_IBD1_scalar, temp1)

  g6p0 <- with( li, cg6p0( pbonzer, snerr))
  g6p1 <- with( li, cg6p1( pbonzer, snerr))

  s6 <- predict_hsp_util( g6p0, g6p1, want_LOD_table, k=k)

  # For the 4-ways, must condense g6p's

  map6to4 <- matrix( 0, 6, 4, dimnames=list( genotypes6, genotypes4_ambig))
  # AB & OO are OK; AAO should receive both AA and AO; etc
  mm <- match( genotypes6, substring( genotypes4_ambig, 1, 2), 0) # the "AA" bit of "AAO"...
  yup <- cbind( which(mm>0), mm[ mm>0])
  map6to4[ yup] <- 1
  mm <- match( genotypes6, substring( genotypes4_ambig, 2, 3), 0) # ... and the "AO" bit
  yup <- cbind( which(mm>0), mm[ mm>0])
  map6to4[ yup] <- 1

  # Really want g4p0[l,i,j] := map6to4[i,k6] %[k6]% g6p0[l,k6,m6] %[m6]% map6to4[m6,j]
  # ... but vecless can't presently handle multi-stages

  A[l,i,k] := g6p0[l,i,j] %[j]% map6to4[j,k]
  g4p0[l,i,j] := map6to4[ k,i] %[k]% A[l,k,j]

  A[l,i,k] := g6p1[l,i,j] %[j]% map6to4[j,k]
  g4p1[l,i,j] := map6to4[ k,i] %[k]% A[l,k,j]

  s4 <- predict_hsp_util( g4p0, g4p1, want_LOD_table, k=k)

  if( want_LOD_table) {
    li$LOD6 <- s6@LOD # matrix
    li$PUP6 <- s6@PUP
    li$LOD4 <- s4@LOD
    li$PUP4 <- s4@PUP
    s6@LOD <- s6@LOD <- s4@LOD <- s4@PUP <- NULL
  }

  li <- cbind( li, s6)
  li[ !li$use6, names( s4)] <- s4[ !li$use6,]

  lociar@locinfo <- li
return( lociar)
}
