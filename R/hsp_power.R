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

  # map the 6 way genos to 4 way
  remap <- map6to4(g6p0, g6p1)

  s4 <- predict_hsp_util( remap$g4p0, remap$g4p1, want_LOD_table, k=k)

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
