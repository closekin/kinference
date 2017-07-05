#' @importFrom atease @ @<-
#' @importFrom mvbutils do.on
#' @export
"find_FSPs_from_HSPs" <- function( snpg, HSPs) {
  # For pairs already picked as HSPs, ie PLOD(HSP,UP) > eta: they might be FSPs

  # Don't need full pairwise screening for FSPs (do post hoc on a few hundred
  # HSPs), hence all in R.

  # HSPs should be M*2 matrix of rows in snpg that are HSPs or FSPs


  define_genotypes()

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


  # what is happening here? Is this magick?! Need to parcel up somewhere else.
  OLOD <- LOD4
  OPUP <- PUP4
  mg <- OLOD@mg

  # Can't quite do this with vecless!
  OLOD[ is.na( OLOD)] <- 0 # set to NA for 4way loci
  n_loci <- nrow( OPUP)
  PUP4 <- LOD4 <- array( 0, c( n_loci, 4, 4))
  for( ig in 1:4) {
    gjseq <- mg[ , ig]
    XXi[ l, gj] := OPUP[ l, gj=gjseq] # Shouldn't work with new vecless syntax... but does !?
    PUP4[ l, {ig}, gj] := XXi[ l, gj]
    # NPUP[,ig,] <- OPUP[ , mg[,ig]]

    XXi[ l, gj] := OLOD[ l, gj=gjseq]
    LOD4[ l, {ig}, gj] := XXi[ l, gj]
  }

  # recover the column/row/etc names
  dimnames(PUP4) <- dimnames(LOD4) <- list(NULL, rownames(mg), rownames(mg))

  ## calculate the P[g1, g2 | k shared alleles]
  # since we save the LOD (for HSP vs UP) and the P[UP] calculate
  # P[g1 g2 | HSP] from that
  PHSP4 <- exp( LOD4) * PUP4 # Pr[gg|HSP] <- 0.5 * PUP4 + 0.5 * Pr[gg|kappa=1]
  P_k0 <- PUP4
  P_k1 <- 2*PHSP4 - PUP4
  P_k1[ P_k1 < 0] <- 0 # rounding error
  # P_k2 <- sqrt(diag(P_k0)) <- we are doing this-ish
  P_k2 <- 0 * P_k0 # get the shape right
  P_k2[l, i, i] := sqrt(P_k0[l, i, i])


  nsib <- nrow( HSPs)
  nloci <- ncol(sibg)

  # if only R were zero-indexed...
  kappa_fsp <- c(1/4, 1/2, 1/4)
  kappa_hsp <- c(1/2, 1/2, 0)

  p12fsp <- p12hsp <- matrix(NA, nloci, nsib)


  # split sibg into g1 and g2 for the two parts of the pairs
  g1 <- sibg[1:nsib, ]
  g2 <- sibg[(nsib+1):(2*nsib), ]

  # for loop version of the code
  g1 <- as.character(g1)
  g2 <- as.character(g2)
  for(i in 1:nsib){
    for(l in 1:nloci){
      # P[g_1l g_2l | FSP]
      p12fsp[l, i] <- kappa_fsp[1] * P_k0[l, g1[i, l], g2[i, l]] +
                      kappa_fsp[2] * P_k1[l, g1[i, l], g2[i, l]] +
                      kappa_fsp[3] * P_k2[l, g1[i, l], g2[i, l]]
      # P[g_1l g_2l | HSP]
      p12hsp[l, i] <- kappa_hsp[1] * P_k0[l, g1[i, l], g2[i, l]] +
                      kappa_hsp[2] * P_k1[l, g1[i, l], g2[i, l]] +
                      kappa_hsp[3] * P_k2[l, g1[i, l], g2[i, l]]
    }
  }
  OD_FH <- p12fsp/p12hsp
  LOD_FH <- log(OD_FH)
  PLOD_FH <- colSums(LOD_FH)


  # some vecless hieroglyphs
  #P_k <- array(0, c(3, dim(P_k0)))
  #P_k[ 1, l, gi, gj] := P_k0[ l, gi, gj]
  #P_k[ 2, l, gi, gj] := P_k1[ l, gi, gj]
  #P_k[ 3, l, gi, gj] := P_k2[ l, gi, gj]

  #Pr_FSP[ l, gi, gj] := kappa_fsp[ k] %[k]% P_k[ k, l, gi, gj]
  #Pr_HSP[ l, gi, gj] := kappa_hsp[ k] %[k]% P_k[ k, l, gi, gj]

  #LOD_FH <- log( Pr_FSP / Pr_HSP)

  ## Would like but can't have yet:
  ## PLOD_FH[ ipair] := SUM_ %[l]% LOD_FH[ l, g1[ ipair, l], g2[ ipair, l] ]

  ## Instead do the lookup as an ugly dot product...
  ## jseq kseq these must cover the 2nd and 3rd index-range of thingo
  #jseq <- 1:4
  #kseq <- 1:4
  #ig1 <- g1
  #storage.mode( ig1) <- 'integer' # otherwise attributes get lost
  ##ig1[] <- match( g1@diplos, rownames( mg))

  #ig2 <- g2
  #storage.mode( ig2) <- 'integer'
  ##ig2[] <- match( g2@diplos, rownames( mg))

  #PLOD_FH[i]:= SUM_ %[l]% (LOD_FH[ l,j,k] %[j,k]% (jseq[j]==ig1[i,l] &
  #                                                 kseq[k]==ig2[i,l]))

  ## unnecessary
  # EPLOD_FH_F_by_locus[ l] := Pr_FSP[ l, gi, gj] %[gi,gj]% LOD_FH[ l, gi, gj]
  # EPLOD_FH_F[] := SUM_ %[l]% (Pr_FSP[ l, gi, gj] %[gi,gj]% LOD_FH[ l, gi, gj])

  #EPLOD_FH_F <- sum( Pr_FSP * LOD_FH)
  #EPLOD_FH_H <- sum( Pr_HSP * LOD_FH)


  # format a return object
  ret <- list()
  ret$bigs <- data.frame(PLOD_FH = PLOD_FH,
                         i       = HSPs[,1],
                         j       = HSPs[,2])

  return(ret)
}
