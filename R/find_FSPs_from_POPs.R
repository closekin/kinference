#' Find full sibling pairs from parent-offspring pairs
#'
#' For pairs already picked as likely parent-offspring pairs (POPs), i.e., those with a weighted pseudo-exclusion (WPSEX) statistic less than some threshold, they might be full sibling pairs (FSPs). This function checks potential POPs with very low WPSEX values for their potential to be FSPs.
#'
#' @param snpg a \code{snpgeno} object
#' @param candiPOPs candidate parent-offspring pairs M*2 matrix of rows in \code{snpg} that have very low WPSEX statistic
#' @param SDwt_POP ???
#' @export
"find_FSPs_from_POPs" <-
function( snpg, candiPOPs, SDwt_POP=0.5) { # shouldn't have default--- hardwire or leave!

  # Don't need full pairwise screening for FSPs (do post hoc on a few hundred
  # candidate POPs), hence all in R.

  define_genotypes()

  # Transform to 4way genotypes
  # based on code in find_duplicates
  # careful, since "factor level" of AB and OO is different in 4way vs 6way
  snpg <- snpg[ c( candiPOPs),]

  # Yet to write 'recode_geno'...
  just_snpg <- snpg

  snpg@diplos <- genotypes4_ambig
  snpg[ just_snpg==AO] <- AAO
  snpg[ just_snpg==AA] <- AAO
  snpg[ just_snpg==BO] <- BBO
  snpg[ just_snpg==BB] <- BBO
  snpg[ just_snpg==OO] <- OO # need to do OO & AB too, since codes are different in 4way vs 6way
  snpg[ just_snpg==AB] <- AB

  n_pairs <- nrow( candiPOPs)
  g1 <- snpg[ 1 %upto% n_pairs,]
  g2 <- snpg[ n_pairs + (1 %upto% n_pairs),]

  # Yet to write 'make_prgeno'...
  # extract.named( make_prgeno( snpg, genotypes4_ambig)) # pA pB pO pgeno[,'AB'] pgeno[,'AAO'] etc
  pA <- snpg@locinfo$pbonzer[,'A']
  pB <- snpg@locinfo$pbonzer[,'B']
  pO <- 1-pA-pB
  n_loci <- length( pA)

  pgeno <- matrix( 0, n_loci, 4, dimnames=list( NULL, genotypes4_ambig))
  pgeno[,AB] <- 2*pA*pB
  pgeno[,OO] <- sqr( pO)
  pgeno[,AAO] <- sqr( pA) + 2*pA*pO
  pgeno[,BBO] <- sqr( pB) + 2*pB*pO

  off <- 1 # until vecless has arbitrary-base arrays
  Pr_same_given_k <- array( 0, c( n_loci, 3))
  Pr_same_given_k[,{off+1}] <- pA * (1-2*pB*(pA+pO)) +
      pB * (1-2*pA*(pB+pO)) +
      pO * ( pA*pA + pB*pB + pO*pO)
  Pr_same_given_k[,{off+2}] <- 1
  Pr_same_given_k[ l, {off+0}] := pgeno[l,g] %[g]% pgeno[l,g]

  Pr_nsame_FSP <- c( 1/4, 1/2, 1/4)
  Pr_same_FSP[ l]:= Pr_nsame_FSP[ k] %[k]% Pr_same_given_k[ l, k]
  Pr_same_POP[ l]:= Pr_same_given_k[ l, {off+1}]

  SD_FSP <- sqrt( Pr_same_FSP * (1-Pr_same_FSP))
  SD_POP <- sqrt( Pr_same_POP * (1-Pr_same_POP))
  SD_denom <- SDwt_POP * SD_POP + (1-SDwt_POP) * SD_FSP
  wt <- (Pr_same_FSP - Pr_same_POP) / sqr( SD_denom)

#  ig1 <- g1
#  storage.mode( ig1) <- 'integer' # otherwise attributes get lost
#  ig1[] <- match( g1@diplos, rownames( mg))
#
#  ig2 <- g2
#  storage.mode( ig2) <- 'integer'
#  ig2[] <- match( g2@diplos, rownames( mg))

  wtsame[i]:= wt[l] %[l]% (g1[i,l]==g2[i,l])
return( wtsame)
}

