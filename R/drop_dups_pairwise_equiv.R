#' Equivalence class sifter
#'
#' Constructs equivalence classes from pairwise equivalences, and returns the "surplus" elements; if you then drop those elements, only one element from each eq-class will be retained. Requires 2-col matrix showing equivalent pairs.
#'
#' Code is taken from Numerical Recipes so I should rewrite it perhaps (original algorithm is by Knuth).
#'
#' @param ij 2-column matrix or data.frame; probably "row numbers" in a dataset, though might work with character strings too
#' @param want_groups if \code{TRUE}, also return the equivalence-classes themselves, as attribute \code{groups}.
#'
#' @return Surplus elements in \code{ij}, perhaps plus attributes \code{groups} if \code{want_groups=TRUE}. You can look at that to figure out which elements are being retained (one "representative" from each equiv class).
#'
#' @examples
#' pairs <- matrix( c(
#' 294, 289,
#' 328, 294,
#' 904, 857,
#' 905, 904),
#'     ncol=2, byrow=TRUE)
#' drop_dups_pairwise_equiv( pairs, TRUE)
#' #[1] 289 328 857 905
#' #attr(,"groups")
#' #attr(,"groups")$`5`
#' #[1] 294 328 289
#' #
#' #attr(,"groups")$`6`
#' #[1] 904 905 857
#' @importFrom mvbutils do.on %except% FOR
#' @export
"drop_dups_pairwise_equiv" <- function( ij, want_groups=FALSE) {
  ij <- as.matrix( ij) # in case it was a data.frame
  uij <- unique( c( ij))
  ij[] <- match( ij, uij)

  n <- max( ij)
  m <- nrow( ij)
  nf <- 1:n # Initialize each element its own class

  for( l in 1:m) {
    j <- ij[l,1]
    while( nf[ j] != j) {
      j <- nf[ j]
    }
    k <- ij[l,2]
    while( nf[ k] != k) {
      k <- nf[ k]
    }
    if( j != k) {
      nf[ j] <- k
    }
  }

  for( j in 1:n) {
    while( nf[ j] != nf[ nf[ j]]) {
      nf[ j] <- nf[ nf[ j]]
    }
  }

  # Keep first member of each nf-group
  groups <- split( 1:n, nf)
  keeps <- do.on( groups, .[1])
  drops <- (1:n) %except% keeps
  drops <- sort( uij[ drops])

  if( want_groups) {
    drops@groups <- FOR( groups, uij[.])
  }

return( drops)

#for (k=1;k<=n;k++) nf[k]=k; Initialize each element its own class.
#for (l=1;l<=m;l++) { For each piece of input information...
#  j=lista[l];
#  while (nf[j] != j) j=nf[j]; Track 1st element up to its ancestor.
#  k=listb[l];
#  while (nf[k] != k) k=nf[k]; Track second element up to its ancestor.
#  if (j != k) nf[j]=k; // If they are not already related, make them so
#}
#for (j=1;j<=n;j++) Final sweep up to highest ancestors.
#while (nf[j] != nf[nf[j]]) nf[j]=nf[nf[j]];
#}
#Alternatively, we may be able to construct a function
#

}
