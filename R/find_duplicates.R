#' @rdname find_POPs
#' @export
#' @importFrom gbasics sqr
#' @importFrom atease @
#' @importFrom mvbutils cq %upto% %that.are.in% my.all.equal %without.name%
"find_duplicates" <-
function(snpg, subset1=1 %upto% nrow( snpg),
         subset2=subset1, max_diff_genos){

  # Sanity...
stopifnot( is.numeric( subset1) && is.numeric( subset2))
stopifnot( all( !duplicated( subset1)) && all( !duplicated( subset2)))
stopifnot( my.all.equal( subset1, subset2) || !length( intersect( subset1, subset2)))

  # Count #loci with different 4way genos. Errors in 4ways should be low.
  define_genotypes()
  temp_snpg <- snpg
  temp_snpg@diplos <- genotypes4_ambig
  temp_snpg[ snpg==AO] <- AAO
  temp_snpg[ snpg==AA] <- AAO
  temp_snpg[ snpg==BO] <- BBO
  temp_snpg[ snpg==BB] <- BBO
  temp_snpg[ snpg==OO] <- OO # need to do OO & AB too, since codes are different in 4way vs 6way
  temp_snpg[ snpg==AB] <- AB


  # Sort loci to get most informative/random ones first
  # use snpg1 for this, arbitrarily
  gtab <- matrix( 0, 4, ncol( snpg), dimnames=list( genotypes4_ambig, NULL))
  for( ig in genotypes4_ambig) {
    gtab[ ig,] <- colSums( temp_snpg==ig)
  }
  gtab <- gtab / nrow( snpg)
  pid <- colSums( sqr( gtab))
  o <- order( pid)
  pid <- pid[ o]

  # Remove extranea
  temp_snpg <- temp_snpg[ ,o]
  attributes( temp_snpg) <- attributes( temp_snpg)[ 'dim']
  temp_snpg <- t( temp_snpg)

  # Trying special-cases here to minimize copying
  if( my.all.equal( subset1, subset2)) {
    if( !my.all.equal( subset1, 1 %upto% ncol( temp_snpg))) {
      temp_snpg <- temp_snpg[, subset1]
    }

    result <- DUP_paircomps_lots(
        geno1= temp_snpg,
        geno2= temp_snpg,
        symmo= TRUE,
        max_diff_genos = max_diff_genos
      )
  } else { # different subsets
    result <- DUP_paircomps_lots(
        geno1= temp_snpg[ , subset1],
        geno2= temp_snpg[ , subset2],
        symmo= FALSE,
        max_diff_genos = max_diff_genos
      )
  }


  result$bigs <- with( result, data.frame( ndiff=big_similar, i=big_i, j=big_j))
  result <- result %without.name% cq( big_similar, big_i, big_j)
  result$call <- sys.call()

return( result)
}
