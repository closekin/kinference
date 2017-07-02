#' Check individual genotypes for aggregate typicality
#'
#' Compute lglk of entire 4way genotype of each individual, ie sum log Pr[ g(i,l)]; and compare the distro of lglks across individuals with its predicted shape given ALFs. Significant mismatch is bad. Can also detect outliers. Up to you what criterion to use for that. NB running \code{locator(2)} lets you click the histogram to figure out where you'd like to cut.
#'
#' Currently, the SPA calcs are a wee bit slow because of heavy use of \code{vecless} which in version 1.0 is sluggish. The lglks themselves are computed in C and are blisteringly fast.
#'
#' @param snpg a \code{snpgeno} (6way genotype)
#' @param indiv_lglk_hist_pars list like in \code{dump_badhetz_fish}, for controlling histogram
#'
#' @return Vector of lglk for each individual. I haven't added any formal uh-oh criteria yet; that could be done via the SPA, as in \code{dump_badhetz_fish}. But, reading off from the graph is probably fine...
#'
#' @importFrom atease @
#' @importFrom gbasics snpgeno
#' @importFrom mvbutils cq extract.named
#' @importFrom gbasics sqr
#' @import vecless
#' @export
"ilglk_geno" <- function(snpg, indiv_lglk_hist_pars=list(), quick=TRUE) {
  define_genotypes()
  extract.named( snpg@locinfo[ cq( pbonzer)])

  p0 <- pbonzer[,'O'] + pbonzer[,'C']
  pA <- pbonzer[,'A']
  pB <- pbonzer[,'B']

  n_samps <- nrow( snpg)
  n_loci <- ncol( snpg)
  snpg4 <- snpgeno( n_samps, n_loci, genotypes4_ambig,
      info=snpg@info[,cq( Our_plate, Our_sample)],
      locinfo=snpg@locinfo[,cq( Locus), drop=FALSE])
  snpg4[ snpg==OO] <- OO
  snpg4[ snpg==AB] <- AB
  snpg4[ snpg==AA] <- AAO
  snpg4[ snpg==AO] <- AAO
  snpg4[ snpg==BB] <- BBO
  snpg4[ snpg==BO] <- BBO

  pgeno <- matrix( 0, n_loci, 4, dimnames=list( NULL, genotypes4_ambig))
  pgeno[ , OO] <- sqr( p0)
  pgeno[ , AB] <- 2*pA*pB
  pgeno[ , AAO] <- 2*pA*p0 + pA*pA
  pgeno[ , BBO] <- 2*pB*p0 + pB*pB
  lpgeno <- log( pgeno) # lambda in doco


  # These are slowish in vecless 1.0
  # ... but less error-prone to write
  K <- function( tt) {
      ttp1 <- tt+1
      KK[ j] := SUM_ %[l]% log( SUM_ %[g]% exp( ttp1[ j] * lpgeno[ l, g]))
    return( c( KK))
    }

  dK <- function( tt) {
      ttp1 <- tt+1
      etp1l[ j, l, g] := exp( ttp1[ j] * lpgeno[ l, g])
      num[ j, l] := lpgeno[ l, g] %[g]% etp1l[ j, l, g]
      denom[ j, l] := SUM_ %[g]% etp1l[ j, l, g]
      dKK[ j] := SUM_ %[l]% (num[ j, l] / denom[ j, l])
    return( c( dKK))
    }


  ddK <- function( tt) {
      ttp1 <- tt+1
      etp1l[ j, l, g] := exp( ttp1[ j] * lpgeno[ l, g])
      num[ j, l] := lpgeno[ l, g] %[g]% etp1l[ j, l, g]
      denom[ j, l] := SUM_ %[g]% etp1l[ j, l, g]
      num2[ j, l] := sqr( lpgeno[ l, g]) %[g]% etp1l[ j, l, g]
      dKK[ j] := SUM_ %[l]% ( num2[ j, l] / denom[ j, l] -  sqr( num[ j, l] / denom[ j, l]))
    return( c( dKK))
    }

  if( quick) {
    K <- compile_vecless( K( -1))
    dK <- compile_vecless( dK( -1))
    ddK <- compile_vecless( ddK( -1))
  }

  if( FALSE) { # Checks: do manually in mtrace
    ntest <- 1000
    Ktest <- function( tt) { # scalar
      Ksim <- meansim <- rep( 0, ntest)
      for( l in 1:n_loci) {
        # Can't directly sample from genotypes4_ambig since can't matrix-subscript mixed int and char
        genos <- rsample( ntest, seq_along( genotypes4_ambig), prob=pgeno[l,], replace=TRUE)
        lp <- lpgeno[ cbind( l, genos)]
        Ksim <- Ksim + ( tt * lp) # actually log( exp( t*lp))
        meansim <- meansim + lp # though see below for better way to check!
      }
      returnList( Ksim, meansim)
    }

    dK( 0)
    sum( pgeno * lpgeno) # should be the same
  }

  ilglk <- indiv_lglk_geno(
      lpgeno= lpgeno,
      geno= snpg4)

  if( FALSE) { # "manual" check on calcs
    lp <- lpgeno[ cbind( rep( 1 %upto% n_loci, n_samps), snpg)]
    dim( lp) <- c( n_loci, n_samps)
    ilglk_manual <- colSums( lp)
  }

  # inv_CDF <- renorm_SPA_cumul( K, dK, ddK)$inv_CDF
  dens_SPA <- renorm_SPA( K, dK, ddK, 'func')

  indiv_lglk_hist_pars <- add_list_defaults( indiv_lglk_hist_pars,
      main='Geno lglk by FISH', #sprintf( 'Geno lglk by FISH: multhresh=%5.2f', method, multhresh_indiv_lglk_fish),
      xlim= range( ilglk),
      xlab='', nclass=50)
  lv <- do.call( 'hist', c( list( x=ilglk), indiv_lglk_hist_pars))
  with( lv, # then plot predicted density. Slowish with vecless 1.0
    lines( mids, diff( breaks) * dens_SPA( mids) * n_samps, col='green')
  )

return( ilglk)
}
