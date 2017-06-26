#' @importFrom mvbutils cq %except% %not.in%
#' @importFrom atease @ @<-
#' @importFrom gbasics make_genopairer sqr
#' @export
"prepare_PLOD_SPA" <- function( geno6, n_pts_SPA_renorm=201) {
# To be run after hsp_power( ..., want_LOD_table=TRUE)
# n_pts_SPA_renorm should really be as big as R can handle without running out memory
# ... but 201 should be OK I guess. If 201 and 301 give almost-identical results then all well!

stopifnot( all( cq( LOD4, LOD6, use6) %in% names( geno6@locinfo)))

  og <- options( vecless.print=FALSE)
  on.exit( options( og))

  # Combine 4way and 6way stuff into overall LOD and PUP
  extract.named( geno6@locinfo[ cq( use6, LOD6, LOD4, PUP6, PUP4)])
  use4 <- !use6

  # DO NOT change actual genos though; they will be changed on-the-fly prior to kin-finding
  # ... code WOULD be this:
  # temp_snpg <- snpg
  # recode4to6temp <- function( x) { x[ x=='AO'] <- AA; x[ x=='BO'] <- BB; x}
  # temp_snpg[ , use4] <- recode4to6temp( snpg[, use4]) # (AA,AO) -> AA; (BB,BO) -> BB

  LOD <- LOD6
  PUP <- PUP6
  cn6 <- colnames( LOD6)
  cn4 <- colnames( LOD4)

  # Change only the entries with "full homz" since "single nulls" won't be accessed
  for( ichangio in grep( 'AA|BB', cn6) %except% grep( 'AO|BO', cn6)){
    was1 <- substring( cn6[ ichangio], 1, 2)
    was2 <- substring( cn6[ ichangio], 4, 5)
    now1 <- sub( '(A|B)\\1', '\\1\\1O', was1)
    now2 <- sub( '(A|B)\\1', '\\1\\1O', was2)
    iget4 <- paste( now1, now2, sep='/')
    if( iget4 %not.in% cn4) { # reverse the order
      iget4 <- paste( now2, now1, sep='/')
    }

    LOD[ use4, ichangio] <- LOD4[ use4, iget4]
    PUP[ use4, ichangio] <- PUP4[ use4, iget4]
  }
  # For safety's sake, LOD( XO,...) := NA; should never be accessed
  hasO <- grep( '(A|B)O', cn6)
  LOD[ use4, hasO] <- NA # security in case of wrong access later for real data
  PUP[ use4, grep( '(A|B)O', cn6)] <- 0
  LOD@mg <- make_genopairer( geno6@diplos)

  make_K <- function( PUP, LOD) { # ... while the sun skines

      # this to be sorted with replacing 'vecless' with 'kinference'
      # BUUTT is isn't :(
      # vecless **should** work just exorting := BUT doesn't seem to
      e <- new.env( parent=asNamespace( 'vecless'))
      # add sqr to the environment so that vecless can see it...
      e$sqr <- gbasics::sqr
      e$renorm_SPA_cumul <- renorm_SPA_cumul
      e$PUP <- PUP
      e$LOD <- e$LODOK <- LOD
      e$LODOK[ is.na( LOD)] <- 0 # leaving NAs in would mess up the calcs
      e$n_pts_SPA_renorm <- n_pts_SPA_renorm

      evalq( envir=e, {
        if( !nrow( PUP)) {
          K <- dK <- ddK <- function( tt) 0*tt
          inv_CDF <- CDF <- function( x) NA+x
        } else {
          PUPLOD <- PUP * LODOK
          PUPLOD2 <- PUPLOD * LODOK

          K <- function( tt) {
            ETT[ it, l, g12] := exp( tt[ it] * LODOK[ l, g12])
            S[ it, l] := PUP[ l, g12] %[g12]% ETT[ it, l, g12]
            rowSums( log( S))
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

          extract.named( renorm_SPA_cumul( K, dK, ddK, n_pts=n_pts_SPA_renorm))
          # ... CDF and inv_CDF but *not* for extreeeme values
        } # if nrow PUP
      })

      # K and co will know PUP & co thru enviro magic
    return( e)
    } # function make_K


  Kenv <- make_K( PUP, LOD)
  geno6@Kenv <- Kenv

  # PUP and LOD are in Kenv now, so don't duplicate them in locinfo
  # They are really the "workhorse" versions, and are a bit cheaty, so don't want
  # them too public

return( geno6)
}
