# This is package kinference 

#".onLoad" <-
#function( libname, pkgname) {
#  oa <- base::system.file( sprintf( '/R/load_%s_dll.R', pkgname),
#      package=pkgname, lib.loc=libname)
#  if( !nzchar( oa)) {
#    oa <- base::system.file( sprintf( '/R/load_%s_dll.r', pkgname),
#        package=pkgname, lib.loc=libname)
#  }
#
#  base::source( oa, local=TRUE)
##}

#' kinference
#'
#' Kin-finding pairwisely for close-kin mark-recapture
#'
#' @docType package
#' @author Mark V Bravington, David L Miller, Shane M Baylis
#' @import Rcpp atease mvbutils gbasics vecless debug abind BH
#' @importFrom Rcpp sourceCpp evalCpp
#' @useDynLib kinference
#' @name kinference
NULL

"old.onLoad" <-
function( libname, pkgname) {
  cat( 'Am I loaded?', pkgname %in% loadedNamespaces(), '\n')
  if( !exists( 'onload_autowrap', asNamespace( pkgname), inherits=FALSE, mode='function')) {
    # then source() the R script and create it here...
    my_dll <- getOption( sprintf( '%s_debug_C', FALSE))
    oa <- system.file( 'R/onload_autowrap.R', package=pkgname, lib.loc=libname)
    eval( substitute( source( oa, local=TRUE), list( oa=oa)), asNamespace( pkgname))
  }

  onload_autowrap( pkgname, libname)
}  ## ported from MVB_kinference 25/9/18


#' add_pairprob_error(): Bare documentation
#'
#' This function has only the bare minimum of documentation necessary for roxygen to
#' parse it. We should probably add some proper documentation here.
#' @param nlocal a param
#' @export

"add_pairprob_error" <-
function( nlocal=sys.parent()) mlocal({
## Needs pp_true and snerr

  g_err <- named( genotypes_C) ? 0 # recorded because it makes perr easy to define
  perr <- 0*nchar( g_err) ? 1 # named

  g_1 <- substring( genotypes_C, 1, 1)
  g_2 <- substring( genotypes_C, 2, 2)
  AA_BB <- g_1==g_2 & (g_1==A | g_1==B)
  g_err[ AA_BB] <- g_1[ AA_BB] %&% 'O'
  perr[ AA_BB] <- snerr[ sprintf( '%1$s%1$s2%1$sO', g_1[ AA_BB])] ? 0 # XX2XO

  AO_BO <- (g_1==A | g_1==B) & g_2==O
  g_err[ AO_BO] <- g_1[ AO_BO] %&% g_1[ AO_BO]
  perr[ AO_BO] <- snerr[ sprintf( '%1$sO2%1$s%1$s', g_1[ AO_BO])] ? 0 # XO2XX

  # Original version: not sure it's selfrecordable
#  make_pp_err <- function( g1, g2) {
#      g1_err <- g_err[ g1]
#      g2_err <- g_err[ g2]
#      pp_true[ cbind( g1, g2) ] * (1-perr[ g1]) * (1-perr[ g2]) +
#      pp_true[ cbind( g1_err, g2) ] * perr[ g1_err] * (1-perr[ g2]) +
#      pp_true[ cbind( g1, g2_err) ] * (1-perr[ g1]) * perr[ g2_err] +
#      pp_true[ cbind( g1_err, g2_err) ] * perr[ g1_err] * perr[ g2_err]
#    }
#
  # pp_err[] <- outer( named( genotypes_C), named( genotypes_C), make_pp_err)

  # Need 2-stage setup here so selfrec works: "pure" assignment can't be overloaded in R
  pp_err <- matrix( -1, length( genotypes_C), length( genotypes_C), dimnames=rep( list( genotypes_C), 2)) ? 1

  gg <- list() ? 0 # this is a trick to allow with(gg,...) in playback, since the "real" gg isn't kept

  gg <- expand.grid( g1=genotypes_C, g2=genotypes_C, stringsAsFactors=FALSE)
  gg <- within( gg, {
     g1_err <- g_err[ g1]
     g2_err <- g_err[ g2]
  })
  extract.named( gg)

  pp_err[ cbind( g1, g2)] <-
      pp_true[ cbind( g1, g2) ] * (1-perr[ g1]) * (1-perr[ g2]) +
      pp_true[ cbind( g1_err, g2) ] * perr[ g1_err] * (1-perr[ g2]) +
      pp_true[ cbind( g1, g2_err) ] * (1-perr[ g1]) * perr[ g2_err] +
      pp_true[ cbind( g1_err, g2_err) ] * perr[ g1_err] * perr[ g2_err]        ? 0


  # Merge C with O. Has to be mlocal() to avoid <<- which buggers playback. Things before 'nlocal' are args; things after are temporaries
  add_up <- function( i, result, ..., nlocal=sys.parent(), dimor, other)  mlocal({
      result <- match( result, genotypes_C)
      dimor <- slice.index( pp_err, i)

      for( other in FOR( list( ...), match( ., genotypes_C))) {
        pp_err[ dimor==result] <- pp_err[ dimor==result] + pp_err[ dimor==other] ? 0
        pp_err[ dimor==other] <- 0 # avoid double-use ? 0
      }
    })

  for( i in 1:2) {
    add_up( i, OO, CC, CO)
    add_up( i, AO, AC)
    add_up( i, BO, BC)
  }

# IE:
# pp_err2[,OO] <- pp_err2[,CC] + pp_err2[,CO] + pp_err2[,OO]
# pp_err2[,CC] <- pp_err2[,CO] <- 0
# pp_err2[AO,] <- pp_err2[AC,] + pp_err2[AO,]
# pp_err2[AC,] <- 0

  pp6_err <- pp_err[ genotypes6, genotypes6] ? 0
})

#' calc_g6probs_IBD0_scalar(): Bare documentation
#'
#' This function has only the bare minimum of documentation necessary for roxygen to
#' parse it. We should probably add some proper documentation here.
#'
#' @param P a param
#' @param snerr a param
#' @param record TRUE or FALSE
#' @export
#' @importFrom vecless set_recording
#' @importFrom mvbutils %&%

"calc_g6probs_IBD0_scalar" <-
function( P, snerr, record=FALSE) {
## SCALAR-ONLY VERSION... this is hard enough!
## Though can be called with 1-row matrix args, eg with( x@locinfo[1,], calc_g6probs_IBD1( pbonzer, snerr))

  # snerr = P( misclassifying true XX as XO, and vice versa)

  set_recording( cq( P, snerr, pp_true, pr2, pp_err, perr, pp6_err), record)
  define_genotypes()           ? 0

  P <- drop( P) # for scalar version
  P <- P              ? 0
stopifnot( my.all.equal( names( P), names( ABCO)))

  snerr <- drop( snerr) # for scalar version

  snerr <- snerr            ? 0
  pp_true <- matrix( 0, length( genotypes_C), length( genotypes_C),
                    dimnames=list( genotypes_C, genotypes_C))         ? 1


  g_1 <- substring( genotypes_C, 1, 1)
  g_2 <- substring( genotypes_C, 2, 2)
  pr2 <-  nchar( named( genotypes_C))    ? 1 # named

  is_het <- g_1 != g_2
  pr2[] <- P[ g_1] * P[ g_2]             ? 0
  pr2[ is_het] <- 2 * pr2[ is_het]       ? 0

  pp_true <- matrix( 0, length( genotypes_C), length( genotypes_C),
                    dimnames=rep( list( genotypes_C), 2))        ? 1

  extract.named( expand.grid( gp1=genotypes_C, gp2=genotypes_C,
                             stringsAsFactors=FALSE))
  pp_true[ cbind( gp1, gp2)] <- pr2[ gp1] * pr2[ gp2]                       ? 0

  # NB that for this UP case, XX/XO errors shouldn't change the overall probs because the cutoffs are chosen to do exactly that!
  add_pairprob_error()

return( pp6_err)
}


#' calc_g6probs_IBD1_scalar(): Bare documentation
#'
#' This function has only the bare minimum of documentation necessary for roxygen to
#' parse it. We should probably add some proper documentation here.
#'
#' @param P a param
#' @param snerr a param
#' @param record TRUE or FALSE
#' @export
#' @importFrom vecless set_recording

"calc_g6probs_IBD1_scalar" <-
function( P, snerr, record=FALSE) {
## SCALAR-ONLY VERSION... this is hard enough!
## Though can be called with 1-row matrix args, eg with( x@locinfo[1,], calc_g6probs_IBD1( pbonzer, snerr))

  set_recording( cq( P, snerr, pp_true, pp3, pp_err, perr, pp6_err), record)
  define_genotypes()           ? 0

  P <- drop( P) # for scalar version
  P <- P               ? 0
stopifnot( my.all.equal( names( P), names( ABCO)))

  snerr <- drop( snerr) # for scalar version
  snerr <- snerr        ? 0
  pp_true <- matrix( 0, length( genotypes_C), length( genotypes_C),
                    dimnames=list( genotypes_C, genotypes_C))          ? 1

  # Note: AB/AB can share either A or B... so have to accumulate

  # always write as S,U1 and S,U2 so probs are unambig
  # then accumulate into sorted versions since *several* S,U1,U2 could contribute to one XY/ZW
  # ... or *none*, of course, if there is no shared allele in XY/ZW

  # 3-col DF of (shared, 1st unshared, 2nd unshared) alleles in the pair
  su1u2 <- as.matrix( expand.grid( ABCO, ABCO, ABCO))
  implied <- su1u2[ ,c(1,2,1,3)]
  swap12 <- implied[,2] < implied[,1]
  implied[ swap12, 1:2] <- implied[ swap12, 2:1]
  swap34 <- implied[,4] < implied[,3]
  implied[ swap34, 3:4] <- implied[ swap34, 4:3]
  cimplied <- cbind( implied[,1] %&% implied[,2], implied[,3] %&% implied[,4])

  pp3 <- P[ su1u2[,1]] * P[ su1u2[,2]] * P[ su1u2[,3]]             ? 0

  # Each row can contribute to at most 2 implied XY/ZW combinations
  # eg AB/AB <- ABB or BAA
  # AB/CD <- nothing
  # AB/AC <- ABC

  implstr <- cimplied[,1] %&% cimplied[,2]
  has_impl2 <- which( duplicated( implstr))

  first_impl <- match( implstr[ -has_impl2], implstr)
  second_impl <- match( implstr[has_impl2], implstr, 0)
  cimpl2 <- cimplied[ second_impl,]

  pp_true[ cimplied[ first_impl,]] <- pp3[ -has_impl2]            ? 0
  pp_true[ cimpl2] <- pp_true[ cimpl2] + pp3[ has_impl2]          ? 0

  # Allow for XX <-> XO errors--- hopefully the only ones! (Watch out for scaffoldy version)
  # Assumes errors INDEPENDENT even tho there could be heritability (eg weak grabbing of mutated primer)
  # CC <-> CO is ignored in this version, since C gets scored as O

  add_pairprob_error()

return( pp6_err)
}

#' Find chains in HSPs; summarize sib-groups
#'
#' Find chains of relatives of fish `seed`.
#'
#' For checking veracity of \emph{potential} half-sibs or other kin-pairs.
#' \code{chain_pairwise} organizes them into chains within which each sample
#' can be linked to another by a succession of direct pairwise links. The
#' general idea is that real HSPs will be in clusters of 2 or 3; a spurious
#' sample with a "lucky" genotype that wants to be everybody's mate will appear
#' in a big but incomplete chain of mostly false-positives, where the direct
#' pairwise links between the other chain-members are weak. You'd only run
#' \code{chain_pairwise} for pairs with a PLOD (or whatever statistic is being
#' used) within a particular suspect range, so each chain may have
#' false-negatives (i.e. missing direct links), but the general idea should be
#' clear.
#' 
#' \code{get_chain} finds the chain for one specific sample.
#' 
#' @aliases get_chain chain_pairwise
#' @param thing output from \code{find_HSPs} or \code{find_POPs} etc, or some
#'              subset thereof
#' @param seed one sample ID, interpreted as a row-number in \code{thing}. To
#'             do:also allow names, via \code{info} attr.
#' @importFrom mvbutils %is.not.a% %where%
#' @return \code{chain_pairwise} returns a list of matrices, each for one
#' chain; the rows and columns of each matrix are the samples in that chain. A
#' "+" in the matrix indicates that those two samples have a direct pairwise
#' link (i.e., they appear together in one row of \code{thing}); a "." means
#' not. The rows and columns of each matrix are sorted so that the linkiest
#' samples are on the bottom and right. \code{get_chain} returns the row-subset
#' of \code{thing} that is chained to \code{seed}.
#' @keywords misc
#' @export

"chain_pairwise" <-
function( thing) {
  extract.named( thing[cq(i,j)])
  
  ij <- sort( unique( c( i, j)))
  ijpairs <- sprintf( '%i.%i', c(i,j), c(j,i))
  
  is_pair <- function( x, y) {
      test <- sprintf( '%i.%i', x, y) %in% ijpairs
      xtest <- test
      xtest[ test] <- '+'
      xtest[ !test] <- '.'
    return( noquote( xtest))
    }
  
  # Split into chains
  to.do <- ij
  chains <- pairmats <- list()
  while( length( to.do)) {
    this_chain <- get_chain( thing, to.do[1])
    chains <- c( chains, list( this_chain))
    i <- this_chain$i
    j <- this_chain$j

    ij <- sort( unique( c( i, j)))
    ijpairs <- sprintf( '%i.%i', c(i,j), c(j,i))
    
    to.do <- to.do %except% c( i, j)
    this_pairs <- outer( ij, ij, is_pair)
    dimnames( this_pairs) <- rep( list( as.character( ij)), 2)
    o <- order( rowSums( this_pairs=='+'))
    this_pairs <- this_pairs[ o, o]
    
    pairmats <- c( pairmats, list( this_pairs))
  }
  
  # Biggest chains last--- easiest to see!
pairmats[ order( do.on( pairmats, nrow( .)))]
}  ## chain_pairwise ported from MVB_kinference 24/11/18

#' check_FPosity(): QC for kin-finding
#'
#' Return predicted mean & variance of CLODs for each sample,
#' ie how prone is that sample's particular genotype to yielding
#' unusually high/low PLODs when compared with a random unrelated
#' sample. Then you can turn this into a prediction of each
#' sample's per-comp chance of yielding a False-Positive with
#' another Unrelated sample. This looks like quite a powerful
#' diagnostic, but is not fully explored yet.
#' The document "d:/docs/genetics/Dart/sbt-baits-v3/too-many-plods.lyx"
#' has more info in section 4.1 on "rat CLODs".
#' 
#' @param snpg a 'snpgeno' object.
#' @param nsim currently inactive. A simulation option exists in the code
#'             to check the null distro (not much use for far tails, of
#'             course).
#' @return Dataframe with columns "ECLOD" and "VCLOD". See examples format.
#' @importFrom mvbutils cq %without.name% %&%
#' @importFrom gbasics sqr
#' @importFrom atease @ @<-
#' @examples
#' ## Rough chance of yielding a PLOD>5, say
#' # cloddo <- check_FPosity( snpg = snpg)
#' # Pr_FPos_5 <- pnorm( 5, mean=cloddo$ECLOD, sd=sqrt( cloddo$CLOD), lower=FALSE)
#' # hist( Pr_Fpos_5, nc=50)
#' ## highlight some known suspects
#' # abline( v=Pr_Fpos_5[ suspects], col='red')
#' @export

"check_FPosity" <-
    function( snpg, nsim=0){
## snpg should have been thru 'prepare_PLOD_SPA' so it has @PPS
stopifnot( 'Kenv' %in% names( attributes( snpg)))

  og <- options( vecless.print=FALSE)
  on.exit( options( og))

  kinference::define_genotypes()
  for( iwhat in cq( LOD, PUP, PUPLOD, PUPLOD2)) {
    assign( 'O' %&% iwhat, snpg@Kenv[[ iwhat]])
  }
  mg <- OLOD@mg

  useN <- snpg@locinfo$useN
## use4 <- !use6
  temp_snpg <- snpg
  recode4to6temp <- function( x) { x[ x=='AO'] <- AA; x[ x=='BO'] <- BB; x}
  recode3to6temp <- function( x) { x[ x=='AO'] <- AA; x[ x=='BO'] <- BB; x[ x=='OO'] <- BB; x}
  temp_snpg[ , useN == 4] <- recode4to6temp( snpg[, useN == 4]) # (AA,AO) -> AA; (BB,BO) -> BB
  temp_snpg[ , useN == 3] <- recode3to6temp( snpg[, useN == 3]) # and similar, for 3-way

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

return( data.frame( ECLOD=colSums( my_e_CLOD), VCLOD=colSums( my_v_CLOD)))
  ## everything after the return() is ignored - does useN actually have any impact on this fun?

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

  result <- with( result, data.frame( PLOD=big_PLOD, i=big_i, j=big_j))
  result <- result %without.name% cq( big_PLOD, big_i, big_j)

  result@bins <- bins
  result@binprobs <- binprobs
  result@eta <- eta
  result@keep_thresh <- keep_thresh
  result@call <- sys.call()

return( result)
}


#' crapometer(): splice hetzminoo_fancy and ilglk_geno for ?better crap-detection
#'
#' crapometer is a bivariate combination of hetzminoo_fancy and ilglk_geno, which may be 
#' better at identifying 'suspect' samples. It may or may not end up in the final
#' kinference toolchain.
#'
#' @param snpg a param
#' @param focusees a param
#' @param boring a param, default \code{1 \%upto\% nrow(snpg)}. See mvbutils for \code{ \%upto\% }.
#' @export

"crapometer" <-
function( snpg, focusees, boring=1 %upto% nrow( snpg)) {
## Bivariate combo of hetzminoo_fancy and ilglk_geno; perhaps this will be better for
# finding suspect samples. But seemingly not much (for school shark), so I haven't exported itfi.

  hmfr <- hetzminoo_fancy( snpg, 'rich', hist_pars=list( plot=FALSE)) # still plots...
  lglk <- kinference::ilglk_geno( snpg, indiv_lglk_hist_pars=list( plot=FALSE))
  
  # Could use SPA to "Gaussianize"...
  stat <- cbind( hmfr, lglk)
  mdull <- colMeans( stat[ boring,]) # mean() doesnt' work column-wise, unlike var()...
  vdull <- var( stat[ boring, ])
  sstat[ f, s]:= stat[ f, s] - mdull[ s]
  
  # should do this with cholesky but the doco is SO SHITE
  iv <- solve( vdull)
  # dist <- sstat %*% (iv %*% t( sstat)) god knows what the bloody syntax is supposed to bloody be
  dist[ f]:= sstat[ f, s] %[s]% iv[ s, s1] %[s1]% sstat[ f, s1]
    
  hist( dist[ boring], nc=50)
  abline( v=dist[ focusees], col='red')

returnList( general_dist=dist[ boring], focus_dist=dist[ focusees])
} ## Crapometer imported from MVB2 24/9/18


#' drop_dups_pairwise_equiv(): Equivalence class sifter
#'
#' Constructs equivalence classes from pairwise equivalences, and returns
#' the "surplus" elements; if you then drop those elements, only one
#' element from each eq-class will be retained. Requires 2-col matrix
#' showing equivalent pairs. Code is taken from Numerical Recipes so I
#' should rewrite it perhaps (original algorithm is by Knuth).
#'
#' @param ij 2-column matrix or data.frame; probably "row numbers" in a
#'           dataset, though might work with character strings too
#' @param want_groups if \code{TRUE}, also return the equivalence-classes
#'                    themselves, as attribute \code{groups}.
#'
#' @return Surplus elements in \code{ij}, perhaps plus attributes \code{groups}
#'         if \code{want_groups=TRUE}. You can look at that to figure out which
#'         elements are being retained (one "representative" from each equiv
#'         class).
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


#' est_ALF_6way(): estimation of ALFs given 6-way genotypes and snerr
#'
#' Performs 'straight' estimation of ALFs, given 6-way genotypes and snerr. Won't
#' allow for changes in C-allele frequency from one population to the next. In
#' principle, should just use genocalldart::choose_geno6_thresholds but fix the
#' count-related thresholds and re-estimate ALFs.
#'
#' @param snpg a snpgeno object with 'snerr' included
#' @param control a param. Defaults to an empty list
#' @export

"est_ALF_6way" <-
function( snpg, control=list()) {
#### "Straight" estimation of ALFs given 6way genotypes and precalced snerr
## This won't allow for changes in C-allele freq from one popn to the next
## In principle, should use choose_geno6_thresholds but fix the count-related thresholds and just re-estimate ALFs

  define_genotypes()
  nl <- ncol( snpg)
  nf <- nrow( snpg)
  n <- matrix( 0, nl, 6, dimnames=list( NULL, genotypes6))
  for( ig in genotypes6) {
    n[ , ig] <- colSums( snpg==ig)
  }
  pobs <- 0 * n[1,]

  extract.named( snpg$locinfo[ cq( snerr, pbonzer)])
  A <- 1L
  B <- 2L
  O <- 3L

  perr <- structure( 1:4, names=colnames( snerr))
  extract.named( as.list( perr)) # AA2AO=1  etc


  lglk_l <- function( ppar) {
      p6[ A] <<- 0.9999 * inv.logit( ppar[ 1]) # avoid OOR
      p6[ B] <<- 0.9999 * (1-p6[A]) * inv.logit( ppar[ 2])
      p6[ O] <<- 1 - p6[ A] - p6[ B]
      pAO <- 2 * p6[ A] * p6[ O]
      pAA <- sqr( p6[ A])
      pBO <- 2 * p6[ B] * p6[ O]
      pBB <- sqr( p6[ B])
      pAB <- 2*p6[A] * p6[B]
      pOO <- max( 0, 1-pAA-pAB-pAO-pBB-pBO)
      pobs[ AA] <<- pAA * (1-perr[ AA2AO]) + pAO * perr[ AO2AA]
      pobs[ AO] <<- pAA * perr[ AA2AO]     + pAO * (1-perr[ AO2AA])
      pobs[ AB] <<- pAB
      pobs[ BB] <<- pBB * (1-perr[ BB2BO]) + pBO * perr[ BO2BB]
      pobs[ BO] <<- pBB * perr[ BB2BO]     + pBO * (1-perr[ BO2BB])
      pobs[ OO] <<- pOO
    return( n_l %**% log( pobs))
  }

  Nlglk <- NEG( lglk_l)

  for( il in 1:nl) {
    if( il %% 10==1) {
      cat( il, '\r'); flush.console()
    }
    perr[] <- snerr[ il,]
    n_l <- n[ il,]
    p6 <- rep( -1, 3) # c( A=-1, B=-1, O=-1)
    pstart <- c( pbonzer[ il, 'A'], pbonzer[ il, 'B'])
    pstart <- logit( c( pstart[ 1], pstart[2] / (1-pstart[ 1]) ) )
    fitto <- nlminb( pstart, Nlglk, control=control)
    pbonzer[ il,] <- c( p6[1:2], 0, p6[ 3])
  }

  snpg$locinfo$pbonzer <- pbonzer
return( snpg)
}  ## est_ALF_6way imported from MVB2 24/9/18



#' @rdname find_POPs
#' @export
#' @importFrom gbasics sqr
#' @importFrom atease @
#' @importFrom mvbutils cq %upto% %that.are.in% my.all.equal %without.name%

"find_duplicates" <-
function(snpg, subset1=1 %upto% nrow( snpg),
         subset2=subset1, max_diff_genos, keep_n=0.5*nrow(snpg)){

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

    # Newer versions don't bother defining pure-R wrappers for .Call-ees
#    result <- DUP_paircomps_lots(
#        temp_snpg,
#        temp_snpg,
#        TRUE,
###        max_diff_genos)  ## uncomment result<- to this close, for MVB version
#        geno1= temp_snpg,
#        geno2= temp_snpg,
#        symmo= TRUE,
#        max_diff_genos = max_diff_genos
#      ) ## MVB2 version of 'result', 24/9/18

    result <- DUP_paircomps_lots(
        geno1= temp_snpg,
        geno2= temp_snpg,
        symmo= TRUE,
        max_diff_genos = max_diff_genos,
        keep_n = keep_n
      ) ## DLM version, 24/9/18
#  } else { # different subsets
#    result <- DUP_paircomps_lots(
#        temp_snpg[ , subset1],
#        temp_snpg[ , subset2],
#        FALSE,
#        max_diff_genos)
#  }  ## MVB2 version of 'else result', 24/9/18

  } else { # different subsets
    result <- DUP_paircomps_lots(
        geno1= temp_snpg[ , subset1],
        geno2= temp_snpg[ , subset2],
        symmo= FALSE,
        max_diff_genos = max_diff_genos,
        keep_n = keep_n
      )
  }  ## DLM version, 24/9/18

  # just return the data.frame with 3 columns, everything else goes in
  # the attributes
  result <- with( result, data.frame( ndiff=big_similar, i=big_i, j=big_j))
  result@call <- sys.call()

  # warning if we're running up against storage constraints
  if(length(result$ndiff) == keep_n){
    warning("Number of returned duplicates equals keep_n. There may be more than keep_n duplicates. Increase keep_n to make sure you have them all!")
  }

return( result)
}


#' find_dups_with_missing(): Bare documentation
#'
#' This function has only the bare minimum of documentation necessary for roxygen to
#' parse it. We should probably add some proper documentation here. From the name,
#' it seems like this function might do a similar thing to find_duplicates(), but with
#' special handling for missing data. The structure of the code is certainly similar.
#'
#' @param snpg a param
#' @param subset1 a param, default 1 %upto% nrow(snpg), where %upto% generates an integer
#'                sequence by counting upwards, and hence will not generate c(1,0) if
#'                nrow(snpg) == 0.
#' @param subset2 a param. Defaults to subset1.
#' @param max_diff_ppn a param
#' @param limit a param. Default 10000.
#' @export

"find_dups_with_missing" <-
function( snpg,
  subset1= 1 %upto% nrow( snpg),
  subset2= subset1,
  max_diff_ppn,
  limit= 10000) {
  # Sanity...
stopifnot( is.numeric( subset1) && is.numeric( subset2))
stopifnot( all( !duplicated( subset1)) && all( !duplicated( subset2)))
stopifnot( my.all.equal( subset1, subset2) || !length( intersect( subset1, subset2)))
stopifnot( all( c( subset1, subset2) %in.range% (1:nrow( snpg))))
stopifnot( !missing( max_diff_ppn))

  og <- options( vecless.print=FALSE)
  on.exit( options( og))

  define_genotypes()

  temp_snpg <- snpg
  if( my.all.equal( snpg@diplos, genotypes6)){
    temp_snpg@diplos <- genotypes4_ambig
    temp_snpg[ snpg==AO] <- AAO
    temp_snpg[ snpg==AA] <- AAO
    temp_snpg[ snpg==BO] <- BBO
    temp_snpg[ snpg==BB] <- BBO
    temp_snpg[ snpg==OO] <- OO # need to do OO & AB too, since codes are different in 4way vs 6way
    temp_snpg[ snpg==AB] <- AB
  } else if( !my.all.equal( snpg@diplos, genotypes4_ambig)){
stop( "Don't have yet code to handle this set of allowed genos")
  }

  # Sort loci to improve speed (highest chance of mismatch first)
  # so that exit on non-dups is quicker
  gtab <- matrix( 0, 3, ncol( snpg), dimnames=list( genotypes4_ambig %except% OO, NULL))
  for( ig in genotypes4_ambig %except% OO) {
    gtab[ ig,] <- colSums( temp_snpg==ig)
  }
  gtab <- gtab / rep( colSums( gtab), each=3)
  pid <- colSums( sqr( gtab))
  o <- order( pid)
  pid <- pid[ o]
  temp_snpg <- temp_snpg[ ,o]

  # Remove extranea and recode with NAs
  OO_code <- which( temp_snpg@diplos==OO)
  attributes( temp_snpg) <- attributes( temp_snpg)[ 'dim'] # unclass but more drastic
  if( OO_code != 0) {
    temp_snpg[ temp_snpg==as.raw( 0)] <- as.raw( 9) # safe...
    temp_snpg[ temp_snpg==OO_code | temp_snpg==NA_geno] <- as.raw( 0)
  } else {
    temp_snpg[ temp_snpg==NA_geno] <- as.raw( 0)
  }

  temp_snpg <- t( temp_snpg)

  # Trying special-cases here to minimize copying
  if( my.all.equal( subset1, subset2)) {
    if( !my.all.equal( subset1, 1 %upto% ncol( temp_snpg))) {
      temp_snpg <- temp_snpg[, subset1]
    }

    result <- DUP_paircomps_incomplete_lots(
        geno1= temp_snpg,
        geno2= temp_snpg,
        symmo= TRUE,
        max_diff_ppn= max_diff_ppn,
        limit= limit)
  } else { # different subsets
    result <- DUP_paircomps_incomplete_lots(
        geno1= temp_snpg[ , subset1],
        geno2= temp_snpg[ , subset2],
        symmo= FALSE,
        max_diff_ppn= max_diff_ppn,
        limit= limit)
  }

  if( result %is.not.a% 'list') {
stop( sprintf( 'Hit limit=%i dups by %i-th sample; aborting', limit, result))
  }

  # just return the data.frame with 3 columns, everything else goes in
  # teh attributes
  result <- with( result, data.frame( ppn_diff= big_ndiff / big_ncomp,  i=big_i, j=big_j, ndiff=big_ndiff, ncomp=big_ncomp))
  result@call <- sys.call()

return( result)
} ## find_dups_with_missing ported in during merge 24/9/18



#' @rdname find_FSPs_from_POPs
#' @export
"find_FSPs_from_HSPs" <-
function( snpg, candiHSPs) {
  # For pairs already picked as HSPs, ie PLOD(HSP,UP) > eta: they might be FSPs

  # Don't need full pairwise screening for FSPs (do post hoc on a few hundred
  # HSPs), hence all in R.

  define_genotypes()

  # HSPs normally from 'find_HSPs'; or can be M*2 matrix of rows in snpg that are poss HSPs
  # if former, make latter

  if( candiHSPs %is.a% 'data.frame') {
    candiHSPs <- as.matrix( candiHSPs[ cq( i, j)])
  }
  sibg <- just_sibg <- snpg[ c( candiHSPs),]

  # Transform to 4way genotypes
  # based on code in find_duplicates
  # careful, since "factor level" of AB and OO is different in 4way vs 6way
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
  # Need lookups into the compressed matrix LOD4 (which should really be 3D array but Rcpp can't cope poor baby)
  # mg <- OLOD@mg # doesn't exist; lost when LOD4 gets plonked into locinfo data.frame

  mg <- make_genopairer( genotypes4_ambig)

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

  # These DO NOT sum to 1 by locus, because G1G2 and G2G1 are doubled-up
  # They would, if g1 & g2 were sorted
  # It doesn't matter for calculating *observed* PLODs, but care is needed with the expectations below...

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


  nsib <- nrow( candiHSPs)
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
  evalq(for(i in 1:nsib){
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
  })
  OD_FH <- p12fsp/p12hsp
  LOD_FH <- log(OD_FH)
  PLOD_FH <- colSums(LOD_FH)

  # Expectations
  # Need sum-to-1 here, so either rewrite when vecless2 appears, or work with compressed forms...
  OPHSP <- exp( OLOD) * OPUP # Pr[gg|HSP] <- 0.5 * PUP4 + 0.5 * Pr[gg|kappa=1]
  P_k0 <- OPUP
  P_k1 <- 2*OPHSP - OPUP
  P_k1[ P_k1 < 0] <- 0 # rounding error
  P_k2 <- 0 * P_k0 # get the shape right
  P_k2[ , diag( mg)] <- sqrt( P_k0[ , diag( mg)]) # only the cases where g1==g2

  p12fspa <- kappa_fsp[1] * P_k0 +  kappa_fsp[2] * P_k1 + kappa_fsp[3] * P_k2
  p12hspa <- kappa_hsp[1] * P_k0 + kappa_hsp[2] * P_k1 + kappa_hsp[3] * P_k2

  EPLOD_FH_F <- sum(log(p12fspa/p12hspa) * p12fspa)
  EPLOD_FH_H <- sum(log(p12fspa/p12hspa) * p12hspa)

  # format a return object
  ret <- data.frame(PLOD_FH = PLOD_FH,
                    i       = candiHSPs[,1],
                    j       = candiHSPs[,2])

  # Next 2 are COMPLETELY WRONG !!!
##  ret@E_FSP <- EPLOD_FH_F
##  ret@E_HSP <- EPLOD_FH_H

  ret@EPLOD_FH_F <- EPLOD_FH_F
  ret@EPLOD_FH_H <- EPLOD_FH_H

  ret@call <- sys.call()

  return(ret)
}  ## version taken from MVB2 at 24/9/18




#' Find full-sib pairs among parent-offspring pairs or among half-sib pairs
#'
#' For pairs already picked as likely parent-offspring pairs (POPs),
#' i.e., those with a weighted pseudo-exclusion (WPSEX) statistic
#' less than some threshold, they might be full sibling pairs (FSPs).
#' This function checks potential POPs with very low WPSEX values for
#' their potential to be FSPs.
#' 
#' The general idea of find_FSPs_from_POPs() is that pairs which are
#' _either_ POPs _or_ FSPs should stand out very clearly from
#' everything else, via find_POPs(). Then the job is to pick
#' between those possibilities. The workflow is supposed to be:
#' 
#' \itemize{
#' \item nail POPs/FSPs first with find_POPs()
#' \item pick between them with find_FSPs_from_POPs() (update: this doesn't
#' work very well...)
#' \item look for HSPs and filter out already-known POPs and FSPs
#' }
#' 
#' Hence the other function, find_FSPs_from_HSPs(), is theoretically
#' unnecessary in that you have already run find_POPs() and
#' find_FSPs_from_POPs() so you should know which of your "HSPs" are
#' really something else. But, nevertheless it's handy to have.
#' 
#' Both functions return expected values under different possible kin-types
#' (not variances, since these cannot be predicted for all kin-types).
#' 
#' The statistic for find_FSPs_from_POPs() is based on the weighted sum
#' of the number of exactly-matching 4-way genotypes, with weights chosen to
#' have high power for this particular discrimination. Weighting is optimized
#' for the unlikely scenario that POPs and FSPs are equally likely a priori,
#' but in practice the weights are not sensitive to this. The test is
#' deliberately crude and robust--- e.g. it avoids exclusion-based checks---
#' on the assumption that you have enough loci to pick HSPs, so the
#' more-related kin-types should be slam-dunks. *But* it doesn't seem
#' powerful enough. More worked needed...
#' 
#' find_FSPs_from_HSPs() again uses 4-way genotypes only (to avoid having
#' to worry about errors) but in a properly optimal PLOD designed for FSP/HSP
#' discrimination--- its expectation is positive for FSPs and negative for
#' HSPs. Theoretical means for those are returned as attributes (variances
#' cannot be predicted). Haven't added means for POPs or UPs since you're not
#' "supposed" to have those in the mix by the time you run
#' find_FSPs_from_HSPs(), but maybe I should fix that at some point.
#' 
#' @aliases find_FSPs_from_POPs find_FSPs_from_HSPs find_FSPs_from_candiHSPs
#' @param snpg a \code{snpgeno} object
#' @param candiPOPs candidate kin-pairs--- normally, a dataframe
#'                  with rows being pairs and columns _i_ and _j_ (and possibly
#'                  others) e.g. from find_POPs() or find_HSPs(). Can also be a
#'                  2-column matrix (each row again one pair).
#' @param candiHSPs candidate kin-pairs--- normally, a dataframe
#'                  with rows being pairs and columns _i_ and _j_ (and possibly
#'                  others) e.g. from find_POPs() or find_HSPs(). Can also be a
#'                  2-column matrix (each row again one pair).
#' @keywords misc
#' @examples
#' # pops_or_fsps <- find_POPs( mysnpg, ...)
#' ## do histograms etc to find likely ones
#' # discro <- find_FSPs_from_POPs( mysngp, pops_or_fsps %where% (wpsex < 0.042))
#' # hist( discro, nc=30, col='grey')
#' # abline( v=discro@E_POP, col='red')
#' # text( discro@E_POP, par( 'usr')[4], 'POP', col='red', pos=1) # below
#' # abline( v=discro@E_FSP, col='lightblue')
#' # text( discro@E_FSP, par( 'usr')[4], 'FSP', col='lightblue', pos=1) # below
#' # abline( v=discro@E_HSP, col='pink', lty='dashed')
#' # text( discro@E_HSP, par( 'usr')[4], 'HSP', col='pink', pos=1) # below
#' ###
#' # h_or_f <- find_HSPs( mysnpg, ...)
#' ## do histograms etc to find likely sib pairs
#' # discro2 <- find_FSPs_from_HSPs(  mysngp, h_or_f %where% (PLOD > 55))
#' # hist( discro2$PLOD_FH, nc=20, col='grey') # HSPs "should" be < 0, FSPs > 0
#' # abline( v=discro2@E_HSP, col='orange')
#' # text( discro2@E_HSP, par( 'usr')[4], 'POP', col='orange', pos=1) # below
#' # abline( v=discro2@E_FSP, col='lightblue')
#' # text( discro2@E_FSP, par( 'usr')[4], 'FSP', col='lightblue', pos=1) # below
#' @export

"find_FSPs_from_POPs" <-
function( snpg, candiPOPs) { 
## Don't need full pairwise screening for FSPs (do post hoc on a few hundred
## candidate POPs), hence all in R.

  define_genotypes()

  # 'candiPOPs' normally from 'find_POPs'; or can be M*2 matrix of rows in snpg that are poss POPs
  # if former, make latter

  if( candiPOPs %is.a% 'data.frame') {
    candiPOPs <- as.matrix( candiPOPs[ cq( i, j)])
  }
  snpg <- snpg[ c( candiPOPs),]

  # Transform to 4way genotypes
  # based on code in find_duplicates
  # careful, since "factor level" of AB and OO is different in 4way vs 6way


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

  # Exactly the same as version below, but this one seems less clear  
  # Pr_same_given_k[,{off+1}] <- pA * (1-2*pB*(pA+pO)) +
  #    pB * (1-2*pA*(pB+pO)) +
  #    pO * ( pA*pA + pB*pB + pO*pO)

  Pr_same_given_k[,{off+1}] <- pA * (sqr(pB) + sqr(1-pB)) +
      pB * (sqr(pA) + sqr(1-pA)) +
      pO * ( pA*pA + pB*pB + pO*pO)
      
  Pr_same_given_k[,{off+2}] <- 1
  Pr_same_given_k[ l, {off+0}] := pgeno[l,g] %[g]% pgeno[l,g]

  Pr_nsame_FSP <- c( 1/4, 1/2, 1/4)
  Pr_same_FSP[ l]:= Pr_nsame_FSP[ k] %[k]% Pr_same_given_k[ l, k]
  Pr_same_POP[ l]:= Pr_same_given_k[ l, {off+1}]

  Pr_nsame_HSP <- c( 1/2, 1/2, 0) # might as well...
  Pr_same_HSP[ l]:= Pr_nsame_HSP[ k] %[k]% Pr_same_given_k[ l, k]

  SD_FSP <- sqrt( Pr_same_FSP * (1-Pr_same_FSP))
  SD_POP <- sqrt( Pr_same_POP * (1-Pr_same_POP))
  SDwt_POP <- 0.5 # hard-wire for "prior" of POPs and FSPs equally likely. Doesn't matter.
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

  ret <- data.frame(wtsame = wtsame,
                    i       = candiPOPs[,1],
                    j       = candiPOPs[,2])

  ret@E_FSP <- wt %*% Pr_same_FSP
  ret@E_POP <- wt %*% Pr_same_POP
  ret@E_HSP <- wt %*% Pr_same_HSP
  ret@E_UP <- wt %*% Pr_same_given_k[,off]
  
return( ret)
}


#' @rdname find_POPs
#' @importFrom atease @	
#' @importFrom mvbutils cq %upto% %that.are.in% my.all.equal extract.named %without.name%
#' @export

"find_HSPs" <- ## from DLM
function( snpg, subset1=1 %upto% nrow( snpg), subset2=subset1,
    one_in_X_eta, 
    keep_n=100000,
    eta= NULL,
    keep_thresh= NULL,
    nbins= 50,
    bins= NULL
    ) {
## snpg should have been thru 'prepare_PLOD_SPA' so it has @PPS
stopifnot( 'Kenv' %in% names( attributes( snpg)))

  # Sanity...
stopifnot( is.numeric( subset1) && is.numeric( subset2))
stopifnot( all( !duplicated( subset1)) && all( !duplicated( subset2)))
stopifnot( my.all.equal( subset1, subset2) || !length( intersect( subset1, subset2)))

  keep_thresh_set <- !is.null( keep_thresh)
  if( !keep_thresh_set) {
      keep_thresh <- -1e23 # no PLOD will get that far!!

           }
    
  og <- options( vecless.print=FALSE)
  on.exit( options( og))

  # Here I'm using L-R tail approx SPA for CDF although Kenv$inv_CDF is likely
  #  more accurate for "moderate" tails but I don't quite trust it in the
  #  extremes actually they are pretty similar
  # Possibly, Kenv$inv_CDF should check if arg exceeds the range it was fitted
  #  to, and if so call inv_CDF_SPA2() instead but the range used in fitting is
  #  very goddamn wide (say +/- 10 SD) !

  define_genotypes()

  for( iwhat in cq( K, dK, ddK, inv_CDF)) {
    assign( iwhat, snpg@Kenv[[ iwhat]])
  }
  set_thresholds( keeping='hi')

  # For 4way loci, temporarily treat XO as XX...
  # ... have already adjusted the LOD entries so that new_LOD6( XX/..) <- LOD4( XXO/..)
  # ... use the LOD that's in Kenv, where SPA is calculated

  extract.named( snpg@locinfo[ cq( useN, LOD6, LOD4)])
##  use4 <- !use6
  temp_snpg <- snpg
  recode4to6temp <- function( x) { x[ x=='AO'] <- AA; x[ x=='BO'] <- BB; x}
  recode3to6temp <- function( x) { x[ x=='AO'] <- AA; x[ x=='BO'] <- BB; x[ x=='OO'] <- BB; x}
##   recode3to6temp <- function( x) { x[ x=='AAO'] <- AA; x[ x=='BBOO'] <- BB; x}
    ##  temp_snpg[ , useN] <- recode4to6temp( snpg[, useN]) # (AA,AO) -> AA; (BB,BO) -> BB
    ##  temporarily commented out to cause 3-way genotyping

  temp_snpg[ , useN == 4] <- recode4to6temp( snpg[, useN == 4]) # (AA,AO) -> AA; (BB,BO) -> BB
  temp_snpg[ , useN == 3] <- recode3to6temp( snpg[, useN == 3]) # (AA,AO) -> AA; (BB,BO) -> BB
  temp_LOD <- snpg@Kenv$LOD # already done in prepare_PLOD_SPA

  # Remove extranea
  attributes( temp_snpg) <- attributes( temp_snpg)[ 'dim']
  temp_snpg <- t( temp_snpg)

  if( is.null( bins)) {
    qq <- (2:nbins-1)/nbins
    bins <- snpg@Kenv$inv_CDF( qq)
  }
  binprobs <- snpg@Kenv$CDF( bins)

  # Trying special-cases here to minimize copying
  if( symmo) {
    if( !my.all.equal( subset1, 1 %upto% ncol( temp_snpg))) {
      temp_snpg <- temp_snpg[, subset1]
    }

    xresult <- HSP_paircomps_lots(
        pair_geno= temp_LOD@mg,
        LOD= t( temp_LOD),
        geno1= temp_snpg,
        geno2= temp_snpg,
        symmo= TRUE,
        eta= eta,
        min_keep_PLOD= keep_thresh,
        keep_n = keep_n,
        bins= bins
      )
  } else { # different subsets
    xresult <- HSP_paircomps_lots(
        pair_geno= temp_LOD@mg,
        LOD= t( temp_LOD),
        geno1= temp_snpg[ , subset1],
        geno2= temp_snpg[ , subset2],
        symmo= FALSE,
        eta= eta,
        min_keep_PLOD= keep_thresh,
        keep_n = keep_n,
        bins= bins
      )
  }

  # warning if we're running up against storage constraints
if( keep_thresh_set && (length(result$big_PLOD) == keep_n)){
    warning("Number of returned HSPs equals keep_n, increase keep_n to ensure you got them all")
}  ## This warning will probably need to be modified in find_POPs, find_duplicates, etc. as well.

  result <- with( xresult, data.frame( PLOD=big_PLOD, i=big_i, j=big_j))
  attributes( result) <- c( attributes( result),
      xresult %without.name% cq( big_PLOD, big_i, big_j))

  # assign extra info as attributes
  result@mean_theory <- snpg@Kenv$dK( 0)
  result@var_theory <- snpg@Kenv$ddK( 0)
  result@mean_HSP <- snpg@Kenv$dK( 0) + sum(snpg@locinfo$Ediff) ## 
  result@bins <- bins
  result@binprobs <- binprobs
  result@eta <- eta
  result@keep_thresh <- keep_thresh
  result@call <- sys.call()

return( result)
}


#' Kin-finders for loads-of-SNPs datasets
#'
#' @aliases find_POPs find_HSPs find_duplicates
#'
#' @description
#' These take a \code{snpgeno} dataset that has been processed as
#' far as \code{check6and4} (and for HSPs, \code{prepare_PLOD_SPA})
#' and find various relations between the samples. Relationships
#' include duplicates (DUPs/dupes; \code{find_duplicates}),
#' parent-offspring pairs (POPs; \code{find_POPs_vs}) and half-sibling
#' pairs (HSPs; \code{find_HSPs}). Unrelated pairs are referred to as
#' UPs. One can specify the same or different subsets of the
#' \code{snpgeno} for comparison: e.g., first subset for the adults,
#' second for the juveniles.
#'
#' @section Kinformation:
#' The idea is that kin-finding is based on a statistic and a threshold
#' \code{eta}, where the latter is chosen to keep false-positives down
#' to a user-specified level. Anything "beyond" \code{eta} will be
#' treated as a kin-pair ("beyond" depends on how the statistic is
#' defined, i.e. whether a kin-pair should come out very low or very high).
#' However, you're also likely to want to look post hoc at the distro of
#' computed statistics \emph{near} \code{eta}, to see whether separation
#' is as clean (or otherwise) as expected--- and also very unbeyond
#' \code{eta} into the zone where UPs are entirely dominant, to check that
#' theory is OK. So, as well as returning the "interesting" pairs that have
#' a statistic close to or on the non-UP size of \code{eta}, the POP and HSP
#' versions also return \emph{summaries} of the distribution of the
#' statistic. The thing is that there will be zillions of statistics from
#' UPs--- enough to blow out computer memory--- and they are not individually
#' interesting. Specifically, the main things returned are:
#'
#' \itemize{
#' \item mean and variance of stats. Computation is restricted to those on
#' the UP-side of \code{eta} (which is nearly all of them, usually) in order
#' to avoid distortion from non-UP cases. The latter will often be so rare
#' that distortion would be negligible--- but means and variances are not
#' "robust", so . Almost all will be include
#' \item counts of binned stats, regardless of whether above or below
#' \code{eta}. The bins are set based on SPAs to the theoretical
#' distributions, and chosen so that an equal number of UP-pairs should
#' fall into each bin.
#' \item cases where the stat is "interesting", i.e. on the non-UP side of
#' \code{keep_thresh}, as a \code{data.frame}. See \bold{Value} for details
#' }
#'
#' The process is controlled by three numbers: \code{nbins} for number of
#' bins, \code{eta} itself, and some nearby threshold \code{keep_thresh}
#' on the UP-side of \code{eta} (it will be automatically set to \code{eta}
#' otherwise) to determine which pairs are explicitly retained for your
#' inspection. There are two ways to specify \code{eta} and \code{keep_thresh}.
#' Usually, you would start with the indirect method, where you choose the
#' predicted-false-positive proportion of UP-pairs via the parameter
#' \code{one_in_X_eta}, and \code{rough_n_pairs_to_keep}. The routines then
#' use SPAs to the corresponding values of \code{eta} and \code{keep_thresh};
#' the returned value of \code{eta} is what you can subsequently use to make
#' the actual kin-decisions yourself after the event (by subsetting the
#' "interesting" pairs, comparing the statistic for each pair to
#' \code{eta})--- assuming that observed does match expected.
#'
#' But, sometimes it doesn't. In that case, the predicted values of \code{eta}
#' and \code{keep_thresh} may be way off the mark, and lead to retaining faaar
#' too few or too many pairs. If so, then look at the histogram of retained
#' statistics from an initial run, and try setting \code{eta} and/or
#' \code{keep_thresh} manually, rather than futzing around with the indirect
#' parameters until you get what you were after.
#'
#' @section Duplicates:
#' You have to set the retention threshold manually, via \code{max_diff_genos}
#' (see arguments). Post-processing step needed to to get the indices to remove
#' -- use \code{\link{drop_dups_pairwise_equiv}}, see \bold{Examples}.
#'
#' To avoid running \code{find_duplicates} on large numbers of fish at once,
#' one can split the dataset; see \bold{Examples}. You first need to run on each
#' subset separately (avoiding a quadratic number of comparisons) and reduce it
#' to non-duplicates (again, see \code{\link{drop_dups_pairwise_equiv}}), then
#' check the pair of reduced subsets (this will compare everything in the first
#' to everything in second, as the subsets are different). Note that when the
#' subsets are different, comparisons are made only \emph{between} subsets, not
#' \emph{within} each subset.
#'
#' Uses 4-way genotyping only, since these should be largely error-free. (Looks
#' like the exceptions are from samples with dodgy DNA.)
#'
#' @section Parent-offspring pairs:
#' 4-way genotypes are used to find "pseudo-exclusions" of the form AAO/BBO,
#' which \emph{usually} means AA/BB or AO/BB or AA/BO (a true exclusion), but
#' \emph{could} mean AO/BO (not an exclusion).
#'
#' \code{find_POPs} merely counts these, for loci where \code{Pr[O]+Pr[C]<pOC_max}
#' in order to avoid excess noise from AO/BO cases.
#' \code{find_POPs} uses all loci (by default) but weights them semi-optimally
#' so that pseudo-exclusions from loci with high "false pseudo-exclusion
#' probability" (i.e., high \code{Pr[AO/BO|UP]}) count much less than ones from
#' loci with very low null rates, for which AAO/BBO almost certainly means AA/BB.
#' We call this "Weighted PSeudo-EXclusion" ("WPSEX").
#'
#' The case AB/OO is also a (non-pseudo) exclusion, but is rarer than AAO/BBO
#' (non-existent for loci without nulls, of course). The count of such cases is
#' included in the output for "interesting" pairs in \code{find_POPs}; see
#' \bold{Value}.
#'
#' POP-finding is based on 4-way genotypes (OO, AAO, BBO, AB) to avoid
#' complications from genotyping-error-rates, and uses pseudo-exclusions rather
#' than likelihood-ratios; the latter is very sensitive to
#' false-negative-exclusions arising from typing-error, or even from mutation
#' with so many SNPs. You can get round that by including estimates of
#' typing-error-rate, but that's not necessarily easy to estimate in advance
#' insofar as it applies per-locus to POPs.
#'
#' @section Speed:
#' These are written in C (\code{Rcpp}) for speed, but for big datasets they
#' might still be quite slow. In the first instance, I certainly wouldn't try
#' them on 20,000 fish at once; I'd try with say 1000 then if that's OK 5000
#' etc. Bear in mind that they can always be run on different subsets of the
#' data, and the results patched back together (results will not change by doing
#' that). If you can run jobs in parallel, that could help a lot.
#'
#' @param snpg a \code{snpgeno} object
#' @param subset1 numeric vector of which samples to use (not logical,
#'                not negative). Defaults to all of them. Iff subset1 and subset2
#'                are identical, only half the comparisons are done (i.e., not
#'                _i_ with _j_ then _j_ with _i_). Some sanity checks are done.
#' @param subset2 numeric vector of which samples to use (not logical,
#'                not negative). Defaults to all of them. Iff subset1 and subset2
#'                are identical, only half the comparisons are done (i.e., not
#'                _i_ with _j_ then _j_ with _i_). Some sanity checks are done.
#' @param WPSEX_UP_POP_balance (\code{find_POPs}) loci receive a weight which
#'                             is proportional to (difference in probability of
#'                             pseudo-exclusion between UP and POP) / (variance
#'                             of indicator of pseudo-exclusion). But, should this
#'                             be variance assuming UP or POP?
#'                             \code{WPSEX_UP_POP_balance} sets the balance; bigger
#'                             values make it more UPpity, so placing more emphasis
#'                             on avoiding false-positives (which is probably the
#'                             Right Thing To Do). 0.99 could be completely fine...
#'                             (but hopefully \code{WPSEX_UP_POP_balance} won't
#'                             affect the result much anyway.)
#' @param one_in_X_eta expected number of false-positive UPs you can tolerate.
#'                     Setting this to say \code{1e6} means you'd expect 1 per million
#'                     comparisons. Used to set the threshold \code{eta}, which is
#'                     returned automatically.
#' @param rough_n_pairs_to_keep For checking, you can set this to trap many more
#'                              high-scoring pairs than you expect there to "really"
#'                              be, say a few thousand (NB the number of pairs retained
#'                              won't exactly equal this). You can subsequently look at
#'                              the "lucky losers" with high but sub-'eta' stats, and
#'                              then filter them out yourself by applying a cutoff of
#'                              \code{eta}. If you leave \code{rough_n_pairs_to_keep}
#'                              at its default of NA, the trap will be set at \code{eta},
#'                              so the result will contain exactly the pairs you want.
#'                              Values above \code{eta} will always be kept, even if you
#'                              specify something silly for \code{rough_n_pairs_to_keep}.
#' @param eta,keep_thresh see \bold{Description}. Can specify either or both. These override
#'            \code{one_in_X_eta} and \code{rough_n_pairs_to_keep} respectively.
#' @param keep_n Integer. Defines the ?maximum number of candidate pairs to keep. Will provide
#'               a warning if the number of identified pairs equals keep_n.
#' @param nbins number of bins to group the stats from the sub-'eta' pairs into. The
#'              bins will be set at quantiles of the expected distribution for UPs.
#' @param max_diff_genos (\code{find_duplicates}) max number of discrepant 4-way
#'                       genotypes to tolerate in "identical" fish. Try increasing
#'                       this from say 10 upwards, and hopefully nothing much will
#'                       change (though at some point things will change a lot, as
#'                       you get into the non-duplicate bit of the distribution).
#'                       See \bold{Duplicates} for how to remove duplicates from
#'                       the data.
#' @param quick whether to "compile" the functions for SPA, which use the magic
#'              \code{:=} operator. It speeds up the SPA bit but almost all the time
#'              is spent on actual POP-finding...
#' @param bins binning for PLODs (we throw away ones outside the range and bin them
#'             according to this within)
#'
#' @return a \code{data.frame} with 3 columns: statistic (\code{PLOD} or
#'         \code{wpsex} or \code{ndiff} number of mismatching genotypes),
#'         \code{i} (index in \code{subset1} of the first pair-member),
#'         \code{j} (index in \code{subset2} of the second). Note that
#'         \code{i} and \code{j} refer to the \emph{subsets}, not to the
#'         rows of the original \code{snpg}. Note that, iff you have set
#'         \code{rough_n_pairs_to_keep}, these will include pairs below the
#'         FP cutoff (which is returned as \code{eta}).
#'
#'         \code{find_POPs} adds a column named \code{nABOO}, showing the
#'         number of AB/OO exclusions for that potential POP. This is a useful
#'         additional diagnostic; it should be close to 0 for true POPs (it
#'         can only result from genotyping error or mutation, whereas AAO/BBO
#'         can result from nulls). For UPs, I was seeing values typically in
#'         the low 20s, which is pretty good separation.
#'
#'         For duplicates, not \emph{all} pairwise duplicates are recorded,
#'         unless the subsets are different--- otherwise you could have
#'         quadratic horror of enormous numbers of pairs arising from a cluster
#'         of say 100 identical controls! Since "duplication" is transitive (ie
#'         if i & j are the same, and i & k are the same, then j & k must also
#'         be the same), only the necessary ones are recorded to allow you to
#'         filter out yourself afterwards. e.g., if samples 1, 3, 5, and 6 are
#'         all duplicates, you'll get this:
#'         %#
#'         \item{# without "ndiff" column}{}
#'         \item{  i j}{}
#'         \item{  3 1}{}
#'         \item{  4 3}{}
#'         \item{  6 4}{}
#'
#'         but you won't see the pairings for 1/4, 1/6, 3/6. If you just want
#'         to strip out all duplicates bar one in each group (and you don't
#'         care which one is kept), then you can use the function
#'         \code{\link{drop_dups_pairwise_equiv}} --- see \bold{Examples}.
#'
#'         For POPs and HSPs, the following are also returned as attributes
#'         (that can be accessed by \code{@} if \code{atease} is loaded). The
#'         main point is that the "boring" below-threshold pairs get put into
#'         bins and are not kept individually. The names sometimes change
#'         depending on which statistic is being used.
#'
#'         \item{eta}{false-positive cutoff to be applied to the statistic in
#'         question (automatically done if \code{rough_n_pairs_to_keep==NA},
#'         or up to you if not). Variance of the stat will only be calculated
#'         from values to the "UP side" of \code{eta}. However, the set of
#'         retained pairs/individuals is actually controlled by...}
#'         \item{keep_thresh}{the cutoff used to retain "interesting" pairs.
#'         Usually obvious from the range of statistic values.}
#'         \item{mean_sub_<stat>, var_sub_<stat>}{empirical values for the
#'         statistic when it is below \code{eta} (ie nearly always).}
#'         \item{mean_theory, var_theory}{of the statistic, to compare to
#'         previous.}
#'         \item{n_<stat>_in_bin}{number of pairs whose statstic fell within
#'         the range of each bin}
#'         \item{bins}{cutpoints for the bins. These should be quantiles,
#'         according to the SPA; so if practice matches theory, the
#'         numbers-per-bin should all be similar.}
#'
#' @examples
#' \dontrun{
#' ## duplicate checking. ckmini2 has 6 fish where 1,3,4,6 are all identical (zero differing loci).
#' ## there's only 7 crappy loci and I faked the data for this anyway, so strict identical is needed
#' ## All-in-one
#' #test <- find_duplicates( ckmini2, max=0) # strict identity
#' #test
#' ##  ndiff i j
#' ##1     0 3 1
#' ##2     0 4 3
#' ##3     0 6 4
#' ## To remove them--- subtlety of keeping ONE from each group
#' #droppies <- drop_dups_pairwise_equiv( test[,2:3])
#' #droppies # 1, 4, 6
#' #ckmini2_nodups <- ckmini2[ -droppies, ]
#' ## Two-stage
#' #first_half <- 1:3
#' #second_half <- (1:nrow( ckmini2)) \%except\% first_half
#' #test1 <- find_duplicates( ckmini2, subset1=first_half, subset2=first_half, max=0)
#' #test1
#' ##  ndiff i j
#' ##1     0 3 1
#' #droppies1 <- first_half[ drop_dups_pairwise_equiv( test1[,2:3])] # NB must do lookup in subset
#' #test2 <- find_duplicates( ckmini2, subset1=second_half, subset2=second_half, max=0)
#' #droppies2 <- second_half[ drop_dups_pairwise_equiv( test2[,2:3])] # 4
#' ## Now check 2nd half vs 1st
#' #test2_1 <- find_duplicates( ckmini2,
#' #    subset1=first_half \%except\% droppies1,
#' #    subset2=second_half \%except\% droppies2,
#' #    max=0)
#' ## Simpler since no internal checks. Just remove 2nd-halfers that match something in the 1st-half
#' #droppies2_1 <- (second_half \%except\% droppies2)[ test2_1[,'j']) # 6
#' #droppies <- c( droppies1, droppies2, droppies2_1)
#' #ckmini2_nodups2 <- ckmini2[ -droppies,]
#' ## HSPs: comparing everything with itself (not sensible for real data, should take out adults first)
#' ## set threshold for 1 FP
#' #test <- find_HSPs( ckdata, one_in_X_eta=sqr( nrow( ckdata))/2 )
#' ## POPs: Ad-Ju comps; again 1 FP
#' #test <- find_POPs( ckdata, subset1=adults, subset2=juves,
#' #    one_in_X_eta=length( adults) * length( juves), rough_n_pairs_to_keep=500)
#' }
#' @importFrom gbasics sqr
#' @importFrom atease @
#' @importFrom vecless := compile_vecless
#' @importFrom stats runif
#' @importFrom mvbutils cq %upto% %that.are.in% my.all.equal extract.named %without.name%
#' @export

"find_POPs" <-
function( snpg, subset1=1 %upto% nrow( snpg), subset2=subset1,
    one_in_X_eta, # die
    rough_n_pairs_to_keep= NA, # C code cutoff, but merge w/ keep_thresh
    eta= NULL, # DO NOT CALL THIS THIS, but stats cutoff value needs spec
               # or calc from E/V of UPs
    keep_thresh= NULL, # merge, C code return cutoff
    keep_n=0.5*nrow(snpg),
    nbins,
    quick=TRUE,
    WPSEX_UP_POP_balance=0.99) {
###################
  # Sanity...
stopifnot( is.numeric( subset1) && is.numeric( subset2))
stopifnot( all( !duplicated( subset1)) && all( !duplicated( subset2)))
stopifnot( my.all.equal( subset1, subset2) || !length( intersect( subset1, subset2)))

  # Decide based #apparent exclusions of AA/BB form, using 4way genos, though
  # ... it's really AAO/BBO so not a true exclu but
  # ... close among pop-loci
  # Sticking with 4way genos so that genotyping errors are low

  # Instead of Nexclu, uses a wted sum of "exclus" in 4way genos to max expected diff between POP and UP
  # This version doesn't allow for geno errors, but does realize that AO/BO could happen
  # Doesn't bother with OO-AB

  define_genotypes()
  extract.named( snpg@locinfo[ cq( useN, PUP4, pbonzer)])
  p0 <- pbonzer[,'O'] + pbonzer[,'C']
  pA <- pbonzer[,'A']
  pB <- pbonzer[,'B']

  # "Exclusion" whenever AAO & BBO, but this *could* be AO/BO
  # I'm calling them "ex" for now anyway, hence pex_BLAH
  pex <- cbind(
      POP= 2 * p0 * pA * pB,
      UP= 2 * (2*pA*p0 + pA*pA) * (2*pB*p0 + pB*pB)
    )

  # Want to keep the POP and UP means as far apart as possible on the scale of SDs...
  # but, which SD? Make it WPSEX_UP_POP_balance * SD[UP] + (1-WPSEX_UP_POP_balance) * SD[POP]

  # Optimal wt would depend on p0 and to some extent on pA
  # wt should be 1 if p0==0 and 0 if p0==1
  delta <- pex[,'UP'] - pex[,'POP'] # mathematically I think this *can't* be -ve
  SD <- sqrt( pex * (1-pex))
  SD_combo <- WPSEX_UP_POP_balance * SD[,'UP'] + (1-WPSEX_UP_POP_balance) * SD[,'POP'] # %*% c( WPSEX_UP_POP_balance, 1-WPSEX_UP_POP_balance)
  V_combo <- sqr( SD_combo)
  ww <- delta / V_combo # considerable algebra appears to show this is optimal
  ww <- c( ww / sum( ww) ) # else get 1-col matrix
stopifnot( all( ww>0))

  pop_loci <- which( ww > 0) # all of them, for now
  # pop_loci <- which( pbonzer[,'O'] + pbonzer[,'C'] < pOC_max)

  temp_snpg <- snpg[ , pop_loci]
  recode4to6temp <- function( x) { x[ x=='AO'] <- AA; x[ x=='BO'] <- BB; x}
  temp_snpg <- recode4to6temp( temp_snpg) # (AA,AO) -> AA; (BB,BO) -> BB

  # Remove extranea
  attributes( temp_snpg) <- attributes( temp_snpg)[ 'dim']
  temp_snpg <- t( temp_snpg)

  # Distro of #excl loci for UPs
  opphetz <- c( 'AAO/BBO', 'BBO/AAO') %that.are.in% colnames( PUP4)
  pex_up <- PUP4[ pop_loci, opphetz]
  rr <- pex_up / (1-pex_up)
  log1m_pexup <- log1p( -pex_up)

  K <- function( tt) {
      # vecless2 notation could 1-step this via SUM
      KK[ l, j]:= log1m_pexup[ l] + log1p( rr[ l] * exp( tt[ j] * ww[ l]))
      K[ j] := KK[ +., j]
    return( c( K)) # without the c(), you get a scalar xtensor, and trouble...
    }

  dK <- function( tt) {
      retw[l,j] := rr[ l] * exp( tt[ j] * ww[ l])
      dKK[ l, j] := ww[ l] * retw[ l, j] / (1+retw[ l, j]) # guessing this is more accurate tahn 1-1/1+x
      dK[ j] := dKK[ +., j]
    return( c( dK))
    }

  ddK <- function( tt) {
     retw[l,j] := rr[ l] * exp( tt[ j] * ww[ l])
      ddKK[ l, j] := sqr( ww[ l]) * retw[ l, j] / sqr( 1+retw[ l, j]) # guessing this is more accurate tahn 1-1/1+x
      ddK[ j] := ddKK[ +., j]
   return( c( ddK))
   }

  n_loci <- ncol( snpg)
  n_sim_check <- 1000
  Ktest <- function( tt) {
    x <- matrix( runif( n_sim_check * n_loci) < pex_up, n_loci, n_sim_check)
    ewx <- exp( x*ww*tt)
    colSums( log( ewx))
  }

  if( quick) {
    K <- compile_vecless( K(0))
    dK <- compile_vecless( dK(0))
    ddK <- compile_vecless( ddK(0))
  }

  symmo <- my.all.equal( subset1, subset2)
  inv_CDF <- renorm_SPA_cumul( K, dK, ddK)$inv_CDF

  set_thresholds( keeping='lo') # cf 'hi' for HSPs

  # Prepare for diagnostics of #excl
  # prolly not needed but HSP C code already does it
  qq <- (2:nbins-1)/nbins
  pciles <- inv_CDF( qq) # inv_CDF_SPA2 may struggle with LOWER tail...

  # Trying special-cases here to minimize copying
  if( symmo){
    if( !my.all.equal( subset1, 1 %upto% ncol( temp_snpg))) {
      temp_snpg <- temp_snpg[, subset1]
    }
    result <- POP_wt_paircomps_lots(
        geno1= temp_snpg,
        geno2= temp_snpg,
        w= ww,
        symmo= TRUE,
        eta= eta,
        max_keep_wpsex= keep_thresh,
        keep_n = keep_n,
        bins= pciles,
        AAO= match( 'AA', snpg@diplos), # NB NB: AO has been recoded to AA
        BBO= match( 'BB', snpg@diplos)
      )
  } else { # different subsets
    result <- POP_wt_paircomps_lots(
        geno1= temp_snpg[ ,subset1],
        geno2= temp_snpg[ ,subset2],
        w= ww,
        symmo= FALSE,
        eta= eta,
        max_keep_wpsex= keep_thresh,
        keep_n = keep_n,
        bins= pciles,
        AAO= match( 'AA', snpg@diplos), # NB NB: AO has been recoded to AA
        BBO= match( 'BB', snpg@diplos)
      )
  }

  # warning if we're running up against storage constraints
  if(length(result$big_wpsex) == keep_n){
    warning("Number of returned POPs equals keep_n, increase keep_n to ensure you got them all")
  }


  # construct the result
  result <- with( result, data.frame( wpsex=big_wpsex, i=big_i, j=big_j))

  # calculate nABOO, only for interesting pairs
  snpg_i <- snpg[ subset1[ result$i], pop_loci]
  snpg_j <- snpg[ subset2[ result$j], pop_loci]
  isABOO <- ((snpg_i==OO) & (snpg_j==AB)) + ((snpg_i==AB) & (snpg_j==OO))
  result$nABOO <- rowSums( isABOO)

  # probably uneccessary ?
  result <- result %without.name% cq( big_wpsex, big_i, big_j)

  # add extra info
  result@bins <- pciles
  result@eta <- eta
  result@keep_thresh <- keep_thresh
  result@n_loci <- length( pop_loci)
  result@mean_theory <- dK( 0)
  result@var_theory <- ddK( 0)
  result@call <- sys.call()

return( result)
}

#' @rdname chain_pairwise
#' @export

"get_chain" <-
function( thing, seed) {

# if length(unique(seed)) == length(seed) ???

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

#' Heterozygotes minus "OO" checking
#'
#' This test looks at whether the allele frequencies in a given fish seem
#' right, or if there are discrepancies due to (i) degraded DNA or (ii)
#' sample contamination.
#' 
#' @param snpg a \code{\link[gbasics]{snpgeno}} object
#' @param target which weighting should be used. \code{"rich"} detects
#'               contaminated data (there are too many heterozygotes) and
#'               \code{"poor"} detectes DNA degredation (there are too
#'               few heterozygotes).
#' @param hist_pars parameters to pass to \code{\link{hist}}
#' @param multhresh A param.
#' @param showPlot show the plot? Default TRUE
#' @keywords misc
#' @export

"hetzminoo_fancy" <-
function( snpg, target=c( 'rich', 'poor'), hist_pars=list(), multhresh=1, showPlot = TRUE) {
###################
  define_genotypes()
  extract.named( snpg@locinfo[ cq(  PUP4, pbonzer)]) ## used to also call use6; not used
  p0 <- pbonzer[,'O'] + pbonzer[,'C']
  pA <- pbonzer[,'A']
  pB <- pbonzer[,'B']

  v <- 2*pA*pB + sqr( p0) - sqr( 2*pA*pB-sqr( p0)) ## Pr(AB) + Pr(OO) - (Pr(AB) - Pr(OO))^2
  target <- match.arg( target)
  edash <- if( target=='rich') {
               (2*pA*p0+sqr(pA)) *  ## Pr(AA|AO)
                   (1-sqr(1-pB)) +  ## original; Pr(BB|AB|BO)
                   (2*pB*p0+sqr( pB)) *  ## Pr(BO|BB)
                   (1-sqr(1-pA)) +  ## modified; Pr(AA|AB|AO)
                   sqr( p0) *      ## Pr(OO)
                   (1-sqr( p0))    ## Pr(!OO)
    } else {
      edash <- 2*(pA*pB + pA*p0 + pB*p0) # Equals Pr(heterozygote) under 3-way HWE
    }

  # Numericals have led to single loci getting all the weight when edash and v are both tiny!
  # So, shrink v a tiny bit...
  msqrt_v <- median( sqrt( v))
  v <- sqr( sqrt( v) * 0.99 + msqrt_v * 0.01)
  
  ww <- edash / v
  ww <- c( ww / sum( ww) ) # else get 1-col matrix
  stopifnot( all( ww>0))

  use_loci <- which( ww > 0) # all of them, for now
  temp_snpg <- snpg[ , use_loci]
  recode4to6temp <- function( x) { x[ x=='AO'] <- AA; x[ x=='BO'] <- BB; x}
  temp_snpg <- recode4to6temp( temp_snpg) # (AA,AO) -> AA; (BB,BO) -> BB

  delta <- (temp_snpg==AB) - (temp_snpg==OO)
  whmo[ f]:= delta[ f, l] %[l]% ww[ l]

  # Null distro: P[S==1] = pAB; P[S==-1] = pOO; P[S==0] = 1-pAB - pOO
  # E[ exp( tt*S)] = (e(tt)+1) * pAB + 1 + (e(-tt) -1)* pOO

  pAB <- 2*pA*pB
  pOO <- sqr( p0)
  four_pab_poo <- 4*pAB*pOO
  compaboo <- 1 - pAB - pOO

  # now setup the functions for the SPA
  # use the above as shortcuts

  K <- function( tt) {
      # here a *1 inside exp omitted
      etwab[ l, j] := pAB[l] * exp( tt[j] * ww[ l])
      # here a *-1 inside exp omitted
      etwoo[ l, j] := pOO[l] * exp( -tt[j] * ww[ l])
      # compaboo computed above, P(not AB & not OO)*exp(0)
      KK[ j]:= SUM_ %[l]% log( compaboo[l] + etwab[ l, j] + etwoo[ l, j])
    return( c( KK)) # without the c(), you get a scalar xtensor, and trouble...
    }

  # derivative of K
  dK <- function( tt) {
      etwab[ l, j] := pAB[l] * exp( tt[j] * ww[ l])
      etwoo[ l, j] := pOO[l] * exp( -tt[j] * ww[ l])
      denom[ l, j] := compaboo[l] + etwab[l,j] + etwoo[l,j]
      dKK[ j] := ww[l] %[l]% ((etwab[l,j]-etwoo[l,j])/denom[l,j])
    return( c( dKK))
    }

  # 2nd derivative of K
  ddK <- function( tt) {
      etwab[ l, j] := pAB[l] * exp( tt[j] * ww[ l])
      etwoo[ l, j] := pOO[l] * exp( -tt[j] * ww[ l])
      denom[ l, j] := compaboo[l] + etwab[l,j] + etwoo[l,j]
      ddKK[ j] := sqr( ww[l]) %[l]% ( (compaboo[l] * (etwab[ l, j] + etwoo[l,j]) + four_pab_poo[ l]) / sqr( denom[l, j]) )
    return( c( ddKK))
    }

  K <- compile_vecless( K(0))
  dK <- compile_vecless( dK(0))
  ddK <- compile_vecless( ddK(0))

  dens_SPA <- renorm_SPA( K, dK, ddK, 'func'
                         ##, already_vectorized=TRUE ## SB: not used in
                         ## renorm_SPA as it exists (DLM version). May want to revisit.
                         )

  # optional graphics and/or user-specified outputs
  switch( mode( hist_pars),
    list = {
        hist_pars <- add_list_defaults( hist_pars,
            main=sprintf( '%s: multhresh=%5.2f', target, multhresh),
            xlim= range( whmo), # so cutoff lines show
            xlab='', nclass=50)
        if (showPlot) {
            lv <- do.call( 'hist', c( list( x=whmo), hist_pars))
            with( lv, lines( mids, diff( breaks) * dens_SPA( mids) * sum( counts), col='green'))
        }
        # abline( v=ncuts, col='red')
      },
    expression = eval( badhetz_hist_pars),
    NULL = NULL
  )

return( c( whmo))
}


#' hsp_power(): Preparation for kin-finding
#'
#' \code{hsp_power} computes locus-wide LOD and PUP tables for each possible
#' genopair, and adds them to the "locinfo" attribute. It also calculates mean
#' and variance of LOD, and those can be used to calculate "Ediff" and
#' "SEdiff". It's needed before \code{\link{hetzminoo_fancy}} and
#' \code{prepare_PLOD_SPA}.
#' 
#' \code{prepare_PLOD_SPA} prepares a \code{snpgeno} for "exact" (SPA)
#' calculations of the null distro of PLOD (ie for true UPs). to a
#' \code{snpgeno}, ready for . It's needed before \code{\link{find_HSPs}}.
#' 
#' The \code{use6} field of \code{lociar@locinfo} determines whether 6-way or
#' 4-way genotypes are assumed in calculating LOD tables.
#' 
#' \code{hsp_power} in particular needs a biiiig tidy-up. It's daft to store
#' LODs for only one specific kin; it'd be better to always calculate P1share
#' and P2share as well as P0share (which is PUP), and then compute
#' whatever-is-needed later on-the-fly.
#' 
#' @aliases hsp_power prepare_PLOD_SPA
#' @param lociar \code{snpgeno} objects with the necessary ingredients
#' @param geno6 \code{snpgeno} objects with the necessary ingredients
#' @param want_LOD_table can't think why you'd set this to FALSE
#' @param k target average kinship for LOD; 0.5 for HSPs, 0.25 for HTPs, etc.
#' @param n_pts_SPA_renorm self-explanatory, or leave alone if not! (Could
#'                         reduce if reeeeally slow.)
#' @return \code{snpgeno} object with augmented columns in "locinfo" attr.
#'         \code{prepare_PLOD_SPA} adds an environment attribute "Kenv" used by
#'         SPA calculations
#' @importFrom atease @ @<-
#' @importFrom vecless make_playback
#' @keywords misc
#' @examples
#' ## Need some examples!
#' @export

"hsp_power" <- function( lociar,
    want_LOD_table=TRUE, # T/F
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

  if( exists( 'genotypes4_ambig', inherits=FALSE)) { # TRUE unless overridden sneakily...
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

    ### bring in an s3 as well:
    map6to3 <- matrix( 0, 6, 3, dimnames=list( genotypes6, genotypes3_ambig))
    # AB & OO are OK; AAO should receive both AA and AO; etc
    mm <- match( genotypes6, substring( genotypes3_ambig, 1, 2), 0) # the "AA" bit of "AAO"...
    yup <- cbind( which(mm>0), mm[ mm>0])
    map6to3[ yup] <- 1
    mm <- match( genotypes6, substring( genotypes3_ambig, 2, 3), 0) # ... and the "AO" bit
    yup <- cbind( which(mm>0), mm[ mm>0])
    map6to3[ yup] <- 1
    mm <- match( genotypes6, substring( genotypes3_ambig, 3, 4), 0) # ... and the "OO" bit
    yup <- cbind( which(mm>0), mm[ mm>0])
    map6to3[ yup] <- 1

    # Really want g4p0[l,i,j] := map6to4[i,k6] %[k6]% g6p0[l,k6,m6] %[m6]% map6to4[m6,j]
    # ... but vecless can't presently handle multi-stages

    A[l,i,k] := g6p0[l,i,j] %[j]% map6to3[j,k]
    g3p0[l,i,j] := map6to3[ k,i] %[k]% A[l,k,j]

    A[l,i,k] := g6p1[l,i,j] %[j]% map6to3[j,k]
    g3p1[l,i,j] := map6to3[ k,i] %[k]% A[l,k,j]

    s3 <- predict_hsp_util( g3p0, g3p1, want_LOD_table, k=k)
    
    if( want_LOD_table) {
      li$LOD6 <- s6@LOD # matrix
      li$PUP6 <- s6@PUP
      li$LOD4 <- s4@LOD
      li$PUP4 <- s4@PUP
      li$LOD3 <- s3@LOD
      li$PUP3 <- s3@PUP
      s6@LOD <- s6@PUP <- s4@LOD <- s4@PUP <- s3@LOD <- s3@PUP <- NULL

      li$e0_6way <- s6$e0
      li$v0_6way <- s6$v0
      li$e1_6way <- s6$e1
      li$v1_6way <- s6$v1

      li$e0_4way <- s4$e0
      li$v0_4way <- s4$v0
      li$e1_4way <- s4$e1
      li$v1_4way <- s4$v1

      li$e0_3way <- s3$e0
      li$v0_3way <- s3$v0
      li$e1_3way <- s3$e1
      li$v1_3way <- s3$v1

      s6$e0 <- s6$v0 <- s6$e1 <- s6$v1 <- NULL
      s4$e0 <- s4$v0 <- s4$e1 <- s4$v1 <- NULL
      s3$e0 <- s3$v0 <- s3$e1 <- s3$v1 <- NULL
    }

    # li <- cbind( li, s6, s4)
    li[ names( s6)] <- s6 # instead of cbind--- this overwrites

    li[ li$useN == 4, names( s4)] <- s4[ li$useN == 4,] ## subs in 4-ways where useN == 4
    li[ li$useN == 3, names( s3)] <- s3[ li$useN == 3,] ## subs in 3-ways where useN == 3

  } else { # ... sneaky override, for non-ABCO systems
    # shouldn't really be called "...6" obvs
    s6$useN <- 6
    li[ names( s6)] <- s6 # instead of cbind--- this overwrites
  }

  lociar@locinfo <- li
return( lociar)
}


#' hsp_power2(): Kin-finding power for microhaplotyped loci
#'
#' This is a short-term fudge for checking HSP-finding power of a bunch of loci that
#' (i) can have as many haplotypes as you like, but (ii) have no errors or nulls.
#'
#' At some point in future, \code{kinference} might be changed so that it can handle
#' >2 alleles gracefully. But not yet. So for now this version does some ghastly
#' "live-hacking" of existing code for \code{\link{hsp_power}} to implement no-errors
#' no-nulls multi-allelic case. It will be hard to follow, so use \code{mtrace} if
#' you really want to see what's going on. The guts of the code is in
#' \code{\link{hsp_power}} and \code{predict_hsp_util}.
#' 
#' @aliases hsp_power2
#' @param lociar Usually, a matrix of allele frequencies (Locus * Alleles). Locus names
#'               are set from the rownames, or "L1", "L2" etc if there are no rownames.
#'               Allele names will be set to "A", "B", "C", etc, regardless of colnames;
#'               you do not have a choice there. Will be renormalized so rows sum to unity.
#'               NB \code{lociar} can also be a \code{snpgeno} object, as expected for
#'               \code{hsp_power}. If so, then the allele freqs are assumed to live in
#'               \code{lociar$locinfo$pbonzer}, and \bold{no} nulls or genotyping errors
#'               are allowed for; hence, for a DartCap "ABCO"-style dataset, \code{hsp_power}
#'               and \code{hsp_power2} will give \emph{different} answers.
#' @param want_LOD_table can't think why you'd set this to FALSE
#' @param k target average kinship for LOD; 0.5 for HSPs, 0.25 for HTPs, etc.
#' @return If \code{lociar} is an allele-frequency matrix, then you get a dataframe
#'         with one row per locus and columns "Ediff", "V.UP", and "sdiff". "Ediff" is
#'         "E[LOD|HSP] - E[LOD|UP]"; "V.UP" is "V[LOD|UP]"; "sdiff" is \code{sqrt(V.UP)/Ediff},
#'         useful for ranking locus power. See \bold{Examples} for use.
#' @seealso hsp_power
#' @importFrom atease @ @<-
#' @importFrom vecless make_playback
#' @keywords misc
#' @examples
#' # ALF <- matrix( runif( 15), 3, 5) # 3 loci; 5 alleles
#' # POW <- hsp_power2( ALF)
#' # look at the contents of each...
#' # Now do it for lots of loci
#' # lots <- 500
#' # ALF <- matrix( runif( lots*5), lots, 5)
#' # POW <- hsp_power2( ALF)
#' # Now say we plan 1e6 pairwise comps, and might expect 100 HSPs
#' # Work relative to E[LOD|UP] which is not returned explicitly; treat that as "origin" ie 0
#' # V <- sum( POW$V.UP)  # V[PLOD|UP]
#' # E <- sum( POW$Ediff) # E[PLOD|HSP] - E[PLOD|UP]
#' # E / sqrt( V) # 10.5 SDs--- good ! Mean of HSPs is 10.5 UP-SDs above mean of UPs,
#' # ... so v. unlikely an UP will get as far as _typical_ HSP. But we need to be a bit
#' # ... more stringent than "typical"
#' # bigUP <- qnorm( 1e-6, mean=0, sd=sqrt( V), lower=FALSE) # most-kinlike UP
#' # smallHSP <- qnorm( 1e-2, mean=E, sd=sqrt( 4*V))      # least-kinlike HSP
#' # ... so that's probably OK...
#' # this is a test.
#' @export

hsp_power2 <- function( lociar,
    want_LOD_table=TRUE, # T/F
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

  if( exists( 'genotypes4_ambig', inherits=FALSE)) { # TRUE unless overridden sneakily...
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

    # li <- cbind( li, s6, s4)
    li[ names( s6)] <- s6 # instead of cbind--- this overwrites
    li[ names( s4)] <- s4 # instead of cbind--- this overwrites

    li[ !li$use6, names( s4)] <- s4[ !li$use6,]
  } else { # ... sneaky override, for non-ABCO systems
    # shouldn't really be called "...6" obvs
    s6$use6 <- TRUE
    li[ names( s6)] <- s6 # instead of cbind--- this overwrites
  }

  lociar@locinfo <- li
return( lociar)
}


#' Check individual genotypes for aggregate typicality
#'
#' \code{ilglk_geno} Computes log-likelihood of entire 4-way genotype of
#' each individual, i.e., sum log Pr[ g(i,l)]; and compare the distro of
#' log-likelihoods across individuals with its predicted shape given allele
#' frequencies. Significant mismatch is bad. Can also detect outliers. Up
#' to you what criterion to use for that.
#'
#' #' \code{lglk_loci} compares, for each locus, the average (across individuals)
#' observed lglk with the theoretical mean and variance. The idea is to help
#' figure out when some loci are going wrongish (e.g. you can get decent fits
#' from a subset of loci). Of course, \code{check6and4} pvals should be the
#' main guide here; \code{lglk_loci} can show an overall deviation, as well as
#' any remaining locus-specific misbehaviour (but shouldn't be much
#' locus-specific stuff thx2 \code{check6and4}). Overall
#' too-good-to-be-true-ism (as seen for \code{Glyphis garricki}) \emph{might}
#' come when ALF is estimated from very small datasets.
#'
#' You can use \code{locator(1)} to click the histogram to figure out
#' where to adjust the \code{xlim}/\code{ylim} values to change the range
#' of the data to inspect more closely.
#'
#' Currently, the SPA calcs are a wee bit slow because of heavy use of
#' \code{vecless} which in version 1.0 is sluggish. The lglks themselves
#' are computed in C and are blisteringly fast.
#'
#' Haven't added any formal uh-oh criteria yet; that could be done via the
#' SPA, as in \code{dump_badhetz_fish}. But, reading off from the graph is
#' probably fine...
#'
#' @aliases ilglk_geno lglk_loci
#' @param snpg a \code{\link[gbasics]{snpgeno}} (6-way genotype)
#' @param indiv_lglk_hist_pars list like in \code{dump_badhetz_fish}, for
#'                             controlling histogram
#' @param quick Controls whether K, dK, and ddK are sent through compile_vecless().
#'              Defaults to TRUE.
#' @param showPlot show the histogram? Defaults to TRUE
#' @return Vector of log-likelihood for each individual. Produces a histogram
#'         of log-likelihood values.
#' @importFrom atease @
#' @importFrom gbasics snpgeno
#' @importFrom mvbutils cq extract.named
#' @importFrom gbasics sqr
#' @import vecless
#' @keywords misc
#' @examples
#' # ll <- lglk_loci( snpg)
#' # hist( ll$sdiff, nc=30, col='grey')
#' # abline( v=0, col='green')
#' @export

"ilglk_geno" <- function(snpg, indiv_lglk_hist_pars=list(), quick=TRUE, showPlot = TRUE) {
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
      main   = 'Geno lglk by FISH', #sprintf( 'Geno lglk by FISH: multhresh=%5.2f', method, multhresh_indiv_lglk_fish),
      xlim   = range( ilglk),
      col    = "grey",
      border = NA,
      xlab   = '',
      nclass = 50)
  if( showPlot) {
      lv <- do.call( 'hist', c( list( x=ilglk), indiv_lglk_hist_pars))

      # some dens_SPA can fail -- do lapply and then weed out baddies
      mids_SPA <- lapply(lv$mids, function(x) try(dens_SPA( x), silent=TRUE))
      good_ind <- unlist(lapply(mids_SPA, class) != "try-error")
      mids_SPA <- unlist(mids_SPA[good_ind])
      # plot predicted density. Slowish with vecless 1.0
      lines( lv$mids[good_ind], diff( lv$breaks)[good_ind] * mids_SPA * n_samps,
            col='blue')
  }

return( ilglk)
}

#' lglk_loci(): Check individual genotypes for aggregate typicality
#' @rdname ilglk_geno
#' @export

"lglk_loci" <-
function( snpg) {
# When loci are presumably misbehaving, in that ilglk_geno looks
# ... OK on a "trusted subset" but not on all--- this may help

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

  # elg <- rowSums( pgeno * lpgeno)
  # elg2 <- rowSums( pgeno * sqr( lpgeno))
  # Or:
  elg[ l] := lpgeno[ l, g] %[g]% pgeno[ l, g]
  elg2[ l] := sqr( lpgeno[ l, g]) %[g]% pgeno[ l, g]
  sdelg <- sqrt( elg2 - sqr( elg))  # per individual

  s4 <- as.integer( unclass(snpg4))
  dim( s4) <- dim( snpg4)
  # olg_fl[ f, l] := lpgeno[ l, s4[ f, l]]    # NYI in vecless 1.0

  GG <- 1:4
  # olg_fl[ f, l] := lpgeno[ l, g] %[g]% (s4[f,l]==GG[g])
  # olg[ l] := ( SUM_ %[f]% olg_fl[ f, l]) / {n_f}
  n_f <- nrow( snpg4)
  olg[ l] := ( SUM_ %[f]% ( lpgeno[ l, g] %[g]% (s4[f,l]==GG[g]))) / {n_f}


return( returnList( elg, sdelg, olg, sdiff=sqrt( n_f) * (olg-elg) / sdelg))
} ## ported from MVB_kinference 25/9/18


#' map4todummy6(): Change SNP-genotype encoding
#' 
#' In case you have just a 4way-genotyped dataset but want to use
#' \code{\link{find_HSPs}} etc which require 6way, this will do it... All
#' \code{AAO} is mapped to \code{AA} and all \code{BBO} to \code{BB}. NA-values
#' should be unaffected (not tested). The \code{.$locinfo$use6} field will be
#' set to FALSE for all loci.
#' 
#' Warnings/errors should be issued if the inputs are wrong.
#' 
#' @param sg a \code{snpgeno} class object with \code{diplos} set to
#'           \code{genotypes4_ambig} (or \code{genotypes6}, in which case it
#'           will be returned unchanged).
#' @param quietly if FALSE and \code{sg} is already 6way-format, a warning is
#'                issued.
#' @return Modified version of the dataset, with encoding of \code{genotypes6}
#'         (see the code of \code{define_genotypes}). The \code{locinfo}
#'         attribute is augmented with a dummy \code{snerr} matrix, and the
#'         \code{use6} column is set to \code{FALSE}. Allele freqs are not
#'         estimated, so you'll need to do that separately.
#' @keywords misc
#' @export

"map4todummy6" <-
function( sg, quietly=FALSE) {
  stopifnot( sg %is.a% 'snpgeno')

  define_genotypes()
  dips <- sg$diplos
  if( my.all.equal( dips, genotypes6)) {
    if( !quietly) {
      warn( "already 6way; nothing to do")
    }
return( sg)
  } else if( !my.all.equal( dips, genotypes4_ambig)) {
stop( "Only meant for 'genotypes4_ambig' data")
  }

  # Not space-efficient...
  x <- as.vector( unclass( sg))
  for( i in 1:4) {
      x[ x==as.raw( i)] <- as.raw( i+10)
  }
  # ... and should leave x==FF (NA) alone

  for( i in 1:4) {
    x[ x==as.raw( i+10)] <- as.raw( match( substring( dips[ i], 1, 2), genotypes6)) # ie AAO -> AA
  }

  li <- sg$locinfo
  li$snerr <- matrix( 0.5, nrow=nrow( li), ncol=4, dimnames=list( NULL, cq( AA2AO, AO2AA, BB2BO,  BO2BB)))
  li$useN <- 4
  sg$locinfo <- li
  sg$diplos <- genotypes6

  sg <- unclass( sg)
  sg[] <- x
  class( sg) <- 'snpgeno'
return( sg)
}  ## ported from MVB_kinference 25/9/18


#' map6to4(): Bare documentation
#'
#' This function has only the bare minimum of documentation necessary for roxygen to
#' parse it. We should probably add some proper documentation here.
#'
#' From the name, it seems like this is the reverse operation of map4todummy6() -
#' it allows you to convert 6-way genotypes (ABO, BBO, etc.) to 4-way genotypes
#' (AA, BB, AB, BA), for use in functions that require one or the other.
#'
#' @param g6p0 a param
#' @param g6p1 a param
#' @export

map6to4 <- function(g6p0, g6p1){

  define_genotypes()

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

  return(list(g4p0=g4p0, g4p1=g4p1))

}

#' map6to3(): Bare documentation
#'
#' This function has only the bare minimum of documentation necessary for roxygen to
#' parse it. We should probably add some proper documentation here.
#'
#' From the name, it seems like this is the reverse operation of map4todummy6() -
#' it allows you to convert 6-way genotypes (ABO, BBO, etc.) to 3-way genotypes
#' (AAO, AB, BBOO), for use in functions that require one or the other.
#'
#' @param g6p0 a param
#' @param g6p1 a param
#' @export

map6to3 <- function(g6p0, g6p1){

  define_genotypes()

  # For the 4-ways, must condense g6p's

  map6to3 <- matrix( 0, 6, 3, dimnames=list( genotypes6, genotypes3_ambig))
  # AB & OO are OK; AAO should receive both AA and AO; etc
  mm <- match( genotypes6, substring( genotypes3_ambig, 1, 2), 0) # the "AA" bit of "AAO"...
  yup <- cbind( which(mm>0), mm[ mm>0])
  map6to3[ yup] <- 1
  mm <- match( genotypes6, substring( genotypes3_ambig, 2, 3), 0) # ... and the "AO" bit
  yup <- cbind( which(mm>0), mm[ mm>0])
  map6to3[ yup] <- 1
  mm <- match( genotypes6, substring( genotypes3_ambig, 3, 4), 0) # ... and the "OO" bit
  yup <- cbind( which(mm>0), mm[ mm>0])
  map6to3[ yup] <- 1

  # Really want g4p0[l,i,j] := map6to4[i,k6] %[k6]% g6p0[l,k6,m6] %[m6]% map6to4[m6,j]
  # ... but vecless can't presently handle multi-stages

  A[l,i,k] := g6p0[l,i,j] %[j]% map6to3[j,k]
  g3p0[l,i,j] := map6to3[ k,i] %[k]% A[l,k,j]

  A[l,i,k] := g6p1[l,i,j] %[j]% map6to3[j,k]
  g3p1[l,i,j] := map6to3[ k,i] %[k]% A[l,k,j]

  return(list(g3p0=g3p0, g3p1=g3p1))

}


#' predict_HSP_util(): Bare documentation
#'
#' This function has only the bare minimum of documentation necessary for roxygen to
#' parse it. We should probably add some proper documentation here.
#'
#' @param pIBD0 a param
#' @param pIBD1 a param
#' @param want_LOD_table a param. Defaults to FALSE
#' @param k a param. Defaults to 0.5
#' @export
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
  # expected value of the PLOD for HSP (mean of distn)
  E.HSP[l] := LOD[l,i,j] %[i,j]% Phsp[l,i,j]
  # expected value of the PLOD for UP (mean of distn)
  E.UP[l] := LOD[l,i,j] %[i,j]% Pup[l,i,j]
  E2.UP[l] := (LOD*LOD)[l,i,j] %[i,j]% Pup[l,i,j]
  V.UP <- E2.UP - sqr( E.UP)
    Ediff <- E.HSP - E.UP

    ## Code for mean and variance of LOD difference given coinheritance ---
    ## results needed for var.PLOD.kin
    P1 <- 2 * Phsp - Pup # of genopair, given exactly 1 coinherited allele
    e1[ l]:= LOD[l,i,j] %[i,j]% P1[l,i,j]
    e2.1[ l]:= (LOD*LOD)[l,i,j] %[i,j]% P1[l,i,j]
    v1 <- e2.1 - sqr( e1)
    e0 <- E.UP
    v0 <- V.UP
    matto <- cbind( e0, v0, e1, v1)
    
  # Standardized difference ie locus power: not so useful post hoc,
  #  but possibly interesting for 6 vs 4 comps
  sdiff <- (E.HSP - E.UP) / sqrt( V.UP)

#  Ediff <- unclass( Ediff)
#  V.UP <- unclass( V.UP)
#  sdiff <- unclass( sdiff)

  retval <- data.frame( Ediff, V.UP, sdiff, matto) ## matto comes in as 4 named columns, not as matto
  if( want_LOD_table) {
    retval@LOD <- gpLOD
    retval@PUP <- gpPUP
  }
return( retval)
}


#' @rdname hsp_power
#' @export
#' @importFrom mvbutils cq %except% %not.in%
#' @importFrom atease @ @<-
#' @importFrom gbasics make_genopairer sqr

"prepare_PLOD_SPA" <- function( geno6, n_pts_SPA_renorm=201) {
    ## To be run after hsp_power( ..., want_LOD_table=TRUE)
    
# n_pts_SPA_renorm should really be as big as R can handle without running
# out of memory but 201 should be OK I guess. If 201 and 301 give
# almost-identical results then all well!

stopifnot( all( cq( LOD4, LOD6, useN) %in% names( geno6@locinfo)))

  og <- options( vecless.print=FALSE)
  on.exit( options( og))

  # Combine 4way and 6way stuff into overall LOD and PUP
  extract.named( geno6@locinfo[ cq( useN, LOD6, LOD4, PUP6, PUP4, LOD3, PUP3)])
##  use4 <- !use6

  # DO NOT change actual genos though; they will be changed on-the-fly prior to kin-finding
  # ... code WOULD be this:
  # temp_snpg <- snpg
  # recode4to6temp <- function( x) { x[ x=='AO'] <- AA; x[ x=='BO'] <- BB; x}
  # temp_snpg[ , use4] <- recode4to6temp( snpg[, use4]) # (AA,AO) -> AA; (BB,BO) -> BB

  LOD <- LOD6
  PUP <- PUP6
  cn6 <- colnames( LOD6)
  cn4 <- colnames( LOD4)
  cn3 <- colnames( LOD3)

    ## Change only the entries with "full homz" since "single nulls" won't be accessed
    ## 6-to-4 equivalences
    for( ichangio in grep( 'AA|BB', cn6) %except% grep( 'AO|BO', cn6)){
        was1 <- substring( cn6[ ichangio], 1, 2)
        was2 <- substring( cn6[ ichangio], 4, 5)
        now1 <- sub( '(A|B)\\1', '\\1\\1O', was1)
        now2 <- sub( '(A|B)\\1', '\\1\\1O', was2)
        iget4 <- paste( now1, now2, sep='/')
        if( iget4 %not.in% cn4) { # reverse the order
            iget4 <- paste( now2, now1, sep='/')
        }
        LOD[ useN == 4, ichangio] <- LOD4[ useN == 4, iget4]
        PUP[ useN == 4, ichangio] <- PUP4[ useN == 4, iget4]
    }
    ## 6-to-3 equivalences
    for( ichangio in grep( 'AA|BB|OO', cn6) %except% grep( 'AO|BO', cn6)) {
        was1 <- substring( cn6[ ichangio], 1, 2)
        was2 <- substring( cn6[ ichangio], 4, 5)
        if( was1 == 'OO' | was2 == 'OO') {  ## crop out the double-nulls; par-format them
            if( was1 == 'OO') { now1 <- 'BBO'}
            if( was2 == 'OO') { now2 <- 'BBO'}
        } else { ## handle the AO | BO set
            now1 <- sub( '(A|B)\\1', '\\1\\1O', was1)
            now2 <- sub( '(A|B)\\1', '\\1\\1O', was2)
        } ## add the trailing 'O' for OO | BO | BB, giving BBOO
        if( now1 == 'BBO') { now1 <- 'BBOO' }
        if( now2 == 'BBO') { now2 <- 'BBOO' }
        iget3 <- paste( now1, now2, sep='/')
        if( iget3 %not.in% cn3) { # reverse the order
            iget3 <- paste( now2, now1, sep='/')
        }
        LOD[ useN == 3, ichangio] <- LOD3[ useN == 3, iget3]
        PUP[ useN == 3, ichangio] <- PUP3[ useN == 3, iget3]
    }
    
    ## For safety's sake, LOD( XO,...) := NA; should never be accessed

    hasO_4way <- grep( '(A|B)O', cn6)
    hasO_3way <- grep( '.O', cn6)
        
    LOD[ useN == 4, hasO_4way] <- NA # security in case of wrong access later for real data
    PUP[ useN == 4, hasO_4way] <- 0  
    LOD[ useN == 3, hasO_3way] <- NA # security in case of wrong access later for real data
    PUP[ useN == 3, hasO_3way] <- 0  

    LOD@mg <- make_genopairer( geno6@diplos)

  make_K <- function( PUP, LOD) { # ... while the sun skines

      # vecless **should** work just exporting := BUT doesn't seem to
      e <- new.env( parent=asNamespace( 'vecless'))
      # add sqr to the environment so that vecless can see it...
#      e$sqr <- gbasics::sqr
      e$sqr <- function( x) x*x # gbasics::sqr but it's not in gbasics!
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
#            rowSums( (SLL/S-gbasics::sqr( SL/S)))
            rowSums( (SLL/S-sqr( SL/S)))
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

    class(geno6) <- c("SPAgeno", "snpgeno")  ## new SPAgeno class to auto-update SPA
                                             ## after subset operations.

  # PUP and LOD are in Kenv now, so don't duplicate them in locinfo
  # They are really the "workhorse" versions, and are a bit cheaty, so don't want
  # them too public

return( geno6)
}


#' set_thresholds(): Bare documentation
#'
#' This function has only the bare minimum of documentation necessary for roxygen to
#' parse it. We should probably add some proper documentation here. This
#' function causes a NOTE from R CMD check, pertaining to the undefined global variable
#' 'keep_n'. set_thresholds() operates in mvbutils::mlocal({}) space, and every caller
#' function has 'keep_n' correctly defined.
#'
#' @param keeping a param
#' @param nlocal a param. Defaults to sys.parent().
#' @export
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
          keep_n / (length( subset1) * length( subset2) / (1+symmo)))
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

# Calculate some log-odds

#' calculate_IBD(): Bare documentation
#'
#' This function has only the bare minimum of documentation necessary for roxygen to
#' parse it. We should probably add some proper documentation here.
#'
#' @param lociar a param
#' @export

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


#' simtest_Kstuff(): check for (big) mistakes in SPA calcs
#'
#' This function has only the bare minimum of documentation necessary for roxygen to
#' parse it. We should probably add some proper documentation here.
#'
#' @param ck a 'SPAgeno' object
#' @param n a param. ?The number of animals to simulate?
#' @param nq a param. Defaults to 20. ?The number of quantiles?
#' @export
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

#' sqr(): multiply an input by itself.
#'
#' Does what it says on the tin. Might also be imported from MVButils,
#' depending on the caller.
#'
#' @param x an argument for which multiplication is possible.
#' @export

"sqr" <-
function( x) x*x  ## ported from MVB, 25/9/18


#' Predict variance of PLOD for HCPs and HTPs
#'
#' Aim is to work out how much your putative half-sibling pairs
#' (HSPs) might be contaminated by half-cousin pairs (HCPs) or
#' half-thiatic pairs (HTPs). HSP-selection is presumably based
#' on the pairwise PLODs for HSP:UP, taking all pairs where that
#' PLOD exceeds some threshold. Given the allele freqs, the mean
#' PLOD is predictable when the truth is UP, HCP, HTP, or HSP.
#' The variance is only predictable for UPs, because linkage makes
#' loci non-independent for kin. However, an empirical variance
#' can be estimated for HSPs based on the observed PLODs above
#' some super-high threshold, e.g., the mean PLOD when truth=HSP
#' (this preliminary variance step could be iterated using
#' different super-high thresholds). Based on the empirical
#' variance for HSPs and the analytical variance for UPs, we
#' basically know how much linkage there might be, so we can
#' predict the PLOD variances for the other kin-pair types. The
#' wrinkle is that it also depends somewhat on the finer-scale
#' organization of the genome, i.e. whether it's lots of
#' chromosomes with no crossover, or a few chromosomes with lots of
#' crossover. \code{var.PLOD.kin} therefore calculates two versions,
#' one assuming the genome is entirely made up of equal-sized
#' chromosome with zero crossover, and the other assuming the
#' genome is a single chromosomes with crossover according to a
#' random process.
#'
#' It's assumed that lots of loci are being used, so that the mix
#' of loci on each "chromo", or the splatter of loci along the
#' single "megachromo", always matches the overall population, on
#' law-of-large-numbers grounds.
#'
#' Stuff like uncertainly in allele frequencies, and in the PLOD
#' variance for HSPs, needs to be accounted for externally, by
#' repeatedly drawing from the posteriors and re-calculating the
#' PLODs and re-running this function.
#'
#' If the variance estimates show really good separation between
#' the kin-pair types, then one could refine the "preliminary
#' variance" step by reducing the super-high threshold (and assuming
#' a truncated-Normal distribution). This might be worthwhile if the
#' preliminary variance otherwise has to be based on a very small
#' number of no-brainer HSPs. The "logical conclusion" of That Kind
#' Of Thing is some kind of MLE involving estimating the population
#' of different types of kin, and we really don't want to go there
#' for now (since that should include the population dynamics
#' shebang). In other words, we'd end up linking the genetic
#' kin-finding model to the population dynamics model, which makes
#' life statistically harder. And god knows it's hard enough.
#' Anyway, if we were taking that approach, it might well be better
#' to avoid PLODs altogether and instead go for inferences about
#' the actual ppn of co-inherited loci, from which
#' estimates-of-co-inherited-variance and inferences about kin-ppns
#' can be made.
#'
#' @param linfo locus info, a \code{data.frame} with columns
#'              \code{e0}, \code{e1}, \code{v0}, \code{v1},
#'              \code{count}. Each row is one "type" of locus, i.e.,
#'              with roughly the same values of e/v 0/1, and
#'              \code{count} says how many such loci there are. e/v
#'              0/1 are means and variances of the per-locus LOD (note
#'              no P) when the locus is or isn't co-inherited.
#' @param emp.V.HSP empirical variance of PLOD for deffo HSPs. You're
#'                  supposed to be running this on real data, so that
#'                  \code{emp.V.HSP} is an actual number; however, for
#'                  testing purposes, you can set up an artificial
#'                  version via \code{C.equiv} below.
#' @param kin.true Var[PLOD|kin.true]. The default, \code{"HCP"},
#'                 corresponds to 4 meioses, \code{"HTP"} to 3, and
#'                 \code{"HSP"} to 2; the latter is only for debugging,
#'                 since it should reproduce the original empirical
#'                 variance!
#' @param debug Logical flag. Defaults to FALSE.
#' @param C.equiv for artificial test, with \code{emp.V.HSP} set to the
#'                no-crossover variance from \code{C.equiv} chromos
#'                (need not be integer). Ignored if \code{emp.V.HSP} is
#'                set.
#' @return Vector with names \code{V0}, \code{Vx}, \code{C.hat},
#'         \code{rho.hat}, \code{n.meio.<XXX>}. First two are variances
#'         under the no- and all-crossover scenarios; \code{C.hat} is
#'         estimated equivalent number of chromosomes, and \code{rho.hat}
#'         is per-locus crossover rate, under the same scenarios.
#'         \code{n.meio.<XXX>} (where \code{"XXX"} is set to
#'         \code{kin.true}) shows how many meioses are involved; yes, a
#'         perverse way to return that piece of info.
#' @examples
#' \dontrun{
#' # H1CPs, for "simulation" equiv to 22 no-Xover chromos
#' var.PLOD.kin( data.frame( count=45, e0=-1, e1=2, v0=0.03, v1=0.02), C.equiv=22)
#' #        VUP       V.HSP          V0          Vx       C.hat     rho.hat n.meio.HCP
#' #     1.3500    208.2273     91.9010    101.8270     22.0000      0.2616      4.0000
#'
#' # HTPs
#' var.PLOD.kin( data.frame( count=45, e0=-1, e1=2, v0=0.03, v1=0.02), C.equiv=22, kin='HTP')
#' #       VUP      V.HSP         V0         Vx      C.hat    rho.hat n.meio.HTP
#' #    1.3500   208.2273   156.5642   186.1779    22.0000     0.2616     3.0000
#'
#' # Vars would go up with fewer equiv chromos
#' var.PLOD.kin( data.frame( count=45, e0=-1, e1=2, v0=0.03, v1=0.02), C.equiv=5, kin='HTP')
#' #       VUP      V.HSP         V0         Vx      C.hat    rho.hat n.meio.HTP
#' #   1.35000  912.37500  684.67500  786.47854    5.00000    0.04951    3.00000
#' }
#' @export

"var.PLOD.kin" <-
function( linfo, emp.V.HSP=V.noX( C.equiv, 2),
         kin.true=c( 'HCP', 'HTP', 'HSP'), debug=FALSE, C.equiv=NULL) {

  kin.true <- match.arg( kin.true)
  n.meio <- c( HCP=4, HTP=3, HSP=2)[ kin.true]

# linfo should be a DF with cols e0, e1, v0, v1, count
  names( linfo) <- names( linfo) %&% '.l'
  extract.named( linfo)
  L <- sum( count.l)
  pi <- count.l / L

  e0 <- pi %*% e0.l
  e1 <- pi %*% e1.l
  v0 <- pi %*% v0.l
  v1 <- pi %*% v1.l

  # noX = no crossover, i.e. all in separate equal chromos
  V.noX <- function( C, meioses) {
      p <- 2 ^ (1-meioses) # Marginal prob of coinheritance = 1/2 for HSPs, 1/8 for HCPs
      EofV <- L * v0 + L * p *(v1-v0)
      VofE <- sqr( L * (e1-e0)) * p * (1-p) / C
    return( EofV + VofE)
  }

  # MoM for HSPs:
  if( !is.null( C.equiv)) {
    C.hat <- C.equiv <- min( L, C.equiv)
  } else {
    C.hat <- if( V.noX( 1, 2) < emp.V.HSP) 1 else
        if( V.noX( L, 2) > emp.V.HSP) L else
        find.root( V.noX, target=emp.V.HSP, start=1, step=1, fdirection='decreasing',
            min.x=1, max.x=L, meioses=2)
  }
  V0 <- V.noX( C.hat, meioses=n.meio)

  # Entirely Xover on one single chromo
  # Should these be done separately by locus type, then averaged??
  ell <- 1 %upto% (L-1)
  e2.1 <- pi %*% (sqr( e1.l) + v1.l)
  e2.0 <- pi %*% (sqr( e0.l) + v0.l)
  e1.0 <- e0 # pi %*% e0.l
  e1.1 <- e1

  phi <- function( r, mul, s) 0.5 * (1 + s*exp( -mul*r*ell))
  PSstar <- P0star <- array( 0, c( L-1, 2, 2, 2))
  si2 <- slice.index( PSstar, 2)
  si3 <- slice.index( PSstar, 3)
  si4 <- slice.index( PSstar, 4)

  V.allX <- function( rho, meioses) {
    # PSstar[ L,i,j,k] is Pr[ composite state @ L = 1 | starting states @ 0 are i, j, k]
    # j & k series irrel for HSP
    # slice.index() == 1 or 2
    # so 3 - 2*s.i. == -1 or +1
    # and 2-s.i. == 1 or 0

    # This code could be much more efficient; pre-compute P0star, and much reduce the phi-calcs
    PSstar[] <<- phi( rho, 4, 3 - 2*si2) *
        (if( meioses>2) phi( rho, 2, 3 - 2*si3) else 1) *
        (if( meioses>3) phi( rho, 2, 3 - 2*si4) else 1)

    P0star[] <<-  (2-si2) *
        (if( meioses>2) 2-si3 else 1) *
        (if( meioses>3) 2-si4 else 1)

    p11.8 <- apply( PSstar * P0star, 1, sum)
    p01.8 <- apply( PSstar * (1-P0star), 1, sum)
    p00.8 <- apply( (1-PSstar) * (1-P0star), 1, sum)

    p11 <- p11.8 / 8
    p01.10 <- p01.8 / 4 # both ways
    p00 <- p00.8 / 8

    p <- 2 ^ (1-meioses) # Marginal prob of coinheritance = 1/2 for HSPs, 1/8 for HCPs

    EL <- L * ( e1.1*p + e1.0*(1-p))
    EL2 <- L * ( e2.1*p + e2.0*(1-p)) +
        2 * (L-ell) %*% ( sqr( e1.1)*p11 + e1.1*e1.0*p01.10 +  sqr( e1.0)*p00)
    VL <- EL2 - sqr( EL)
  }

  if( debug) {
    mtrace( V.allX) # surely wanna
    0 # help debugging...
  }

  # If no Odis (e.g. with very few loci!) then no point in going further
  rho.hat <- if( C.hat > L-1) 100 else
      find.root( V.allX, target=emp.V.HSP, start=1/L, step=0.2/L,
          fdirection='decreasing', min.x=0, meioses=2)
  Vx <- V.allX( rho.hat, meioses=n.meio)

return( unlist( returnList( VUP=L*v0, V.HSP=emp.V.HSP, V0, Vx, C.hat, rho.hat, n.meio)))
}

## first two lines here are for Rcpp to work
#' add_list_defaults(): Bare documentation
#'
#' This function has only the bare minimum of documentation necessary for roxygen to
#' parse it. We should probably add some proper documentation here.
#'
#' @param l a param
#' @param ... other params, passed to something.
#' @export
#' @useDynLib kinference
#' @importFrom Rcpp evalCpp
#' @importFrom mvbutils %without.name% ?

"add_list_defaults" <-
    function( l, ...) {
###### Add defaults to list 'l' if not already in 'l'
  defaults <- list(...)
  l <- c( l, defaults %without.name% names( l))
return( l)
}


#' define_genotypes(): Documentation (warning?) ported from gbasics
#'
#' You probably don't want to run this... it's documented only so it
#' can be exported for other packages to use. It lets you refer to a
#' genotype "AB" as just AB etc. It creates them *in the caller*, so
#' a lot of routines run it internally as a first step so that their
#' subsequent code is simpler- and that is the intended use. If you
#' run it from the R prompt, the things will be created in
#' ‘.GlobalEnv’.
#'
#' The full set of possible genotypes are created as objects; e.g.
#' there will be an object ‘AB’ which is really a string "AB", but
#' with a ‘noquote’ attribute so that you just see it without the
#' quotes.
#'
#' It also creates some allowed groupings of genotypes (e.g. if there
#' is no C allele at a locus; whether or not distinguishing single
#' nulls, etc).
#'
#' This particular version sets up 3-allele-plus-null cases (A,B,C,O)
#' and various combinations.
#'
#' Since you are not really meant to call this, I'm not going to
#' improve the documentation.
#'
#' @param nlocal frame number, or environment, to create things in -
#'               see 'mlocal'. Leave this alone unless you _really_
#'               know what you're doing. Defaults to sys.parent()
#' @return None, but various objects are created; see the code.
#' @export
#' @importFrom mvbutils cq extract.named named mlocal
#' @seealso 'mlocal' in package 'mvbutils'.
#' @examples
#' ## Not run:
#'
#' # define_genotypes()
#' # A # thar she
#' ## End(Not run)

"define_genotypes" <-
    function( nlocal=sys.parent()) mlocal({
  ABCO <- named( cq( A, B, C, O))
  extract.named( ABCO) # A, B, C, and O

  genotypes <- cq( OO, AO, BO, AB, AA, BB, AAO, BBO, AC, BC, CO, CC, "CCO")
  genotypes_ambig <- cq( OO, AB, AC, BC, AAO, BBO, "CCO") ## sans quotes, CCO is treated as global
  genotypes4_ambig <- cq( OO, AB, AAO, BBO)
  genotypes6 <- cq( AA, AB, AO, BB, BO, OO)
  genotypes_C <- cq( AA, AB, AC, AO, BB, BC, BO, CC, CO, OO)
  genotypes3_ambig <- cq( AB, AAO, "BBOO") ## sans quotes, BBOO is also treated as global

  NA_geno <- as.raw(255)
  
  for( ig in genotypes) {
    # assign( ig, structure( as.raw( match( ig, genotypes)), class='ABOSNP'))
    assign( ig, structure( ig, class='noquote')) # for nicer printing
  }
})



#' inv_CDF_SPA2(): Bare documentation
#'
#' This function has only the bare minimum of documentation necessary for roxygen to
#' parse it. We should probably add some proper documentation here.
#'
#' @param p a param
#' @param K a param
#' @param dK a param.
#' @param ddK a param.
#' @param tol a param. Defaults to formals( ridder)$tol
#' @importFrom stats pnorm dnorm qnorm
#' @importFrom gbasics logit inv.logit
#' @export
"inv_CDF_SPA2" <-
    function( p, K, dK, ddK, tol=formals( ridder)$tol) {

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


#' calculate_LOD_HSP(): Bare documentation
#'
#' This function has only the bare minimum of documentation necessary for roxygen to
#' parse it. We should probably add some proper documentation here.
#'
#' @param lociar a param
#' @param k a param. Defaults to 0.5
#' @export

"calculate_LOD_HSP" <-
    function(lociar, k=0.5) {

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


#' PLOD_loghisto(): PLOD distro log-frequency plot with UP and HSP expetations
#'
#' Plots a log-frequency histogram for the output of \code{find_HSPs()}, with
#' the expected mean PLOD for unrelated pairs in red, the expected distribution
#' of unrelated pairs in blue, and the expected mean PLOD for HSPs in green.
#'
#' @param hsps the output of a call to \code{find_HSPs()}
#' @param UP plot the mean PLOD for unrelated pairs? Default TRUE
#' @param HSP plot the mean PLOD for HSPs? Default TRUE
#' @param showUP plot the expected density curve for unrelated pairs using the SPA
#'               approximation (default TRUE), Normal approximation (default FALSE),
#'               both, or neither. SPA approximation will plot in blue (or colour 4),
#'               Normal in magenta (or colour 6).
#' @param ... additional pars, passed to \code{plot}
#' @export

"PLOD_loghisto" <-
    function(hsps, UP = TRUE, HSP = TRUE, showUP = c(SPA = TRUE, Normal = FALSE), ...) {

        plot( hsps@bins,log(hsps@n_PLODs_in_bin),type='S', ...,  xlab="PLOD",
             ylab="log(Frequency)")
        if( UP) { abline(v = hsps@mean_theory, col = 2, lwd = 2) }
        if( HSP) { abline(v = hsps@mean_HSP, col = 3, lwd = 2) }
        if( showUP["SPA"]) {
            lines(hsps@bins,log(diff(c(0,hsps@binprobs))*sum(hsps@n_PLODs_in_bin[hsps@bins<0])),
                  lwd=2,col=4)
        }
        if( showUP["Normal"]) {
            lines(hspsa@bins,log(diff(c(0,pnorm(hspsa@bins, mean = hspsa@mean_theory,
                                      sd = sqrt(hspsa@var_theory)) *
                                      sum(hspsa@n_PLODs_in_bin[hspsa@bins<0])))),
                  lwd = 2, col = 6) ## Normal approx
        }
    }


#' HSP_histo(): PLOD distro frequency plot for the HSP and/or HTP bump regions
#'
#' Plots an absolute-frequency histogram for the output of \code{find_HSPs()}, with
#' the lower bound set by the user. Lower bounds should be set to exclude (as much
#' as possible) the UP bump, as this will otherwise swamp the signal from the HSP
#' bump. Users must manually set a lower bound for full-sibling PLODs ('fullsib_cut')
#' on order to exclude full-siblings from the variance estimate for HSP PLODs.
#'
#' @param hsps the output of a call to \code{find_HSPs()}
#' @param lb PLOD lower bound for plot extent. Should exclude the UP bump
#' @param ub PLOD upper bound for plot extent. Defaults to maximum PLOD score plus a
#'           little padding
#' @param fullsib_cut PLOD score above which there are only full-sibs
#' @param bin hist bin width. Default 5. lb, ub, and bin together define 'breaks',
#'            so you can't pass 'breaks' via \code{...}
#' @param HSPmean plot the mean PLOD for HSPs? Default TRUE
#' @param HSPdist plot the distribution of PLOD for HSPs? Default TRUE
#' @param ... additional pars, passed to \code{hist()}
#' @export

"HSP_histo" <-
    function(hsps, lb, ub = max(hsps$PLOD)+10, fullsib_cut, bin = 5,
             HSPmean = TRUE, HSPdist = TRUE, ...) {
        
        hist.plod=hist(hsps$PLOD[hsps$PLOD > lb & hsps$PLOD < ub],breaks=seq(lb, ub, bin),
                       col=8,main="HSP PLOD",xlab="PLOD", ...)
        if( HSPmean) {
            E.hsp = hsps@mean_HSP
            abline(v=E.hsp,lwd=2,col=2)
        }
        if( HSPdist) {
            V.hsp=mean(sqr(hsps$PLOD[hsps$PLOD>E.hsp & hsps$PLOD < fullsib_cut]-E.hsp))
            obs.num <- hist.plod$counts
            exp.num <- 2*sum(hsps$PLOD>E.hsp & hsps$PLOD<fullsib_cut)*
                (pnorm(hist.plod$breaks[-1],E.hsp,sqrt(V.hsp))-
                 pnorm(hist.plod$breaks[-length(hist.plod$breaks)],E.hsp,sqrt(V.hsp)))
            points(hist.plod$mids,exp.num,pch=16,col=4,type='b')
        }
    }


#' HSP_oddness_oneway(): show the distributions of oddness metrics over an HSP histo
#'
#' Plots an absolute-frequency histogram for the output of \code{find_HSPs()}, with
#' additional coloured regions showing individuals with odd-looking CLOD scores
#' (see \code{check_FPosity()} ), ilglk stat (see \code{ilglk_geno()} ), or hetz
#' stat (see \code{hetzminoo_fancy()} ). The user should choose a PLOD range for plotting
#' that excludes the UP bump, as unrelated pairs can swamp the signal from HSPs and HTPs.
#' \code{HSP_oddness_oneway()} shows cases where at least one member of the pair has an
#' unusual score in the oddness metric, whereas \code{HSP_oddness_twoway()} shows cases
#' where both members of the pair have an unusual score in the oddness metric.
#'
#' @param hsps the output of a call to \code{find_HSPs()}
#' @param snpg the 'snpgeno' or 'SPAgeno' object from which 'hsps' was built
#' @param lb PLOD lower bound for plot extent. Should exclude the UP bump
#' @param ub PLOD upper bound for plot extent. Defaults to the maximum PLOD score plus a
#'           little padding
#' @param bin hist bin width. Default 5. lb, ub, and bin together define 'breaks', so you
#'            can't pass 'breaks' via \code{...}
#' @param CLOD_prop the quantile of CLOD below which animals are highlighted. Default 0.001
#' @param ilglk_prop the quantile of ilglk stat below which animals are highlighted. Default 0.001
#' @param hetz_prop the quantile of hetz stat below which animals are highlighted. Default 0.001
#' @param ... additional pars, passed to \code{hist()}
#' @seealso PLOD_oddness_oneway
#' @export

"HSP_oddness_oneway" <-
    function(hsps, snpg, lb, ub = max(hsps$PLOD)+10, bin = 5,
             CLOD_prop = 0.001, ilglk_prop = 0.001, hetz_prop = 0.001,
             ...) {

        cloddo <- check_FPosity(snpg)
        clod.stat <- log(pnorm( 5, mean=cloddo$ECLOD, sd=sqrt( cloddo$VCLOD), lower.tail=FALSE))
        hetz.poor<- hetzminoo_fancy(snpg, 'poor', showPlot = FALSE)
        ilglk<- ilglk_geno(snpg, showPlot = FALSE)

        hist1=hist(hsps$PLOD[ hsps$PLOD > lb], breaks = seq( lb, ub, bin),
                   main="HSP PLOD",xlab="PLOD", ...)
        for(i in 1:3) {
            if(i==1) {
                tf<- clod.stat[hsps$i] < quantile(clod.stat, CLOD_prop) |
                    clod.stat[hsps$j] < quantile(clod.stat, CLOD_prop)
                colour <- rgb(1,0,0,alpha = 0.5)
            }
            if(i==2) {
                tf<- ilglk[hsps$i]< quantile(ilglk, ilglk_prop) |
                    ilglk[hsps$j] < quantile(ilglk, ilglk_prop)
                colour = rgb(0,1,0, alpha = 0.5)
            }
            if(i==3) {
                tf<- hetz.poor[hsps$i] < quantile(hetz.poor, hetz_prop) |
                    hetz.poor[hsps$j]< quantile(hetz.poor, hetz_prop)
                colour = rgb(0,0,1, alpha = 0.5)
            }
            hist2=hist( hsps$PLOD[ hsps$PLOD > lb & tf],breaks=seq( lb,ub,bin), add=T,
                       border = colour, col = colour)
        }
        legend('topright',c('low CLOD stat','low ilglk stat','low hetz stat'),
               title='ONE member of the pair has:',
               fill = c(rgb(1,0,0, 0.5), rgb(0,1,0, 0.5), rgb(0,0,1, 0.5)))
    }


#' @rdname HSP_oddness_oneway
#' @export

"HSP_oddness_twoway" <-
    function(hsps, snpg, lb, ub = max(hsps$PLOD)+10, bin = 5,
             CLOD_prop = 0.001, ilglk_prop = 0.001, hetz_prop = 0.001,
             ...) {
        
        cloddo <- check_FPosity(snpg)
        clod.stat <- log(pnorm( 5, mean=cloddo$ECLOD, sd=sqrt( cloddo$VCLOD), lower.tail=FALSE))
        hetz.poor<- hetzminoo_fancy(snpg, 'poor', showPlot = FALSE)
        ilglk<- ilglk_geno(snpg, showPlot = FALSE)

        hist1=hist(hsps$PLOD[ hsps$PLOD > lb], breaks = seq( lb, ub, bin),
                   main="HSP PLOD",xlab="PLOD", ...)
        for(i in 1:3) {
            if(i==1) {
                tf<- clod.stat[hsps$i] < quantile(clod.stat, CLOD_prop) &
                    clod.stat[hsps$j] < quantile(clod.stat, CLOD_prop)
                colour <- rgb(1,0,0,alpha = 0.5)
            }
            if(i==2) {
                tf<- ilglk[hsps$i]< quantile(ilglk, ilglk_prop) &
                    ilglk[hsps$j] < quantile(ilglk, ilglk_prop)
                colour = rgb(0,1,0, alpha = 0.5)
            }
            if(i==3) {
                tf<- hetz.poor[hsps$i] < quantile(hetz.poor, hetz_prop) &
                    hetz.poor[hsps$j]< quantile(hetz.poor, hetz_prop)
                colour = rgb(0,0,1, alpha = 0.5)
            }
            hist2=hist( hsps$PLOD[ hsps$PLOD > lb & tf],breaks=seq( lb,ub,bin), add=T,
                       border = colour, col = colour)
        }
        legend('topright',c('low CLOD stat','low ilglk stat','low hetz stat'),
               title='BOTH members of the pair have:',
               fill = c(rgb(1,0,0, 0.5), rgb(0,1,0, 0.5), rgb(0,0,1, 0.5)))
    }


#' PLOD_oddness_oneway(): show densities of oddness metrics over the PLOD range
#'
#' Plots the percentage of all pairs in each bin with an 'unusually' low
#' CLOD score, ilglk stat, or hetz stat, across the range of PLOD. \code{PLOD_oddness_oneway()}
#' shows the percentage of cases where one member has a low score, and
#' \code{PLOD_oddness_twoway()} shows the percentage of cases where both members of the pair
#' have a low score.
#' 
#' @param hsps the output of a call to \code{find_HSPs()}
#' @param snpg the 'snpgeno' or 'SPAgeno' object from which 'hsps' was built
#' @param lb PLOD lower bound for plot extent. Should exclude the UP bump
#' @param ub PLOD upper bound for plot extent. Defaults to the maximum PLOD score plus a
#'           little padding
#' @param bin hist bin width. Default 5
#' @param CLOD_prop the quantile of CLOD below which animals are highlighted. Default 0.001
#' @param ilglk_prop the quantile of ilglk stat below which animals are highlighted. Default 0.001
#' @param hetz_prop the quantile of hetz stat below which animals are highlighted. Default 0.001
#' @param ... additional pars, passed to \code{plot()}. 'ylim' and 'breaks' are set internally,
#'            so you cannot pass them via \code{...}.
#' @seealso HSP_oddness_oneway
#' @export

"PLOD_oddness_oneway" <-
    function(hsps, snpg, lb = min(hsps$PLOD)-10, ub = max(hsps$PLOD)+10, bin = 5,
             CLOD_prop = 0.001, ilglk_prop = 0.001, hetz_prop = 0.001, ...) {

        cloddo <- check_FPosity(snpg)
        clod.stat <- log(pnorm( 5, mean=cloddo$ECLOD, sd=sqrt( cloddo$VCLOD), lower.tail=FALSE))
        hetz.poor<- hetzminoo_fancy(snpg, 'poor', showPlot = FALSE)
        ilglk<- ilglk_geno(snpg, showPlot = FALSE)

        hist1 <- hist( hsps$PLOD[hsps$PLOD>lb], breaks=seq( lb,ub,bin), plot=F)
        
        tfA<- clod.stat[hsps$i] < quantile(clod.stat, CLOD_prop) |
            clod.stat[hsps$j] < quantile(clod.stat, CLOD_prop)
        histA <- hist(hsps$PLOD[hsps$PLOD>lb & tfA], breaks=seq(lb,ub,bin), plot=F)
        
        tfB<- ilglk[hsps$i]< quantile(ilglk, ilglk_prop) |
            ilglk[hsps$j] < quantile(ilglk, ilglk_prop)
        histB <- hist(hsps$PLOD[hsps$PLOD>lb & tfB], breaks=seq(lb,ub,bin), plot=F)
        
        tfC<- hetz.poor[hsps$i] < quantile(hetz.poor, hetz_prop) |
            hetz.poor[hsps$j]< quantile(hetz.poor, hetz_prop)
        histC <- hist(hsps$PLOD[hsps$PLOD>lb & tfC], breaks=seq(lb,ub,bin), plot=F)

        yup <- max(na.omit(c(histA$counts/hist1$counts, histB$counts/hist1$counts,
                             histC$counts/hist1$counts)))
        plot(hist1$breaks[-1], histA$counts/hist1$counts, type='b', col=(2),
             pch=(15), xlab="PLOD value", ylab="Percent of fish pairs",
             ylim = c(0, yup), ...)
        points(hist1$breaks[-1], histB$counts/hist1$counts, type='b',col=(3), pch=(16))
        points(hist1$breaks[-1], histC$counts/hist1$counts, type='b',col=(4), pch=(17))
        legend('topright',c('low CLOD stat','low ilglk stat','low hetz stat'),
               title='ONE member of the pair has:', pch=15:17,col=2:4)
    }



#' @rdname PLOD_oddness_oneway
#' @export

"PLOD_oddness_twoway" <-
    function(hsps, snpg, lb = min(hsps$PLOD)-10, ub = max(hsps$PLOD)+10, bin = 5,
             CLOD_prop = 0.001, ilglk_prop = 0.001, hetz_prop = 0.001, ...) {

        cloddo <- check_FPosity(snpg)
        clod.stat <- log(pnorm( 5, mean=cloddo$ECLOD, sd=sqrt( cloddo$VCLOD), lower.tail=FALSE))
        hetz.poor<- hetzminoo_fancy(snpg, 'poor', showPlot = FALSE)
        ilglk<- ilglk_geno(snpg, showPlot = FALSE)

        hist1 <- hist( hsps$PLOD[hsps$PLOD>lb], breaks=seq( lb,ub,bin), plot=F)
        
        tfA<- clod.stat[hsps$i] < quantile(clod.stat, CLOD_prop) &
            clod.stat[hsps$j] < quantile(clod.stat, CLOD_prop)
        histA <- hist(hsps$PLOD[hsps$PLOD>lb & tfA], breaks=seq(lb,ub,bin), plot=F)
        
        tfB<- ilglk[hsps$i]< quantile(ilglk, ilglk_prop) &
            ilglk[hsps$j] < quantile(ilglk, ilglk_prop)
        histB <- hist(hsps$PLOD[hsps$PLOD>lb & tfB], breaks=seq(lb,ub,bin), plot=F)
        
        tfC<- hetz.poor[hsps$i] < quantile(hetz.poor, hetz_prop) &
            hetz.poor[hsps$j]< quantile(hetz.poor, hetz_prop)
        histC <- hist(hsps$PLOD[hsps$PLOD>lb & tfC], breaks=seq(lb,ub,bin), plot=F)

        yup <- max(na.omit(c(histA$counts/hist1$counts, histB$counts/hist1$counts,
                             histC$counts/hist1$counts)))
        plot(hist1$breaks[-1], histA$counts/hist1$counts, type='b', col=(2),
             pch=(15), xlab="PLOD value", ylab="Percent of fish pairs",
             ylim = c(0, yup), ...)
        points(hist1$breaks[-1], histB$counts/hist1$counts, type='b',col=(3), pch=(16))
        points(hist1$breaks[-1], histC$counts/hist1$counts, type='b',col=(4), pch=(17))
        legend('topright',c('low CLOD stat','low ilglk stat','low hetz stat'),
               title='BOTH members of the pair have:', pch=15:17,col=2:4)
    }

#' check6and4(): bare documentation
#'
#' Checks 6-way and 4-way genotype frequencies against HWE expectations, and generates
#' plots of observed / expected frequencies. Recently moved into kinference from
#' genocalldart.
#' @param geno6 a genotype object with 4-way and 6-way genocalls
#' @param thresh_pchisq_6and4 thresholds for 'bad' and 'really bad' p-values
#' @param return_what a character vector, defaulting to c('just_pvals', 'all')
#' @param extra_title a character string to be added to the bottom-right corner of all plots.
#'                    Best if < 25 characters.
#' @seealso geno_deambig_ABC
#' @importFrom grDevices rgb
#' @importFrom graphics legend
#' @importFrom graphics plot
#' @importFrom graphics points
#' @importFrom stats na.omit
#' @importFrom stats quantile
#' @export

"check6and4" <-
function( geno6, 
    thresh_pchisq_6and4, 
    return_what=c( 'just_pvals', 'all'),
    extra_title = "") {
##########    
  n_fish <- nrow( geno6)
  n_loci <- ncol( geno6)
  define_genotypes()
  diplos <- geno6@diplos
stopifnot( my.all.equal( sort( genotypes6), sort( diplos)))
  
  p6 <- with( geno6@locinfo, 
      calc_g6probs( pbonzer[,'A'], pbonzer[,'B'], pbonzer[,'C'], 
          snerr=snerr))[,genotypes6]
  # At some point, will also need p6_IBD... later...

  # Chi-sq check: requires gpred & gobs attr on geno6
  g6pred <- p6 * n_fish
  g6obs <- matrix( 0, n_loci, 6, dimnames=list( NULL, genotypes6))
  
  for( ig in genotypes6) {
    g6obs[,ig] <- colSums( geno6==match( ig, diplos))
  }

  pval6 <- chisq_genofreq_check( geno6, gobs=g6obs, gpred=g6pred, test='G', 
      thresh_pchisq_loci=thresh_pchisq_6and4, trim=FALSE, extra_title = extra_title)@locinfo$pval

  g4obs <- gtab6to4( g6obs)
  g4pred <- gtab6to4( g6pred)
  pval4 <- chisq_genofreq_check( geno6, gobs=g4obs, gpred=g4pred, test='G', 
      thresh_pchisq_loci=thresh_pchisq_6and4, trim=FALSE, extra_title = extra_title)@locinfo$pval
  
  return_what <- match.arg( return_what)
  if( return_what=='all') {
    geno6@locinfo$pval6 <- pval6
    geno6@locinfo$pval4 <- pval4
return( geno6)
  } else {
return( returnList( pval6, pval4))
  }
}


#' chisq_genofreq_check(): check genotype frequencies against HWE expectations
#'
#' Checks observed genotype frequencies against expected frequencies, presumably with expectation
#' defined by HWE.
#' @param lociar a snpgeno object
#' @param gpred predicted allele frequencies
#' @param gobs observed allele frequencies
#' @param thresh_pchisq_loci a param. Presumably, a threshold p-val for flagging loci with
#'                           suspicious-looking allele frequencies.
#' @param test a character string, either "Pearson" or "G"
#' @param trim TRUE or FALSE. TRUE will keep only above max thresh_pchisq_loci. Arguably better as
#'             to be done post-hoc.
#' @param seq_paxis numeric#'
#' @param extra_title a character string to be added to the bottom-right corner of all plots.
#'                    Best if < 25 characters.
#' @importFrom stats pchisq
#' @importFrom graphics mtext
#' @importFrom graphics par
#' @seealso kinference::check6and4
#' @export

"chisq_genofreq_check" <-
function( lociar, 
    gpred= lociar@gpred,
    gobs= lociar@gobs,
    thresh_pchisq_loci,  # NULL to not worry; 1 value a threshold; 2 vals to inspect "iffy" ones
    test,  # 'Pearson' or 'G'
    trim, # TRUE to keep only above max thresh_pchisq_loci. Arguably better done post hoc...
    seq_paxis=0.025,
    extra_title = "") {
##########    
# Either 6- or 4-geno version should work
# 1 DoF in either case
# Assumed null distro of chisq(1) is pretty approximate
  n_loci <- nrow( gobs)

  # After ML, should never happen that gobs>0 & gpred==0... but we'll check
stopifnot( !any( gobs>0 & gpred==0))  

  chistat <- if( test=='Pearson') 
      rowSums( sqr( gobs - gpred) / gpred) 
    else if( test=='G') # must handle 0log0 which is 0 but R doesn't know that (cf nlogp func in my Pascal armoury)
      2 * rowSums( gobs * log( ifelse( gobs>0, gobs/gpred, 1)))
    else 
stop( 'test must be "Pearson" or "G"')
      
  DoF <- ncol( gpred)-3 # ... maybe ... !?
  pval <- pchisq( chistat, df=DoF, lower.tail=FALSE) ### df = ???
  lociar@locinfo$pval <- pval
  liffies <- length( thresh_pchisq_loci)
  keep_loci <- if( liffies) pval > max( thresh_pchisq_loci) else rep( TRUE, n_loci)
  iffy_loci <- !keep_loci & (pval > min( thresh_pchisq_loci))
  
  # Plot histogram of pvals - should be approximately uniformly distributed if the loci are behaving as we would like
  # [ should follow recordo paradigm as per geno_deambig ]
  par( mfrow=c(1,1)) # just one plot on first page
  hist(pval, 
      main=sprintf( "%s GoF of %i-genotypes: pval from chisq( %i)", test, ncol( gpred), DoF), 
      xlab="P-value: LOW == BAD", 
      breaks=seq(0,1,seq_paxis),
      xlim=c(0,1))
  mtext(extra_title, side = 1, adj = 1, padj = 5) ## SB
  abline( v=thresh_pchisq_loci, col='red')
  
  # Locussy fits
  opar <- par(mfrow=c(2, ncol( gpred) %/% 2), 
      mar=c( 3, 3, 0, 0)+0.1, oma=c( 2, 2, 3, 1))
  # omi=c(0,0,.6,0),mai=c(.8,.8,.2,.2)) Paige uses absolute margins
  on.exit( par( opar))
  
  gtypes <- colnames( gobs)
  for( g in 1:ncol( gpred)) {
    # all of 'em
    plot(gpred[,g], gobs[,g], pch=16, cex=0.4, col='lightblue', # pch='.' is too small
        xlab='', ylab='', main='', 
        xlim=c( 0, max( c( gpred[,g], gobs[,g]))),  
        ylim=c( 0, max( c( gpred[,g], gobs[,g]))))          
    abline(0,1,col=8,lwd=2)
    points( gpred[ keep_loci,g], gobs[ keep_loci,g], pch=16, cex=0.4, col='green') # overplot to make visible against the line
    mtext( side=3, gtypes[ g], line=-1)
    # bad ones with X
    points( gpred[ !keep_loci, g], gobs[ !keep_loci, g], col='magenta', pch=16, cex=0.6)
    if( any( iffy_loci)) {
      points( gpred[ iffy_loci, g], gobs[ iffy_loci, g], col='orange', pch=4, cex=1)
    }
  }
  mtext( 'Green = "good"', side=3, cex=1.5, outer=TRUE, line=1.5)
  mtext( 'Expected', side=1, cex=1.5, outer=TRUE)
  mtext( 'Observed', side=2, cex=1.5, outer=TRUE)
  mtext(extra_title, side = 1, adj = 1, padj = 4) ## SB
    
  if( liffies) {
    legend( 'bottomright', pch=c( 16, rep( 4, liffies)), col=c( 'magenta', if( liffies>1) 'orange', 'green'), pt.cex=c( 0.6, if( liffies>1) 1, 0.6),
        legend=c( sprintf( '%4.1e < Pr ...', thresh_pchisq_loci[ 1]), if( liffies>1) sprintf( '... < %4.1e', thresh_pchisq_loci[2]), '... < 1') )
  }
  
  if( trim) { 
    lociar <- lociar[ , keep_loci, ,drop=FALSE]
  }
return( lociar)
}


#' gtab6to4(): condense 6-way genotype counts to 4-way
#'
#' Condenses 6-way genotype counts to 4-way (e.g., gives 'AAO' instead of AA or AO)
#' @param gt6 a genotype object scored as 6-way.
#' @return a genotype object scored as 4-way.
#' @seealso gtab4
#' @export

"gtab6to4" <-
function( gt6) {
######### Condense 6way-genotype counts to 4way (AAO instead of AA, AO)
  define_genotypes()
  gt4 <- matrix( 0, nrow( gt6), 4, dimnames=list( dimnames( gt6)[[1]], genotypes4_ambig))
  gt4[,AB] <- gt6[,AB]
  gt4[,OO] <- gt6[,OO]
  gt4[,AAO] <- gt6[,AA] + gt6[,AO]
  gt4[,BBO] <- gt6[,BB] + gt6[,BO] 
return( gt4)
}



#' @export
"[.SPAgeno" <-
function( x, i, j, ..., drop=FALSE){
    x <- NextMethod( '[') 
    x <- prepare_PLOD_SPA( x)
    ## should only call PREPARE_PLOD_SPA if loci have changed
return( x) 
}

#' @export
"[<-.SPAgeno" <-
function( x, i, j, value) {

    x <- NextMethod( '[<-')
    x <- prepare_PLOD_SPA( x)
    ## should only call PREPARE_PLOD_SPA if loci have changed
return( x)  

}




# MVB's workaround for futile CRAN 'no visible blah' check:
globalVariables( package="kinference",
                names=c( ".Traceback"
                       ,"oa"
                       ,"pkgname"
                       ,"libname"
                       ,"g_err"
                       ,"genotypes_C"
                       ,"perr"
                       ,"g_1"
                       ,"g_2"
                       ,"AA_BB"
                       ,"A"
                       ,"B"
                       ,"snerr"
                       ,"AO_BO"
                       ,"O"
                       ,"pp_err"
                       ,"gg"
                       ,"g1_err"
                       ,"g1"
                       ,"g2_err"
                       ,"g2"
                       ,"pp_true"
                       ,"add_up"
                       ,"result"
                       ,"dimor"
                       ,"i"
                       ,"other"
                       ,"..."
                       ,"."
                       ,"OO"
                       ,"CC"
                       ,"CO"
                       ,"AO"
                       ,"AC"
                       ,"BO"
                       ,"BC"
                       ,"pp6_err"
                       ,"genotypes6"
                       ,"P"
                       ,"pr2"
                       ,"record"
                       ,"ABCO"
                       ,"is_het"
                       ,"gp1"
                       ,"gp2"
                       ,"pp3"
                       ,"su1u2"
                       ,"implied"
                       ,"swap12"
                       ,"swap34"
                       ,"cimplied"
                       ,"implstr"
                       ,"has_impl2"
                       ,"first_impl"
                       ,"second_impl"
                       ,"cimpl2"
                       ,"thing"
                       ,"j"
                       ,"ij"
                       ,"ijpairs"
                       ,"is_pair"
                       ,"test"
                       ,"x"
                       ,"y"
                       ,"xtest"
                       ,"to.do"
                       ,"chains"
                       ,"pairmats"
                       ,"this_chain"
                       ,"this_pairs"
                       ,"o"
                       ,"snpg"
                       ,"og"
                       ,"iwhat"
                       ,"LOD"
                       ,"PUP"
                       ,"PUPLOD"
                       ,"PUPLOD2"
                       ,"Kenv"
                       ,"mg"
                       ,"OLOD"
                       ,"use6"
                       ,"locinfo"
                       ,"use4"
                       ,"temp_snpg"
                       ,"recode4to6temp"
                       ,"AA"
                       ,"BB"
                       ,"n_loci"
                       ,"OPUP"
                       ,"NPUP"
                       ,"ig"
                       ,"gjseq"
                       ,"XXi"
                       ,"l"
                       ,"gj"
                       ,"Pg"
                       ,"gi"
                       ,"e_CLOD"
                       ,"%[gj]%"
                       ,"e2_CLOD"
                       ,"geno"
                       ,"g1seq"
                       ,"my_e_CLOD"
                       ,"%[g1]%"
                       ,"my_e2_CLOD"
                       ,"my_v_CLOD"
                       ,"my_rat_CLOD"
                       ,"%[l]%"
                       ,"SUM_"
                       ,"nsim"
                       ,"gsim"
                       ,"il"
                       ,"sim_e_CLOD"
                       ,"sim_e2_CLOD"
                       ,"sim_v_CLOD"
                       ,"sim_rat_CLOD"
                       ,"K"
                       ,"ETT"
                       ,"it"
                       ,"g"
                       ,"tt"
                       ,"log_S"
                       ,"%[g]%"
                       ,"dK"
                       ,"g12"
                       ,"LODOK"
                       ,"S"
                       ,"%[g12]%"
                       ,"SL"
                       ,"ddK"
                       ,"SLL"
                       ,"gbasics"
                       ,"bins"
                       ,"qq"
                       ,"nq"
                       ,"mean_theory"
                       ,"var_theory"
                       ,"symmo"
                       ,"subset1"
                       ,"big_PLOD"
                       ,"big_i"
                       ,"big_j"
                       ,"binprobs"
                       ,"eta"
                       ,"keep_thresh"
                       ,"hmfr"
                       ,"lglk"
                       ,"kinference"
                       ,"stat"
                       ,"mdull"
                       ,"boring"
                       ,"vdull"
                       ,"var"
                       ,"sstat"
                       ,"f"
                       ,"s"
                       ,"iv"
                       ,"dist"
                       ,"%[s1]%"
                       ,"%[s]%"
                       ,"s1"
                       ,"hist"
                       ,"abline"
                       ,"focusees"
                       ,"expr"
                       ,"obj"
                       ,"name"
                       ,"expand_dim"
                       ,"prefixdims"
                       ,"dimorlen"
                       ,"d"
                       ,"e"
                       ,"where"
                       ,"whoami"
                       ,"uij"
                       ,"n"
                       ,"m"
                       ,"nf"
                       ,"k"
                       ,"groups"
                       ,"keeps"
                       ,"drops"
                       ,"want_groups"
                       ,"nl"
                       ,"pobs"
                       ,"pbonzer"
                       ,"lglk_l"
                       ,"p6"
                       ,"ppar"
                       ,"pAO"
                       ,"pAA"
                       ,"pBO"
                       ,"pBB"
                       ,"pAB"
                       ,"pOO"
                       ,"AA2AO"
                       ,"AO2AA"
                       ,"AB"
                       ,"BB2BO"
                       ,"BO2BB"
                       ,"n_l"
                       ,"Nlglk"
                       ,"flush.console"
                       ,"pstart"
                       ,"fitto"
                       ,"nlminb"
                       ,"control"
                       ,"subset2"
                       ,"diplos"
                       ,"genotypes4_ambig"
                       ,"AAO"
                       ,"BBO"
                       ,"gtab"
                       ,"pid"
                       ,"DUP_paircomps_lots"
                       ,"max_diff_genos"
                       ,"big_similar"
                       ,"max_diff_ppn"
                       ,"OO_code"
                       ,"NA_geno"
                       ,"DUP_paircomps_incomplete_lots"
                       ,"limit"
                       ,"big_ndiff"
                       ,"big_ncomp"
                       ,"candiHSPs"
                       ,"sibg"
                       ,"just_sibg"
                       ,"PUP4"
                       ,"LOD4"
                       ,"PHSP4"
                       ,"P_k0"
                       ,"P_k1"
                       ,"P_k2"
                       ,"nsib"
                       ,"nloci"
                       ,"kappa_fsp"
                       ,"kappa_hsp"
                       ,"p12fsp"
                       ,"p12hsp"
                       ,"OD_FH"
                       ,"LOD_FH"
                       ,"PLOD_FH"
                       ,"OPHSP"
                       ,"p12fspa"
                       ,"p12hspa"
                       ,"EPLOD_FH_F"
                       ,"EPLOD_FH_H"
                       ,"ret"
                       ,"E_FSP"
                       ,"E_HSP"
                       ,"candiPOPs"
                       ,"just_snpg"
                       ,"n_pairs"
                       ,"pA"
                       ,"pB"
                       ,"pO"
                       ,"pgeno"
                       ,"off"
                       ,"Pr_same_given_k"
                       ,"Pr_nsame_FSP"
                       ,"Pr_same_FSP"
                       ,"%[k]%"
                       ,"Pr_same_POP"
                       ,"Pr_nsame_HSP"
                       ,"Pr_same_HSP"
                       ,"SD_FSP"
                       ,"SD_POP"
                       ,"SDwt_POP"
                       ,"SD_denom"
                       ,"wt"
                       ,"wtsame"
                       ,"E_POP"
                       ,"E_UP"
                       ,"inv_CDF"
                       ,"LOD6"
                       ,"temp_LOD"
                       ,"nbins"
                       ,"CDF"
                       ,"HSP_paircomps_lots"
                       ,"bits"
                       ,"make_CLOD"
                       ,"PHSP"
                       ,"Pg2_g1_H"
                       ,"e_CLOD_HSP"
                       ,"p0"
                       ,"pex"
                       ,"delta"
                       ,"SD"
                       ,"SD_combo"
                       ,"WPSEX_UP_POP_balance"
                       ,"V_combo"
                       ,"ww"
                       ,"pop_loci"
                       ,"AAish"
                       ,"BBish"
                       ,"HSP"
                       ,"UP"
                       ,"POP"
                       ,"FSP"
                       ,"rr"
                       ,"log1m_pexup"
                       ,"KK"
                       ,"retw"
                       ,"dKK"
                       ,"ddKK"
                       ,"n_sim_check"
                       ,"Ktest"
                       ,"runif"
                       ,"ewx"
                       ,"quick"
                       ,"pciles"
                       ,"POP_wt_paircomps_lots"
                       ,"n_in_bin"
                       ,"n_wpsex_in_bin"
                       ,"big_wpsex"
                       ,"snpg_i"
                       ,"snpg_j"
                       ,"isABOO"
                       ,"nABOO"
                       ,"oset"
                       ,"set"
                       ,"seed"
                       ,"newj"
                       ,"newi"
                       ,"v"
                       ,"target"
                       ,"edash"
                       ,"msqrt_v"
                       ,"median"
                       ,"use_loci"
                       ,"whmo"
                       ,"four_pab_poo"
                       ,"compaboo"
                       ,"etwab"
                       ,"etwoo"
                       ,"denom"
                       ,"dens_SPA"
                       ,"hist_pars"
                       ,"multhresh"
                       ,"lv"
                       ,"lines"
                       ,"mids"
                       ,"breaks"
                       ,"counts"
                       ,"badhetz_hist_pars"
                       ,"li"
                       ,"lociar"
                       ,"li1"
                       ,"temp0"
                       ,"cg6p0"
                       ,"temp1"
                       ,"cg6p1"
                       ,"g6p0"
                       ,"g6p1"
                       ,"s6"
                       ,"want_LOD_table"
                       ,"map6to4"
                       ,"mm"
                       ,"yup"
                       ,"%[j]%"
                       ,"g4p0"
                       ,"g4p1"
                       ,"s4"
                       ,"PUP6"
                       ,"n_samps"
                       ,"snpg4"
                       ,"info"
                       ,"Our_plate"
                       ,"Our_sample"
                       ,"Locus"
                       ,"lpgeno"
                       ,"ttp1"
                       ,"etp1l"
                       ,"num"
                       ,"num2"
                       ,"ntest"
                       ,"Ksim"
                       ,"meansim"
                       ,"genos"
                       ,"lp"
                       ,"ilglk"
                       ,"indiv_lglk_geno"
                       ,"ilglk_manual"
                       ,"indiv_lglk_hist_pars"
                       ,"mids_SPA"
                       ,"good_ind"
                       ,"elg"
                       ,"elg2"
                       ,"sdelg"
                       ,"GG"
                       ,"n_f"
                       ,"olg"
                       ,"%[f]%"
                       ,"sg"
                       ,"dips"
                       ,"quietly"
                       ,"warn"
                       ,"my_dll"
                       ,"onload_autowrap"
                       ,"pIBD1"
                       ,"Phsp"
                       ,"pIBD0"
                       ,"Pup"
                       ,"ngp"
                       ,"wanted"
                       ,"gpLOD"
                       ,"gpPUP"
                       ,"LOD_as_2D"
                       ,"PUP_as_2D"
                       ,"omg"
                       ,"double_wanted"
                       ,"what"
                       ,"E.HSP"
                       ,"%[i,j]%"
                       ,"E.UP"
                       ,"E2.UP"
                       ,"V.UP"
                       ,"Ediff"
                       ,"sdiff"
                       ,"retval"
                       ,"geno6"
                       ,"cn6"
                       ,"cn4"
                       ,"ichangio"
                       ,"was1"
                       ,"was2"
                       ,"now1"
                       ,"now2"
                       ,"iget4"
                       ,"hasO"
                       ,"make_K"
                       ,"n_pts_SPA_renorm"
                       ,"keeping"
                       ,"hi"
                       ,"lo"
                       ,"probinverts"
                       ,"one_in_X_eta"
                       ,"rough_n_pairs_to_keep"
                       ,"XX"
                       ,"ck"
                       ,"g12_code"
                       ,"this_g12"
                       ,"LOD_obs"
                       ,"PLOD"
                       ,"bin"
                       ,"kin.true"
                       ,"n.meio"
                       ,"linfo"
                       ,"L"
                       ,"count.l"
                       ,"e0"
                       ,"e0.l"
                       ,"e1"
                       ,"e1.l"
                       ,"v0"
                       ,"v0.l"
                       ,"v1"
                       ,"v1.l"
                       ,"V.noX"
                       ,"p"
                       ,"meioses"
                       ,"EofV"
                       ,"VofE"
                       ,"C"
                       ,"C.equiv"
                       ,"C.hat"
                       ,"emp.V.HSP"
                       ,"V0"
                       ,"ell"
                       ,"e2.1"
                       ,"e2.0"
                       ,"e1.0"
                       ,"e1.1"
                       ,"phi"
                       ,"mul"
                       ,"r"
                       ,"PSstar"
                       ,"P0star"
                       ,"si2"
                       ,"si3"
                       ,"si4"
                       ,"V.allX"
                       ,"rho"
                       ,"p11.8"
                       ,"p01.8"
                       ,"p00.8"
                       ,"p11"
                       ,"p01.10"
                       ,"p00"
                       ,"EL"
                       ,"EL2"
                       ,"VL"
                       ,"rho.hat"
                       ,"Vx"
                       ,"useN"
                       ,"g3p0"
                       ,"g3p1"
                       ,"genotypes3_ambig"
                       ,"keep_n"
                       ,"LOD3"
                       ,"PUP3"
                       ,"mean_HSP"
                       ,"n_PLODs_in_bin"
                        )
                )

