# This is package kinference 
#' @rawNamespace import( Rcpp)
#' @rawNamespace import( atease)
#' @rawNamespace import( mvbutils)
#' @rawNamespace import( gbasics)
#' @rawNamespace import( vecless)

".onLoad" <-
function( libname, pkgname) {
  ## This part should only kick in when debugging C code with VSCode--
  ## should do nothing otherwise
  oa <- system.file( sprintf( 'R/load_%s_dll.R', pkgname),
      package=pkgname, lib.loc=libname)
  if( !nzchar( oa)) {
    oa <- system.file( sprintf( 'R/load_%s_dll.r', pkgname),
        package=pkgname, lib.loc=libname)
  }

  if( nzchar( oa)) {
    source( oa, local=TRUE)
  }

  # This should always happen; defines R wrappers for dot-calls
  # The function is defined by Cloaders_kinference.R
  run_Cloaders_kinference()
}


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




#' PLOD threshold for HSPs
#' 
#' This function proposes a PLOD threshold for excluding almost all 3rd-order
#' kin, and computes the associated False-Negative Probability (i.e., that a
#' true HSP will have a PLOD below that threshold).
#' 
#' The rationale comes from fitting a mixture distribution to observed PLODs
#' within some range that is expected to contain only 2nd, 3rd, and
#' \emph{perhaps} a few 4th order kin. The threshold is then "chosen" (or
#' proposed; it's really up to you) so that the expected number of
#' false-positives from 3rd-order kin-pairs (i.e., with PLODs above the
#' threshold) matches whatever you decide. A histogram with expected values is
#' plotted (unless you tell it not to).
#' 
#' Means and variances of the mixture components are automatically set in
#' advance, so the mixture-fit only has to estimate the proportion of kin-pairs
#' of each type. The means are easily calculated from kinship coefficients,
#' based on allele frequencies. Variances, however, also depend strongly on the
#' degree of linkage between loci, and to some extent on the \emph{nature} of
#' the linkage (more chromosomes, or more crossovers?). This is handled
#' internally by the function \code{\link{var_PLOD_kin}} (qv), which uses the
#' observed "overdispersion" of PLODs for a subset of \emph{definite} 2nd-order
#' kin to place bounds on the variances of 3rd and 4th order kin (based on two
#' extreme assumptions about the \emph{nature} of linkage). The code of
#' \code{autopick_threshold} then explores different variances within those
#' bounds and
#' 
#' Despite the name, \emph{you} still have to supply sensible values for a
#' couple of parameters, based on looking at your data and understanding what
#' you are trying to do. So it's not \emph{completely} automated--- and never
#' will be! Choosing a threshold is \emph{not} an "optimization process" with
#' explicit bias/variance tradeoffs; rather, it's about ensuring that you have
#' adequate "engineering tolerance" in the next stage of CKMR. The
#' False-Negative Probability will, to a great extent, compensate for the
#' choice of threshold (i.e. removing any bias in the fitted CKMR model)
#' \emph{unless} you set the threshold too low, and thus qend up with some
#' 3rd-order kin-pairs in your set of "definite 2nd-orders".
#' \subsection{Fitrange and use4thIf you are using an HSP-oriented PLOD, then
#' the range of PLODs you fit to should extend from somewhere above 0 (which is
#' verrry close to the expected PLOD for 3rd-order kin), up to the RHS of the
#' HSP bump, but clearly not so far as to include any FSPs and POPs. If you set
#' the lower limit high enough, then you don't need to worry about 4th-order
#' kin intruding (their contribution would be negligible), so you can get away
#' with fitting a 2-component mixture (simpler, less to go wrong...) by setting
#' \code{use4th=FALSE}. But if you push the lower range closer to 0 (which does
#' give you a larger sample size for fitting), then you might need to set
#' \code{use4th=TRUE}. The (substantial) downside of doing that, is that there
#' are often more PLODs close to 0 than near the HSP mean, so the mixture-fit
#' (which has make more assumptions when also using 4th-orders--- and those
#' assumptions may not be perfect) will "concentrate its efforts" on getting a
#' good fit near 0, rather than near the 2nd-order mean which is what we really
#' need. It is worth experimenting. }
#' 
#' 
#' @param x a \code{\link{snpgeno}} or its \code{locinfo} attribute. Must
#' already have been prepared by running \code{hsp_power} (qv) and
#' \code{prepare_PLOD_SPA} (qv).
#' @param kin a dataframe of "close-ish" kin-pairs and their PLODs, presumably
#' from running \code{find_HSPs} (qv); must have a column "PLOD".
#' @param fitrange_PLOD two numbers, specifying the range of PLODs from
#' \code{kin} to use in \emph{fitting} (though all are \emph{plotted}, by
#' default)
#' @param FPtol_pairs how many expected False-Positive 3rd-order kin should the
#' threshold exclude?
#' @param use4th whether to allow for 4th-order kin when fitting.
#' @param selecto whether to choose the threshold based on the best mixture fit
#' ("ML"), or the most conservative ("paranoid").
#' @param NVAR how many variances to try, between the limits set by
#' \code{var_PLOD_kin}
#' @param plot_bins bin-width for histogram plotting. Default NULL means no
#' plot.
#' @param shading_density By default, all PLODs in \code{kin} will be included
#' in the histogram, even though only a subset are used in fitting. The
#' histogram bars \emph{not} used in fitting (i.e., below
#' \code{fitrange_PLOD[1]}) will be lightened in colour, according to this
#' parameter. Results are graphics-device-dependent, so you may need to
#' experiment away from the default; larger numbers usually mean lighter
#' shading. Setting \code{shading_density=NA} should result in a light
#' transparent rectangle covering the entire LHS of the graph, which you might
#' prefer. You can also set \code{xlim} as usual, to remove those left-hand
#' bars altogether.
#' @param want_all_results if TRUE, return dataframe(s) containing results for
#' each variance explored. This lets you examine "sensitivity".
#' @param ... other parameters passed to \code{hist}, eg \code{xlim},
#' \code{ylim}, \code{col}. Many others will be ignored, and some will cause
#' problems.
#' @return The proposed threshold, with lots of attributes. You can use those
#' to calculate False-Neg Probabilities for \emph{other} possible thresholds,
#' as per \bold{Examples}; you \emph{don't} have to accept the one that is
#' proposed here. Threshold choice is \bold{up to you} (and not the fault of
#' \code{\link{kinference}})!
#' @seealso \code{\link{hsp_power}}, \code{\link{var_PLOD_kin}}
#' @keywords misc
#' @examples
#' 
#' # Better have one...
#' 
#' @export autopick_threshold
"autopick_threshold" <-
function(
  x,
  kin,
  fitrange_PLOD,
  FPtol_pairs,
  use4th,
  selecto= c( 'ML', 'paranoid'),
  NVAR= 10,
  plot_bins= NULL,
  shading_density= 10,
  want_all_results= FALSE,
  ... # for plot
){
stopifnot(
    length( fitrange_PLOD)==2,
    all( is.finite( fitrange_PLOD))
  )
  fitrange_PLOD <- sort( fitrange_PLOD)
  min_PLOD <- fitrange_PLOD[ 1]
  max_PLOD <- fitrange_PLOD[ 2]

  li <- if( x %is.a% 'snpgeno') x$locinfo else x
  E_HSP <- sum( li$E.HSP)
  E_UP <- sum( li$E.UP)

  E2 <- E_HSP
  E3 <- (E_UP + E_HSP) / 2 # 0 is good-enuf approx *IF* stat...
  # ...is actually HSP::UP PLOD, but it _might_ be something else eg HSP::HTP
  E4 <- (3*E_UP + 1*E_HSP) / 4

  if( E3 > min_PLOD){
    warning( sprintf( 'fitrange_PLOD goes below 3rd-order mean, which is %5.2f; ' %&%
        'probably a bad idea', E3))
  }
  if( E2 > max_PLOD){
stop( sprintf( 'fitrange_PLOD exceeds 2nd-order mean, which is %5.2f; noooo!', E2))
  }

  kin_PLOD <- kin$PLOD %such.that% (. %in.range% fitrange_PLOD)
  overmean <- kin_PLOD %such.that% (. > E2)
  emp_V_HSP <- mean( sqr( overmean-E2))

  vpk <- var_PLOD_kin( x, emp_V_HSP, n_meio=c( 3, 4))

  CDF <- function( plod, P234, SD234){
      P234[1] * pnorm( plod, mean=E2, sd=SD234[1]) +
      P234[2] * pnorm( plod, mean=E3, sd=SD234[2]) +
      P234[3] * pnorm( plod, mean=E4, sd=SD234[3])
    }
  mixlglk3 <- function( P3, SD234){
      P234 <- c( 1-P3, P3, 0)
      pdf <<- (1-P3) * pdf2 +  P3 * pdf3
      Pr_in_range <<- diff( CDF( fitrange_PLOD, P234, SD234))
      pdf <- pdf / Pr_in_range
    return( sum( log( pdf)))
    }

  mixlglk4 <- function( Ppar, SD234){
    P3 <<- Ppar[1]
    P4 <<- (1-P3) * Ppar[2] # Ppar[2] is ppn of not-3s that are 4s
    P2 <<- 1 - P3 - P4
    pdf <<- P2*pdf2 + P3*pdf3 + P4*pdf4
    Pr_in_range <<- diff( CDF( fitrange_PLOD, c( P2, P3, P4), SD234))
    pdf <- pdf / Pr_in_range
  return( sum( log( pdf)))
  }

  set_stuff <- function( nlocal=sys.parent()) mlocal({
    # NB Evaluated directly in caller

    # We observe Nobs kin that are either 2nd or 3rd order
    # but iff they're within fitrange_PLOD. So the total number of 2nd+3rd(+4th) would be...
    Nall <- Nobs / Pr_in_range
    N3 <- Nall * P3 # just 3rds

    # Threshold is where "FPtol_pairs" of 3rds would be expected above it; 4ths irrel
    thresh <- qnorm( FPtol_pairs / N3, mean=E3, sd=SD3, lower.tail=FALSE)
    Pr_FNeg <- pnorm( thresh, mean=E2, sd=SD2)    
  })
  
  Nobs <- length( kin_PLOD)
  SD2 <- sqrt( emp_V_HSP)
  pdf2 <- dnorm( kin_PLOD, mean=E2, sd=SD2)
  pdf <- 0*pdf2 # placeholder
  Pr_in_range <- P2 <- P3 <- P4 <- (-999) # placeholder

  V3_range <- sort( vpk[ ,'M3'])
  SD3i <- seq( from=sqrt( V3_range[1]), to=sqrt( V3_range[2]), length=NVAR)
  
  # We may not use 4th-order, but prior to loop the compus are cheap
  V4_range <- sort( vpk[ ,'M4'])
  SD4i <- seq( from=sqrt( V4_range[1]), to=sqrt( V4_range[2]), length=NVAR)
  
  allvals3 <- allvals4 <- NULL

  for( ivar in seq_along( SD3i)){
    SD3 <- SD3i[ ivar]
    SD4 <- SD4i[ ivar]
    this_SD234 <- c( SD2, SD3, SD4)
    
    pdf3 <- dnorm( kin_PLOD, mean=E3, sd=SD3)
    pdf4 <- dnorm( kin_PLOD, mean=E4, sd=SD4)
    
    silly <- c( 0.01, 1-0.01)
    # If we hit either of these, it's not sensible...
    # min_silly means it's ALL HSPs, NO HTPs--- in which case, do by eye...
    # max_silly means NO HSPs
    bestio <- optimize( mixlglk3, silly, maximum=TRUE, SD234=this_SD234)
    P3 <- bestio$maximum
    if( min( abs( P3 - silly)) < 0.01){
      warning( sprintf( 'Silly "best" fit for SD3=%5.1f', SD3))
    }
    lglk <- mixlglk3( P3, this_SD234) # sets Pr_in_range and pdf
    set_stuff() # N3, thresh, etc

    # Add spurious (and silly) SD4: won't matter cos P4==0, but lets CDF calc work OK
    newvals3 <- returnList( SD2, SD3, SD4=0.1, P2=1-P3, P3, P4=0, Nall, lglk, thresh, Pr_FNeg)
    allvals3 <- if( is.null( allvals3)) data.frame( newvals3) else
        rbind( allvals3, newvals3)
        
    if( use4th){
      startio <- c( 0.1, 0.7)
      # 4ths likely 2--4X 2nds 
      bestio <- nlminb( startio, NEG( mixlglk4), 
          lower=c( 0, 1.5/(1+1.5)), upper=c( 1, 5/(1+5)),
          SD234= this_SD234)
      lglk <- mixlglk4( bestio$par, this_SD234)
      set_stuff()
      newvals4 <- returnList( SD2, SD3, SD4, P2, P3, P4, N3, Nall, lglk, thresh, Pr_FNeg)
      allvals4 <- if( is.null( allvals4)) data.frame( newvals4) else
          rbind( allvals4, newvals4)
    }
  }

  selecto <- match.arg( selecto)
  
  allvals <- if( use4th) allvals4 else allvals3
  
  i_highest <- which.max( allvals$thresh)
  i_fittest <- which.max( allvals$lglk)
  picki <- if( selecto=='paranoid') i_highest else i_fittest

  # Graph?
  if( !is.null( plot_bins)){
    # Taken from HSP_histo (before its renaming)
    # X-range goes from lowest *observed* PLOD (in kin), to *chosen* max_PLOD
    # CDF
    histo <- hist( kin$PLOD %such.that% (. < max_PLOD),
        breaks=seq( from= min( kin$PLOD), to= max_PLOD, by= plot_bins),
        col="lightgrey",xlab="PLOD", 
        main = sprintf( 'Autothresh: %s, #FP=%5.1f', selecto, FPtol_pairs),
        ...)

    # Too confusing to have >1 fit on same plot, so...
    # with() next allows direct use of that row of vals
    with( as.list( allvals[ picki,]), {
      # Composite distro
      abline( v=thresh, lty=2, col='blue')
      probblies <- Nall * diff( CDF( histo$breaks, c( P2, P3, P4), c( SD2, SD3, SD4)))
      lines( histo$mids, probblies, col='violet', lty=1)

      # Might as well see the components...
      for( ord in c( 2, 3, if( use4th) 4 else NULL)){
        Eo <- get( 'E' %&% ord)
        SDo <- get( 'SD' %&% ord)
        Po <- get( 'P' %&% ord)
        
        probbly <- diff( pnorm( histo$breaks, mean=Eo, sd=SDo)) * Nall * Po
        lines( histo$mids, probbly, col='violet', lty=2) 
      } # for ord
    }) # with picki

    # Shade out region _not_ used in fitting.
    if( dev.capabilities( 'semiTransparency')[[1]] && is.na( shading_density)){
      DENSITY <- NA
      COL <- gray( 0.9, alpha=0.8) # almost white
    } else {
      DENSITY <- shading_density * par( 'pin')[1]
      COL <- 'white'
    }
    rect( par( 'usr')[1], 0, min_PLOD, par( 'usr')[4], density=DENSITY, border=NA,
        col= COL)
        
    # Means: show them last, so not covered by shading
    abline( v=c( E2, E3, E4), col='black', lwd=3) # E4 *shouldn't* appear; should be off LHS!        
    
    legend("topright", legend = c( 'Theory means HSP & HTP', 'Overall', 'Component'),
      lwd= c( 3, 1, 1), lty= c( 1, 1, 2),
      col= c( 'black', 'violet', 'violet'), bg='white')
  }

  flatto <- function( matto){
      # Add row & col names to matrix that is getting vectorized
      m <- c( matto)
      names( m) <- outer( rownames( matto), colnames( matto), 
          function( x, y) sprintf( '%s_%s', x, y))
    return( m)
    }

  threshold <- allvals$thresh[ picki]
  attributes( threshold) <-  c(
      returnList( fitrange_PLOD, selecto, use4th, FPtol_pairs),
      info=list( c(
         flatto( vpk[,-1]),
         vpk@info %without.name% "n_meio",
         unlist( returnList( E2, E3, E4)),
         unlist( allvals[ picki, ] %without.name% 'lglk'))
        )
    )
  if( want_all_results){
    # Express lglks relto best
    MAX_lglk <- max( allvals$lglk)
    allvals$lglk <- allvals$lglk - MAX_lglk
    allvals3$lglk <- allvals3$lglk - MAX_lglk
    if( use4th){
      allvals4$lglk <- allvals4$lglk - MAX_lglk
    }
  
    threshold@allvals3 <- allvals3
    threshold@allvals4 <- allvals4
  }
  
return( threshold)
}


"calc_g6probs" <-
function( pA, pB, pC, snerr) {
  # snerr = P( misclassifying true XX as XO, and vice versa)
  define_genotypes()
  p0really <- 1 - (pA+pB+pC)
  pO <- p0really + pC

  pAA <- pA^2
  pAB <- 2*pA*pB
  pBB <- pB^2
  pAC <- 2*pA*pC
  pBC <- 2*pB*pC

  pA0 <- 2*pA*p0really # really true null called 0
  pB0 <- 2*pB*p0really
  pOO <- pO^2

  # Calcs assume matrix (ie Xtuple loci); coerce if not, then undim later
  if( no_dim <- is.null( dim( snerr))) {
    snerr <- matrix( snerr, 1, ncol=length( snerr), dimnames=list( NULL, names( snerr)))
  }

  # Error possible iff (i) true AA or (ii) A present C not

  pAO <- pAC + pA0 * (1-snerr[,'AO2AA']) + pAA * snerr[,'AA2AO']
  pAA <- pAA * (1-snerr[,'AA2AO']) + pA0 * snerr[,'AO2AA']

  pBO <- pBC + pB0 * (1-snerr[,'BO2BB']) + pBB * snerr[,'BB2BO']
  pBB <- pBB * (1-snerr[,'BB2BO']) + pB0 * snerr[,'BO2BB']

  p6 <- do.call( 'cbind', FOR( genotypes6, get( 'p' %&% .)))
  perr <- cbind(
      AAO=pA0 * snerr[,'AO2AA'] + pAA * snerr[,'AA2AO'],
      BBO=pB0 * snerr[,'BO2BB'] + pBB * snerr[,'BB2BO'])
  colnames( p6) <- genotypes6
  if( no_dim) {
    p6 <- p6[1,]
    perr <- perr[1,]
  }

  p6@perr <- perr
return( p6)
}


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

  pp6_err <- pp6_err # MVB was paranoid that this might be needed for some reason. Harmless...
return( pp6_err)
}


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

  pp6_err <- pp6_err # MVB was paranoid that this might be needed for some reason. Harmless...
return( pp6_err)
}


"calc_g6probs_IBD2_scalar" <-
function( P, snerr, record=FALSE) {
## SCALAR-ONLY VERSION... this is hard enough!
## Though can be called with 1-row matrix args, eg with( x@locinfo[1,],
##    calc_g6probs_IBD2( pbonzer, snerr))

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

  # UP case was:
  # pp_true[ cbind( gp1, gp2)] <- pr2[ gp1] * pr2[ gp2]                       ? 0
  pp_true[ cbind( genotypes_C, genotypes_C)] <- pr2[ genotypes_C]             ? 0

  # NB that for this UP case, XX/XO errors shouldn't change the overall probs because the cutoffs are chosen to do exactly that!
  add_pairprob_error()
  pp6_err <- pp6_err # MVB was paranoid that this might be needed for some reason. Harmless...

return( pp6_err)
}


"calc_hspPower_checksum" <-
function( li){
  Nseq <- 1 %upto% nrow( li)
return(  with( li, round( sum(
    Nseq %**% cbind( pbonzer[,1], snerr, useN) ),
    digits=10)
  ))
}


"calc_PLODSPA_checksum" <-
function( li){
  Nseq <- 1 %upto% nrow( li)
return(  with( li, round( sum(
    Nseq %**% cbind( useN, LOD6, LOD4, PUP6, PUP4, LOD3, PUP3)),
    digits=10)
  ))
}


"calculate_IBD" <-
function(lociar){

  define_genotypes()
  li <- lociar@locinfo
  li1 <- li[1,]

  temp0 <- with( li1, calc_g6probs_IBD0_scalar( pbonzer, snerr, record=TRUE))
  cg6p0 <- make_playback( calc_g6probs_IBD0_scalar, temp0)

  temp1 <- with( li1, calc_g6probs_IBD1_scalar( pbonzer, snerr, record=TRUE))
  cg6p1 <- make_playback( calc_g6probs_IBD1_scalar, temp1)

  temp2 <- with( li1, calc_g6probs_IBD2_scalar( pbonzer, snerr, record=TRUE))
  cg6p2 <- make_playback( calc_g6probs_IBD2_scalar, temp2)

  pIBD0 <- with( li, cg6p0( pbonzer, snerr))
  pIBD1 <- with( li, cg6p1( pbonzer, snerr))
  pIBD2 <- with( li, cg6p2( pbonzer, snerr))

  return(list(pIBD0 = pIBD0,
              pIBD1 = pIBD1,
              pIBD2 = pIBD2
              ))
}


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




#' Find chains in HSPs; summarize sib-groups
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
#' @aliases chain_pairwise get_chain
#' @param thing output from \code{find_HSPs} or \code{find_POPs} etc, or some
#' subset thereof
#' @param seed one sample ID, interpreted as a row-number in \code{thing}. To
#' do:also allow names, via \code{info} attr.
#' @return \code{chain_pairwise} returns a list of matrices, each for one
#' chain; the rows and columns of each matrix are the samples in that chain. A
#' "+" in the matrix indicates that those two samples have a direct pairwise
#' link (i.e., they appear together in one row of \code{thing}); a "." means
#' not. The rows and columns of each matrix are sorted so that the linkiest
#' samples are on the bottom and right. \code{get_chain} returns the row-subset
#' of \code{thing} that is chained to \code{seed}.
#' @keywords misc
#' @export chain_pairwise
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
}




#' QC for kin-finding; private for now
#' 
#' Incomplete! Suppose to return predicted mean & variance of CLODs for each
#' sample, ie how prone is that sample's particular genotype to yielding
#' unusually high/low PLODs when compared with a random unrelated sample. Then
#' you can turn this into a prediction of each sample's per-comp chance of
#' yielding a False-Positive with another Unrelated sample. This looks like
#' quite a powerful diagnostic, but is not fully explored yet. (Update: we used
#' this when "chasing Weird Johnny and the Shoulder of Doom". It worked... but
#' entirely \emph{not} in the way we expected, revealing a rather different
#' cause..!) The document
#' "d:/docs/genetics/Dart/sbt-baits-v3/too-many-plods.lyx" has more info in
#' section 4.1 on "rat CLODs".
#' 
#' There is a bunch of code in the function connected with simulations and more
#' elaborate calculations, currently commented out. So I'm "internalling" this
#' for now. Let's hope you're not able to see this in the documentation...
#' 
#' @param snpg a \code{snpgeno} object.
#' @param nsim currently inactive. A simulation option exists in the code to
#' check the null distro (not much use for far tails, of course).
#' @return Dataframe with columns "ECLOD" and "VCLOD". See examples format.
#' @keywords misc
#' @examples
#' 
#' ## Rough chance of yielding a PLOD>5, say
#' # cloddo <- check_FPosity( snpg = snpg)
#' # Pr_FPos_5 <- pnorm( 5, mean=cloddo$ECLOD, sd=sqrt( cloddo$CLOD), lower=FALSE)
#' # hist( Pr_Fpos_5, nc=50)
#' ## highlight some known suspects
#' # abline( v=Pr_Fpos_5[ suspects], col='red')
#' 
#' @export check_FPosity
"check_FPosity" <-
function( snpg, nsim=0){
## snpg should have been thru 'prepare_PLOD_SPA' so it has @PPS
stopifnot( 'Kenv' %in% names( attributes( snpg)))

  og <- options( vecless.print=FALSE)
  on.exit( options( og))

  define_genotypes()
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

  attributes( result) <- c( attributes( result),
      returnList( bins, binprobs, eta, keep_thresh))
  result@call <- sys.call()

return( result)
}




#' Locus QC check
#' 
#' Checks 6-way and 4-way genotype frequencies against HWE expectations, and
#' generates plots of observed / expected frequencies. Recently moved into
#' kinference from genocalldart.
#' 
#' 
#' @param geno6 a \code{snpgeno} object with 4-way and 6-way genocalls
#' @param thresh_pchisq_6and4 thresholds for \code{bad} and \code{really bad}
#' p-values
#' @param return_what one of \code{just_pvals} or \code{all}; see value
#' @param extra_title a character string to be added to the bottom-right corner
#' of all plots. Best if < 25 characters.
#' @return Creates per-locus vectors \code{pval6} and \code{pval4} for 6-way
#' and 4-way genotypes respectively. If \code{return_what="just_pvals"}, these
#' are returned in a list; if \code{return_what="all"}, they are added as
#' columns to \code{geno6$locinfo}.
#' @keywords misc
#' @export check6and4
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




#' Check observed genotypes against HWE expectations
#' 
#' Checks observed genotype frequencies against expected frequencies,
#' presumably with expectation defined by HWE.
#' 
#' 
#' @param lociar a snpgeno object
#' @param gpred predicted allele frequencies
#' @param gobs observed allele frequencies
#' @param thresh_pchisq_loci a param. Presumably, a threshold p-val for
#' flagging loci with suspicious-looking allele frequencies.
#' @param test a character string, either "Pearson" or "G"
#' @param trim TRUE or FALSE. TRUE will keep only above max thresh_pchisq_loci.
#' Arguably better as to be done post-hoc.
#' @param seq_paxis numeric#'
#' @param extra_title a character string to be added to the bottom-right corner
#' of all plots. Best if < 25 characters.
#' @seealso kinference::check6and4
#' @keywords misc
#' @export chisq_genofreq_check
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




#' Grouping pairwise duplicates
#' 
#' Constructs equivalence classes to show from pairwise equivalences, and
#' returns the "surplus" elements; if you then drop those elements, only one
#' element from each eq-class will be retained.
#' 
#' Input should be row numbers in a \code{snpgeno} objects of duplicates, as a
#' two-column data.frame or matrix with each row being a pair of duplicates, or
#' the output from \code{find_duplicates()} (a 3-col matrix). Identifies
#' \code{groups} of equivalent observations (e.g., if i and j are duplicates,
#' and j and k are duplicates, then i, j, and k are all equivalent). Outputs a
#' vector of the row numbers for all-but-one of each group.
#' 
#' @param ij 2-column matrix or data.frame; probably "row numbers" in a
#' dataset, though might work with character strings too
#' @param want_groups if \code{TRUE}, also return the equivalence-classes
#' themselves, as attribute \code{groups}.
#' @return Surplus elements in \code{ij}, perhaps plus attributes \code{groups}
#' if \code{want_groups=TRUE}. You can look at that to figure out which
#' elements are being retained (one "representative" from each equiv class).
#' @keywords misc
#' @examples
#' 
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
#' 
#' @export drop_dups_pairwise_equiv
"drop_dups_pairwise_equiv" <-
function( ij, want_groups=FALSE) {
  if(ncol(ij) == 3) {
    ij <- ij[,2:3]
  }  ## so that users don't have to specify it every time

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


"DUP_paircomps_incomplete_lots" <-
function(geno1, geno2, symmo, max_diff_ppn, limit) {
    .Call(`_kinference_DUP_paircomps_incomplete_lots`, geno1, geno2, symmo, max_diff_ppn, limit)
}


"DUP_paircomps_lots" <-
function(geno1, geno2, symmo, max_diff_loci, keep_n, nbins, binterval) {
    .Call(`_kinference_DUP_paircomps_lots`, geno1, geno2, symmo, max_diff_loci, keep_n, nbins, binterval)
}




#' Estimate ALFs from 6-way genotypes and snerr
#' 
#' Performs "straight" estimation of ALFs, given 6-way genotypes and snerr.
#' Won't allow for changes in C-allele frequency from one population to the
#' next. In principle, should just use
#' \code{genocalldart::choose_geno6_thresholds} but fix the count-related
#' thresholds and re-estimate ALFs. (What??)
#' 
#' 
#' @param snpg a \code{snpgeno} object with \code{snerr} included
#' @param control as per \code{nlminb}
#' @keywords misc
#' @export est_ALF_6way
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
}




#' Estimate allele frequencies, including nulls
#' 
#' Uses "ABCO" genotypes, ie up to 3 scorable alleles plus possible nulls, eg
#' from \code{geno_deambig_ABC} (qv). NALF (null allele frequency) is estimated
#' from HWE deviations, so this requires a decent sample size. But, it doesn't
#' require the elaborate 6-way genotyping and the massive read-depths needed
#' for that.
#' 
#' 
#' @param lociar a \code{loc.ar} object with the @geno_amb attribute.
#' @return Returns the input, adding a 4-column matrix \code{pambig} to the
#' "locinfo" attribute, plus attributes \code{gobs} and \code{gpred} showing
#' observed and expected counts of each genotype per locus.
#' @seealso geno_dembig_ABC
#' @keywords misc
#' @export est_ALF_ABCO
"est_ALF_ABCO" <-
function( lociar, geno_amb = lociar@geno_amb) {
########## Taken largely from "pipeline_for_SBT_baits.r"
########## MVB: I'd like to clean this up
########## Careful "parallel Newton-Raphson" could allow vectorization and whoosh-factor, but NFN I guess

  define_genotypes() # AAO etc
  ## geno_amb <- lociar@geno_amb
  if( is.null( geno_amb)) {
stop( "No 'geno_amb' attribute :(")
  }
    stopifnot( is.character( geno_amb) || my.all.equal( geno_amb@diplos, genotypes_ambig))

    "inv.logit" <- function(x) plogis(x)

  expected <- NULL
  lglk <- function( params, nobs, return_expected=FALSE) {
      has_C <- length( params)==3

      # Reparamed for with-C case to logit scale, to avoid probs when pC~=0
      pA <<- inv.logit( params[1])
      pB <<- (1-pA) * inv.logit( params[2])
      pC <<- if( has_C) (1-pA-pB) * inv.logit( params[ 3]) else 0
      pO <<- max( 0, 1 - pA - pB - pC) # rounding error guard

      phat <- make_pgeno( pA, pB, pC, which_genotypes=genotypes_ambig)
      expected <<- n_fish * phat
      lglk <- nobs %*% log(phat + (nobs==0))        # avoid 0log0 gotcha
      pen <- penscale * sum( log( cosh( params-start_par)))
    return( lglk - pen)
    }
  pA <- pB <- pC <- pO <- (-1) # overwritten when lglk runs

  n_fish <- nrow( geno_amb)
  n_loci <- ncol( geno_amb)
  gobs <- gpred <- matrix( 0, n_loci, length( genotypes_ambig), dimnames=list( NULL, genotypes_ambig))
  for( g in genotypes_ambig) {
    gobs[,g] <- colSums( geno_amb==g)
  }

  # MVB: Old code looks pretty slow.
  # Calc g.freq for all loci **BUT ONLY ACCEPTABLE FISH** to preserve matrix size
  # g.freq <- apply(gABO.obs[!iamb.f, ], 2, function(x) table(factor(x,levels=c("AA","AB","BB","OO"))))

  pambig_est <- matrix(NA, n_loci, 4, dimnames=list( NULL, cq( A, B, C, O)))   # nloci rows, 2 cols (pA, pB) (p0 = 1-rowSums(p.est))
  conv <- rep( NA, n_loci) # convergence diagnostic

  tiny <- 2^-12 # avoid rounding error

  scatn( 'Starting ambig-geno_amb MALF/NALF estimation on %i loci:\n', n_loci)
  evalq( # for debug speed
  for( ll in 1:n_loci)  {
    if( ll %% 50 == 0) { cat( '\r', ll); flush.console() }
    has_C <- sum( gobs[ll,cq( AC, BC, CCO)])>0

    # Rough ests based on presence: mild overflow guard
    pA <- 1 - sqrt( 1- sum( gobs[ll, cq( AAO, AB, AC)]) / (1 + n_fish))
    pB <- 1 - sqrt( 1- sum( gobs[ll, cq( BBO, AB, BC)]) / (1 + n_fish))
    pC <- if( !has_C) 0 else 1 - sqrt( 1- sum( gobs[ll, cq( CCO, AC, BC)]) / (1 + n_fish))

    # Really, loci with rubbish pB (or pA) should have been chucked by now... but just in case...
    if( pA==0) pA <- 1/(2*n_fish)
    if( pB==0) pB <- 1/(2*n_fish)


    # Hard overflow guard...
    duhhh <- pA + pB + pC
    if( duhhh > 0.99) {
      duhhh <- duhhh + 0.011 # push it over 1
      pA <- pA / duhhh
      pB <- pB / duhhh
      pC <- pC / duhhh
    }

    "logit" <- function(p) qlogis(p)

    start_par <- c( logit( pA), logit( pB / (1-pA)), if( has_C) logit( pC / (1-pA-pB)))

    # Set reasonable penalty scale
    penscale <- 0
    testo <- numeric( 3)
    for( i in 1:3) {
      testo[ i] <- lglk( start_par+c( -1, 0, 1)[i], nobs=gobs[ll,])
    }
    penscale <- max( abs( diff( testo))) / 1e4

    fit <- nlminb( start_par, NEG( lglk), nobs=gobs[ll,])
    conv[ll] <- fit$convergence
    besto <- fit$par

    # Try refit with reduced and recentred penalty
    start_par <- besto
    penscale <- penscale / 10
    fit <- nlminb( start_par, NEG( lglk), nobs=gobs[ll,])

    # Make sure 'expected' is up-to-date
    lglk( fit$par, nobs=gobs[ll,])
    gpred[ll,] <- expected
    pambig_est[ll,] <- c( pA, pB, pC, pO)
  }) # for, evalq

  scatn( 'Convergence result table (0 is ideal)')
  print( table( conv))

  # Gene frequency for subsequent analysis (pA, pB, p0)
  lociar@locinfo$pambig <- pambig_est
  lociar@gobs <- gobs
  lociar@gpred <- gpred # for chi-sq checks
  lociar@subset_like_loci <- c( lociar@subset_like_loci, 'gobs', 'gpred')

return( lociar)
}


"est_ALF_ABO_quick" <-
function( x=NULL, AB, AAO, BBO, OO, tol=1e-7, EMtol=1e-3, quietly=FALSE){
## Multilocus A/B/O freq estimation from 4way genotypes, 
## with nulls obvs but assuming neglig geno _error_.
## Either from a 'snpgeno' or similar (currently must be 4-way), 
## in which case 'pbonzer' gets added to 'locinfo';
## or as direct entries of totals.
## EM algo is very simple here, allowing vectorization
## To speed up the notorious EM, an outer iteration of Aitken accel
## Vectorizing makes it look ugly (in R anyway) but it is bloody fast...

  if( !is.null( x)){ 
stopifnot( missing( AB), missing( AAO), missing( BBO), missing( OO))  
    gbasics::define_genotypes()
stopifnot( my.all.equal( x@diplos, genotypes4_ambig))
    AB <- colSums( x=='AB')
    AAO <- colSums( x=='AAO')
    BBO <- colSums( x=='BBO')
    OO <- colSums( x=='OO')
  }
  
  n <- AB + AAO + BBO + OO
  LLL <- length( AB)

  # Need some starting values...
  omega <- 0.05 + 0.95 * sqrt( OO/n) # just tame it a bit
  # alpha/beta we will do as if no nulls... then scale to non-null total
  # Terrible, but bounded!
  gamma <- (2*AAO + AB) / (2*BBO + AB)
  beta <- (1-omega) / (1+gamma)
  alpha <- 1 - beta - omega

  # Scale the lot...
  i2n <- 1/(2*n)
  AB <- AB * i2n
  AAO <- AAO * i2n
  BBO <- BBO * i2n
  OO <- OO * i2n

  nitsi <- AO <- AA <- BO <- BB <- 0*AB
  
  # only iterate unconverged ones. tRickeRy..!
  prev_prev_omega <- prev_omega <- 0*AB -1 # two converged iterations in case lucky start

  # Next 2 not used in convergence checks: just for Aitken accel below
  prev_prev_alpha <- prev_alpha <- alpha
  prev_prev_beta <- prev_beta <- beta
  Aitken <- function( x0, x1, x2){
      d0 <- x1 - x0
      d1 <- x2 - x1
      dd0 <- d1 - d0

      res <- x0 - sqr( d0) / dd0
      res[ !is.finite( res)] <- x0[ !is.finite( res)]
    return( res)
    }

  iAit <- seq_along( AB)
  nits <- 0
  n_super_its <- 0 # for Aitkening
  MAX_AITKEN <- 20 # few will take more than 10
  while( n_super_its < MAX_AITKEN){
    aitA <- alpha
    aitB <- beta
    aitO <- omega # ... for superconvergence

    # Which elements will be worked on?
    i <- iAit # to start with
    i <- i[ ((AB+AAO)[iAit]>0) & ((AB+BBO)[iAit]>0)] # other cases have no alleles for A and/or B!
  
    repeat{
      nits <- nits+1
      nitsi[ i] <- nitsi[ i] + 1
      AA[ i] <- AAO[ i] * sqr( alpha[ i]) / (sqr( alpha[ i]) + 2*alpha[ i]*omega[ i])
      AO[ i] <- AAO[ i] - AA[ i]

      BB[ i] <- BBO[ i] * sqr( beta[ i]) / (sqr( beta[ i]) + 2*beta[ i]*omega[ i])
      BO[ i] <- BBO[ i] - BB[ i]

      alpha[ i] <- AB[ i] + 2*AA[ i] + AO[ i]
      beta[ i] <- AB[ i] + 2*BB[ i] + BO[ i]
      omega[ i] <- AO[ i] + BO[ i] + 2*OO[ i]

      # Check 1 step AND 2 steps back... in case of bouncing!
      is_dun <- EMtol > pmax( 
          abs( omega[ i] - prev_omega[ i]),
          abs( omega[ i] - prev_prev_omega[ i])
        )
      if( (nits>3) && all( is_dun)){
    break
      }

      i <- i[ !is_dun]
      prev_prev_omega[ i] <- prev_omega[ i]
      prev_omega[ i] <- omega[ i]
      prev_prev_alpha[ i] <- prev_alpha[ i]
      prev_alpha[ i] <- alpha[ i]
      prev_prev_beta[ i] <- prev_beta[ i]
      prev_beta[ i] <- beta[ i]    
    } # normal EM iteration

    # Not gonna over-subscript here; these are so quick
    ACC_alpha <- Aitken( alpha, prev_alpha, prev_prev_alpha)[ iAit]
    ACC_beta <- Aitken( beta, prev_beta, prev_prev_beta)[ iAit]
    ACC_omega <- Aitken( omega, prev_omega, prev_prev_omega)[ iAit]
    
    # Force a few regular EM steps in next superit
    prev_omega[] <- (-1)
    prev_prev_omega[] <- (-2) 
    
    # That accelerates the terms *individually*, but scaling is crucial, so...
    OK <- (ACC_alpha>=0) & (ACC_beta>=0) & (ACC_omega>=0)
    iACC_sum <- 1/(ACC_alpha + ACC_beta + ACC_omega)
    # scatn( 'Nits: %i, Nsuper: %i', nits, n_super_its)
    oa <- alpha
    alpha[ iAit[ OK]] <- (ACC_alpha * iACC_sum)[OK]
    # scatn( 'Alpha change')
    # print( (oa-alpha)[ iAit])
    
    beta[ iAit[OK]] <- (ACC_beta * iACC_sum)[OK]
    omega[ iAit[OK]] <- (ACC_omega * iACC_sum)[OK]
    
    # Next will not check OK during update; but that's OK :) cos then it relies entirely on EM...
    # ... if EM itself has converged, then we are done
    iAit <- which( abs( alpha - aitA) + abs( beta - aitB) + abs( omega - aitO) > tol)
    
    if( !length( iAit)){
  break
    }
    
    n_super_its <- n_super_its + 1
  } # Aitken iteration...

  if( !quietly){
    scatn( "Summary of iterations:")
    print( summary( nitsi))
  }


  if( length( iAit)){
warning( sprintf( "Still %i not-fully-converged loci after maximum [%i] Aitken superloops",
      length( iAit), MAX_AITKEN))
  }

  if( missing( x)){
    mat <- cbind( alpha, beta, omega)
    mat@nits <- nits
return( mat)
  } else {
    x@locinfo$pbonzer <- matrix( c( alpha, beta, 0*beta, omega), ncol=4, 
        dimnames=list( NULL, cq( A, B, C, O)))
return( x)
  }
}




#' Kin-finders for loads-of-SNPs datasets
#' 
#' These take a \code{snpgeno} dataset that has been processed as far as
#' \code{\link{check6and4}} (and for HSPs, \code{\link{prepare_PLOD_SPA}}) and
#' find various relations between the samples. Relationships include duplicates
#' (DUPs/dupes), parent-offspring pairs (POPs) and half-sibling pairs (HSPs),
#' plus of course unrelated pairs (UPs). You can specify the same or different
#' subsets of the \code{snpgeno} for comparison: e.g., first subset for the
#' adults, second for the juveniles.
#' 
#' Some categories will "catch" others (eg \code{find_HSPs} will certainly
#' include any POPs too), so you may need the splitter routines such as
#' \code{split_POPs_from_HSPs} afterwards. The safest general-purpose strategy-
#' but often \emph{not} the most sensible, if your data is nicely organized and
#' you know what you want- would be:
#' 
#' \code{find_duplicates} and then get rid of them
#' 
#' \code{find_HSPs} to get \emph{all} kin (though you will usually have to
#' sacrifice some HSPs to false-neg because you'll need a threshold)
#' 
#' \code{split_POPs_from_HSPs} to split HSPs from POPs/FSPs
#' 
#' \code{split_POPs_from_FSPs} to split the latter.
#' 
#' The non-splitter functions, i.e. \code{find_XXX}, might be run on huge
#' numbers of samples, entailing a huge^2 number of comparisons. You don't want
#' all those individual comparison results, and your computer certainly
#' wouldn't enjoy trying to keep them! So the general idea is to set a
#' threshold for what constitutes "maybe worth keeping individually" (that you
#' expect will be generous enough to contain everything you \emph{do} want,
#' plus some dross), and then to retain just binned counts of the relevant comp
#' statistic for all comps (usually, the vast majority) which don't make your
#' threshold.
#' 
#' In addition, the \code{limit_pairs} argument is there to prevent your
#' computer locking out with bazillions of unwanted pairs (in case you guess
#' the bin limit inapproriately); the comparisons will be stopped if
#' \code{limit_pairs} is hit, with a warning. In that case, you probably need
#' to change a threshold, or re-run with larger \code{limit_pairs}. The default
#' isn't meant to correspond to any biomathematical logic, it's just to stop
#' blue smoke coming out your USB ports.
#' 
#' For \code{find_duplicates}, there are at least two different use-cases.
#' First, you might want an initial run on a non-too-large subset of your data,
#' to check that dups \emph{can} be clearly distinguished and to look at
#' typical extent of genotyping errors (based on clear duplicates that don't
#' match at every locus). For that, you can set \code{nbins} and choose some
#' reasonable guess as to \code{max_diff_loci} (say, 5\% of the number of
#' loci). Because you set \code{nbins>0}, \emph{every} pair (almost...) gets
#' checked at \emph{all} loci, so it can be slow. Thus, if you have done this
#' before and have a good sense of "how bad can a real duplicate be?", then set
#' \code{nbins=0} (and \code{max_diff_geno} to a small but safe value that
#' won't miss any realistic duplicate-with-genotyping-error) so it will abort a
#' comparison early as soon as it reaches \code{max_diff_geno} differing loci.
#' That saves a \emph{lot} of time on big datasets! You won't get a histo of
#' number-of-diffs, but you don't need one for that use-case. The "almost" is
#' that \code{find_duplicates} uses "transitivity" (if A is a dup of B and of
#' C, then we don't need to check B vs C), so it only counts differences for
#' not-yet-known duplicates \emph{based on} \code{max_diff_loci}. To discard
#' duplicates and to find entire equivalence-classes of duplicates, e.g. from a
#' control specimen included in numerous plates, see
#' \code{drop_dup_pairwise_equiv}.
#' 
#' @aliases find_duplicates find_HSPs find_POPs %upto%
#' @param snpg a \code{snpgeno} object
#' @param subset1,subset2 numeric vectors of which samples to use (not logical,
#' not negative). Defaults to all of them. Iff \code{subset1} and
#' \code{subset2} are identical, only half the comparisons are done (i.e., not
#' \emph{i} with \emph{j} then \emph{j} with \emph{i}). Some sanity checks are
#' made.
#' @param gerr (\code{find_POPs_lglk}) Genotyping error rate (apart from any
#' AA/AO-type errors)--- which had better be a small number. You have to pick
#' it yourself, but it is only used to "robustify" the lglk-based (PO)PLOD for
#' testing POPs vs UPs, and thus can be a rough guesstimate. FWIW we have used
#' 0.01 (i.e. 1\%) for some datasets, which is rather higher than replicate
#' analysis suggests, but is "safe" while still being small enough not to muck
#' up overall statistical performance. You should really do the same thing
#' yourself, and if you are very paranoid then try sensitivity analyses; but in
#' practice, the results of \code{find_POPs_lglk} are liable to be so clear-cut
#' that you may not feel it necessary to try more than one small value...
#' @param max_diff_loci (\code{find_duplicates}) max number of discrepant 4-way
#' genotypes to tolerate in "identical" fish. Only the pairs with fewer than
#' \code{max_diff_loci} discrepancies will be retained. Try increasing this
#' from say 10 upwards, and hopefully nothing much will change (though at some
#' point things will change a lot, as you get into the non-duplicate bit of the
#' distribution). See \bold{Duplicates} for how to remove duplicates from the
#' data.
#' @param limit_pairs Integer. Defines the \emph{maximum} number of candidate
#' pairs to keep. Will provide a warning if the number of identified pairs
#' equals limit_pairs.
#' @param nbins,minbin,maxbin \code{find_XXX} functions summarise their
#' pairwise comparison statistics into bins (in the part of the range where
#' exact values are uninteresting), as well as returning specific pairs that
#' pass the "interesting" threshold. \code{nbins} sets the number of bins,
#' \code{minbin} sets the top value of the lowest bin (so that bin stretches
#' from -Inf to \code{minbin} for HSPs); \code{maxbin} sets the highest. For
#' HSPs, the minimum is 3 bins (-Inf:minbin),[minbin:maxbin),[maxbin:Inf).
#' \code{minbin} is not used for duplicates or (at present) for POPs, since the
#' statistics there are defined so that the lowest possible value is 0. The
#' defaults for \code{minbin} and/or \code{maxbin} may not be what you need in
#' all cases, so be prepared to select manually and then re-run. For
#' duplicates, where calculations can be slow for big datasets, you can set
#' \code{nbins=0} to disable binning and focus instead on just finding the
#' pairs with fewer than \code{max_diff_loci} discrepancies. Each pairwise
#' calculation normally loops over all the loci, but is aborted when the
#' running total of discrepant loci reaches \code{maxbin} (or, if
#' \code{nbins=0}, when it reaches \code{max_diff_loci}), thus saving
#' considerable time. It is therefore not sensible to have
#' \code{maxbin<max_diff_loci} (think about it!).
#' @param keep_thresh (\code{find_HSPs} and \code{find_POPs}) is the analog of
#' \code{max_diff_loci} for \code{find_duplicates}. It determines which pairs
#' to retain for individual inspection. For \code{find_HSPs}, this is the
#' lowest retained PLOD; for \code{find_POPs} (currently), it's the highest
#' retained \code{wpsex}. Set it with the aim of including "anything
#' interesting" (ie not \emph{missing} any interesting pairs) and expect false
#' positives- ie be willing to have some weaker kin in there, and to
#' subsequently filter those out yourself "manually", as per vignette. For
#' HSPs, values like 0 (near the HTP mean) or -5 are a good start. For POPs,
#' experiment (at a total guess, try 0.1). You may have to re-run the function
#' a coupla times if you have been too brutal or too generous here- though "too
#' generous" can be fixed post hoc just by filtering the result, as long as you
#' haven't generated tooooo many pairs (see \code{limit_pairs}).
#' @param eta (\code{find_HSPs} and \code{find_POPs}) Not essential; limit for
#' calculating empirical mean and var PLOD, to compare with theoretical
#' \code{mean_UP} and \code{var_UP}. If you care about this (and you might not,
#' since for \code{find_HSPs} the observed/expected binwise comparison is
#' perhaps clearest), then set it to somewhere above 0 that should include
#' almost all UPs and exclude most strong kin; that's an \emph{upper} limit for
#' HSPs, and a \emph{lower} limit (currently) for POPs. The defaults are
#' guesses, and might not do that adequately, so an apparent mismatch may not
#' actually matter; be prepared to look at the histograms and think. The
#' general idea is that the number of UPs should dominate any other kin-type in
#' large sparsely-sampled datasets, so "contaminating" the empirical UP
#' statistics with a few weakish kin at the top end shouldn't mess up mean and
#' variance too much.
#' @param WPSEX_UP_POP_balance (\code{find_POPs}) loci receive a weight which
#' is proportional to (difference in probability of pseudo-exclusion between UP
#' and POP) / (variance of indicator of pseudo-exclusion). But, should this be
#' variance assuming UP or POP? \code{WPSEX_UP_POP_balance} sets the balance;
#' bigger values make it more UPpity, so placing more emphasis on avoiding
#' false-positives (which is probably the Right Thing To Do). 0.99 could be
#' completely fine... (but hopefully \code{WPSEX_UP_POP_balance} won't affect
#' the result much anyway.)
#' @param show_plot whether to plot log histogram. Regardless, plot will not be
#' shown if other arguments would lead to stupid result (e.g. no bins...).
#' @return A \code{data.frame} with extra attributes (see below) and at least 3
#' columns: statistic (\code{PLOD} or \code{wpsex} or \code{ndiff} number of
#' mismatching genotypes), then \code{i} and \code{j} (index in \code{subset1}
#' and \code{subset2} respectively of the first pair-member). Note that
#' \code{i} and \code{j} refer to the \emph{subsets}, not to the rows of the
#' original \code{snpg}- so if you are comparing subsets within the same
#' overall \code{snpg}, then you have to adjust accordingly. The attributes in
#' all cases include \code{bins} (upper boundaries), some kind of count
#' statistic for number of comparisons in each bin (names vary),
#' \code{binprobs} (theoretical CDF for UPs in \code{find_HSPs}; should exist
#' for \emph{POPs} (not UPs) in \code{find_POPs} (the \code{wpsex} version) but
#' currently doesn't), some of the input parameters, and the \code{call} that
#' invoked the function. \code{find_POPs} adds a column named \code{nABOO},
#' showing the number of AB/OO exclusions for that potential POP. This is a
#' useful additional diagnostic; it should be close to 0 for true POPs (it can
#' only result from genotyping error or mutation, whereas AAO/BBO can result
#' from nulls). For UPs, I was seeing values typically in the low 20s, which is
#' pretty good separation. \code{find_HSPs} and \code{find_POPs} have a bunch
#' of extra attributes which should be reeeeasonably clear. For
#' \code{find_HSPs}, \code{mean_sub_PLOD} and \code{var_sub_PLOD} are the
#' empirical means & var below \code{eta}, andy they should be close to
#' \code{mean_UP} and \code{var_UP} \emph{iff} \code{eta} has been chosen
#' sensibly. For \code{find_POPs}, the same goes for \code{mean_wpsex_hi} and
#' \code{var_wpsex_hi}. For duplicates, not \emph{all} pairwise duplicates are
#' recorded, unless the subsets are different- otherwise you could have
#' quadratic horror of enormous numbers of pairs arising from a cluster of say
#' 100 identical controls! Since "duplication" is transitive (ie if i & j are
#' the same, and i & k are the same, then j & k must also be the same), only
#' the necessary ones are recorded to allow you to filter out yourself
#' afterwards. e.g., if samples 1, 3, 5, and 6 are all duplicates, you'll get
#' this: \item{ # without "ndiff" column::}{} \item{ i j::}{} \item{ 3 1::}{}
#' \item{ 4 3::}{} \item{ 6 4:}{but you won't see the pairings for 1/4, 1/6,
#' 3/6. If you just want to strip out all duplicates bar one in each group (and
#' you don't care which one is kept), then you can use the function
#' \code{drop_dups_pairwise_equiv} - see \bold{Examples}. For POPs and HSPs,
#' the following are also returned as attributes (that can be accessed by
#' \code{@} if \code{atease} is loaded). The main point is that the "boring"
#' below-threshold pairs get put into bins and are not kept individually. The
#' names sometimes change depending on which statistic is being used.} \item{
#' eta}{false-positive cutoff to be applied to the statistic in question
#' (automatically done if \code{rough_n_pairs_to_keep==NA}, or up to you if
#' not). Variance of the stat will only be calculated from values to the "UP
#' side" of \code{eta}. However, the set of retained pairs/individuals is
#' actually controlled by...} \item{ keep_thresh}{the cutoff used to retain
#' "interesting" pairs. Usually obvious from the range of statistic values.}
#' \item{ mean_sub_<stat>, var_sub_<stat>}{empirical values for the statistic
#' when it is below \code{eta} (ie nearly always).} \item{
#' mean_theory,}{var_theory:of the statistic, to compare to previous.} \item{
#' n_<stat>_in_bin}{number of pairs whose statstic fell within the range of
#' each bin} \item{ bins}{cutpoints for the bins. These should be quantiles,
#' according to the SPA; so if practice matches theory, the numbers-per-bin
#' should all be similar.}
#' @section Kinformation: The idea is that kin-finding is based on a statistic
#' and a threshold \code{eta}, where the latter is chosen to keep
#' false-positives down to a user-specified level. Anything "beyond" \code{eta}
#' will be treated as a kin-pair ("beyond" depends on how the statistic is
#' defined, i.e. whether a kin-pair should come out very low or very high).
#' However, you're also likely to want to look post hoc at the distro of
#' computed statistics \emph{near} \code{eta}, to see whether separation is as
#' clean (or otherwise) as expected- and also very unbeyond \code{eta} into the
#' zone where UPs are entirely dominant, to check that theory is OK. So, as
#' well as returning the "interesting" pairs that have a statistic close to or
#' on the non-UP size of \code{eta}, the POP and HSP versions also return
#' \emph{summaries} of the distribution of the statistic. The thing is that
#' there will be zillions of statistics from UPs- enough to blow out computer
#' memory- and they are not individually interesting. Specifically, the main
#' things returned are:
#' 
#' \itemize{ \item mean and variance of stats. Computation is restricted to
#' those on the UP-side of \code{eta} (which is nearly all of them, usually) in
#' order to avoid distortion from non-UP cases. The latter will often be so
#' rare that distortion would be negligible- but means and variances are not
#' "robust", so . Almost all will be include \item counts of binned stats,
#' regardless of whether above or below \code{eta}. The bins are set based on
#' SPAs to the theoretical distributions, and chosen so that an equal number of
#' UP-pairs should fall into each bin. \item cases where the stat is
#' "interesting", i.e. on the non-UP side of \code{keep_thresh}, as a
#' \code{data.frame}. See \bold{Value} for details }
#' 
#' The process is controlled by three numbers: \code{nbins} for number of bins,
#' \code{eta} itself, and some nearby threshold \code{keep_thresh} on the
#' UP-side of \code{eta} (it will be automatically set to \code{eta} otherwise)
#' to determine which pairs are explicitly retained for your inspection. There
#' are two ways to specify \code{eta} and \code{keep_thresh}. Usually, you
#' would start with the indirect method, where you choose the
#' predicted-false-positive proportion of UP-pairs via the parameter
#' \code{one_in_X_eta}, and \code{rough_n_pairs_to_keep}. The routines then use
#' SPAs to the corresponding values of \code{eta} and \code{keep_thresh}; the
#' returned value of \code{eta} is what you can subsequently use to make the
#' actual kin-decisions yourself after the event (by subsetting the
#' "interesting" pairs, comparing the statistic for each pair to \code{eta})-
#' assuming that observed does match expected.
#' 
#' But, sometimes it doesn't. In that case, the predicted values of \code{eta}
#' and \code{keep_thresh} may be way off the mark, and lead to retaining faaar
#' too few or too many pairs. If so, then look at the histogram of retained
#' statistics from an initial run, and try setting \code{eta} and/or
#' \code{keep_thresh} manually, rather than futzing around with the indirect
#' parameters until you get what you were after.
#' @keywords misc
#' @examples
#' 
#' \dontrun{
#' ## And out-of-date!
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
#' #second_half <- (1:nrow( ckmini2)) \
#' #test1 <- find_duplicates( ckmini2, subset1=first_half, subset2=first_half, max=0)
#' #test1
#' ##  ndiff i j
#' ##1     0 3 1
#' #droppies1 <- first_half[ drop_dups_pairwise_equiv( test1[,2:3])] # NB must do lookup in subset
#' #test2 <- find_duplicates( ckmini2, subset1=second_half, subset2=second_half, max=0)
#' #droppies2 <- second_half[ drop_dups_pairwise_equiv( test2[,2:3])] # 4
#' ## Now check 2nd half vs 1st
#' #test2_1 <- find_duplicates( ckmini2,
#' #    subset1=first_half \
#' #    subset2=second_half \
#' #    max=0)
#' ## Simpler since no internal checks. Just remove 2nd-halfers that match something in the 1st-half
#' #droppies2_1 <- (second_half \
#' #droppies <- c( droppies1, droppies2, droppies2_1)
#' #ckmini2_nodups2 <- ckmini2[ -droppies,]
#' ## HSPs: comparing everything with itself (not sensible for real data, should take out adults first)
#' ## set threshold for 1 FP
#' #test <- find_HSPs( ckdata, one_in_X_eta=sqr( nrow( ckdata))/2 )
#' ## POPs: Ad-Ju comps; again 1 FP
#' #test <- find_POPs( ckdata, subset1=adults, subset2=juves,
#' #    one_in_X_eta=length( adults) * length( juves), rough_n_pairs_to_keep=500)
#' }
#' 
#' @export find_duplicates
"find_duplicates" <-
function(
  snpg,
  subset1= 1 %upto% nrow( snpg),
  subset2= subset1,
  max_diff_loci,
  limit_pairs= 0.5*nrow(snpg),
  nbins= 50,
  maxbin= ncol( snpg)/2,
  show_plot = TRUE
){
  # Sanity...
stopifnot(
    is.numeric( subset1),
    is.numeric( subset2),
    all( !duplicated( subset1)),
    all( !duplicated( subset2)),
    my.all.equal( subset1, subset2) || !length( intersect( subset1, subset2)),
    !missing( max_diff_loci),
    max_diff_loci <= maxbin,
    nbins >= 0 # 0 is OK; just return actual pairs below max_diff_loci
  )

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

  # MVB: I don't trust SB's code here because rounding error may mean that the number of bins would not match nbins
  # SB had: bins <- seq(minbin+binterval, (minbin + (nbins*binterval)), binterval)
  # nbins=0 is OK cos it saves heaps of time, but obvs you don't get binned results
  binterval <- maxbin / max( nbins-1, 1) # bin *starting* at maxbin counts all bigger
  bins <- seq( from= 0, to= maxbin, length= max( 1, nbins-1))[-1] # avoid woe when nbins=0
  # binprobs not set

  # Trying special-cases here to minimize copying
  if( my.all.equal( subset1, subset2)) {
    if( !my.all.equal( subset1, 1 %upto% ncol( temp_snpg))) {
      temp_snpg <- temp_snpg[, subset1]
    }

    result <- DUP_paircomps_lots(
        geno1= temp_snpg,
        geno2= temp_snpg,
        symmo= TRUE,
        max_diff_loci= max_diff_loci,
        keep_n= limit_pairs,
        nbins= nbins,
        binterval= binterval,
        maxbin = maxbin
      )
  } else { # different subsets
    result <- DUP_paircomps_lots(
        geno1= temp_snpg[ , subset1],
        geno2= temp_snpg[ , subset2],
        symmo= FALSE,
        max_diff_loci = max_diff_loci,
        keep_n = limit_pairs,
        nbins = nbins,
        binterval = binterval,
        maxbin= maxbin
      )
  }

  n_ndiff_in_bin <- result$n_ndiff_in_bin

  # just return the data.frame with 3 columns, everything else goes in
  # the attributes
  result <- with( result, data.frame( ndiff=big_similar, i=big_i, j=big_j))
  result@call <- sys.call()
  result@bins <- bins
  result@n_ndiff_in_bin <- n_ndiff_in_bin

  # warning if we're running up against storage constraints
  if(length(result$ndiff) == limit_pairs){
    message("Returning the ", limit_pairs, " most similar pairs, increase limit_pairs if more are required")
## warning("Number of returned duplicates equals limit_pairs. There may be more than limit_pairs duplicates. Increase limit_pairs to make sure you have them all!")
  }

  if( show_plot && (length( n_ndiff_in_bin) > 2)) {
    loggo <- log( n_ndiff_in_bin[1:(length( n_ndiff_in_bin)-2)] + 0.001,
        base=10)
    if( all( loggo < 0)){
      warning( "No pair below 'maxbin'--- nothing to plot")
    } else {
      plot(c(result@bins), loggo,
          ylim= c( 0, max( loggo) + 1),
          type = "S",
          xlab = "n different genos", ylab = "log10(#pairs)")
    }
  }

return( result)
}




#' Dup-finding with some missing genotypes
#' 
#' This documentation seems to have disappeared... though I \emph{think} this
#' is the function Paige relies on in GT for detecting recaptures! So it's
#' quite important... The basic point is that GT data \emph{does} sometimes
#' have missing loci, whereas data for CKMR kinference is not allowed to. All
#' pairwise comparisons will be done willy-nilly (ie regardless of the precise
#' missingness pattern for the pair), so subset your data beforehand to remove
#' Bad Eggs.
#' 
#' 
#' @param snpg a \code{snpgeno} object
#' @param subset1,subset2 numeric vectors of which samples to use (not logical,
#' not negative). Defaults to all of them. Iff \code{subset1} and
#' \code{subset2} are identical, only half the comparisons are done (i.e., not
#' \emph{i} with \emph{j} then \emph{j} with \emph{i}). Some sanity checks are
#' made.
#' @param max_diff_ppn What \emph{proportion} of non-missing (ie
#' scored-in-both) loci to treat as the threshold for duplicity?
#' @param limit if you hit this many "duplicates", it will stop, to avoid
#' blowing out memory. It means you set \code{max_diff_ppn} too high. For
#' consistency, we should probably have called this \code{keep_n} as per other
#' \code{find_XXX} functions.
#' @return A dataframe with columns \code{ppn_diff}, then \code{i} and \code{j}
#' (index in \code{subset1} and \code{subset2} respectively of the first
#' pair-member), then \code{ndiff} and \code{ncomp} (numbers of loci differing
#' and actually compared for that pair). Note that \code{i} and \code{j} refer
#' to the \emph{subsets}, not to the rows of the original \code{snpg}- so if
#' you are comparing subsets within the same overall \code{snpg}, then you have
#' to adjust accordingly.
#' @keywords misc
#' @export find_dups_with_missing
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

  # just return the data.frame with 3 columns, everything else goes into atts
  result <- with( result, data.frame( ppn_diff= big_ndiff / big_ncomp,  i=big_i, j=big_j, ndiff=big_ndiff, ncomp=big_ncomp))
  result@call <- sys.call()

return( result)
}


"find_HSPs" <-
function(
    snpg,
    subset1=1 %upto% nrow( snpg),
    subset2=subset1,
    limit_pairs= 0.5 * nrow( snpg),
    keep_thresh,
    eta= NULL,
    nbins= 50,
    minbin= NULL,
    maxbin= NULL
){
## snpg should have been thru 'prepare_PLOD_SPA' so it has @PPS
stopifnot( 'Kenv' %in% names( attributes( snpg)))

  # Sanity...
stopifnot( is.numeric( subset1) && is.numeric( subset2))
stopifnot( all( !duplicated( subset1)) && all( !duplicated( subset2)))
stopifnot( my.all.equal( subset1, subset2) || !length( intersect( subset1, subset2)))

  # SB has:
  hspPower_change <- snpg@hspPower_checksum != calc_hspPower_checksum( snpg$locinfo)
  # MVB: though that's not the best sum to check! pbonzer should sum to 1 per row... ditto PUP<X> below
  PLODSPA_change <- snpg@PLODSPA_checksum != calc_PLODSPA_checksum( snpg$locinfo)

  if(hspPower_change | PLODSPA_change) {
      warning("snpg$locinfo appears to have been modified after hsp_power() and/or prepare_PLOD_SPA() were last called. I sure hope you know what you're doing...")
  }

  define_genotypes()
  og <- options( vecless.print=FALSE)
  on.exit( options( og))

  # For 4way loci, temporarily treat XO as XX...
  # ... have already adjusted the LOD entries so that new_LOD6( XX/..) <- LOD4( XXO/..)
  # ... use the LOD that's in Kenv, where SPA is calculated

  ### Do we need LOD6/4/3 ?
  # extract.named( snpg@locinfo[ cq( useN, LOD6, LOD4, LOD3)])
  useN <- snpg@locinfo$useN
  temp_snpg <- snpg
  recode4to6temp <- function( x) { x[ x=='AO'] <- AA; x[ x=='BO'] <- BB; x}
  recode3to6temp <- function( x) { x[ x=='AO'] <- AA; x[ x=='BO'] <- BB; x[ x=='OO'] <- BB; x}

  temp_snpg[ , useN == 4] <- recode4to6temp( snpg[, useN == 4])
  temp_snpg[ , useN == 3] <- recode3to6temp( snpg[, useN == 3])
  temp_LOD <- snpg@Kenv$LOD # already done in prepare_PLOD_SPA, and based on useN to switch 6/4/3

  # Remove extranea
  attributes( temp_snpg) <- attributes( temp_snpg)[ 'dim']
  temp_snpg <- t( temp_snpg)

  li <- snpg@locinfo # convenient

  # MVB: don't trust length of seqs with non-integer steps!
  # bins <- seq(minbin+binterval, (minbin + (nbins*binterval)), binterval)

  EHSP <- sum( li$E.HSP) / 3 # who knows?!
  if( is.null( eta)) {
    eta <- EHSP / 3 # who knows?!
  }

  if( is.null( minbin)) {
    VUP <- sum( li$V.UP)
    EUP <- sum( li$E.UP)
    minbin <- EUP - 4 * sqrt( VUP) # 0.003% of UPs below that
  }

  if( is.null( maxbin)) {
    EPOP <- sum( li$E.POP)
    maxbin <- EPOP + (EPOP - EHSP) * 0.1
  }
  bins <- seq( from=minbin, to=maxbin, length=nbins)
  binterval <- bins[2] - bins[1]

  # Distro of PLOD|UP via SPA
  binprobs <- c( snpg@Kenv$CDF( bins), 1)

  # Trying special-case "all vs all" here to minimize copying
  if( all(subset2 %in% subset1) & all(subset1 %in% subset2)) {
    if( !my.all.equal( subset1, 1 %upto% ncol( temp_snpg))) { # NB transpose!
      temp_snpg <- temp_snpg[, subset1]
    }

    xresult <- HSP_paircomps_lots(
        pair_geno = temp_LOD@mg,
        LOD = t(temp_LOD),
        geno1 = temp_snpg,
        geno2 = temp_snpg,
        symmo = TRUE,
        eta = eta,
        min_keep_PLOD = keep_thresh,
        minbin = minbin,
        binterval = binterval,
        nbins = nbins,
        keep_n = limit_pairs)
  } else { # different subsets
    xresult <- HSP_paircomps_lots(
        pair_geno= temp_LOD@mg,
        LOD= t( temp_LOD),
        geno1= temp_snpg[ , subset1],
        geno2= temp_snpg[ , subset2],
        symmo= FALSE,
        eta= eta,
        min_keep_PLOD= keep_thresh,
        minbin = minbin,
        binterval = binterval,
        nbins = nbins,
        keep_n = limit_pairs
      )
  }

  # warning if we're running up against storage constraints
  if( !missing(keep_thresh) && (length(xresult$big_PLOD) == limit_pairs)){
    message( "Returning just the ", limit_pairs, " pairs with the highest PLOD scores; increase 'limit_pairs' if more are required")
  }  ## This warning will probably need to be modified in find_POPs, find_duplicates, etc. as well.

  result <- with( xresult, data.frame( PLOD=big_PLOD, i=big_i, j=big_j))
  attributes( result) <- c( attributes( result),
      xresult %without.name% cq( big_PLOD, big_i, big_j))

  # assign extra info as attributes
  result@mean_UP <- snpg@Kenv$dK( 0)   # was called mean_theory
  result@var_UP <- snpg@Kenv$ddK( 0)
  result@mean_HSP <- snpg@Kenv$dK( 0) + sum(snpg@locinfo$Ediff) ##
  ##   result@mean_HSP <- sum(snpg@locinfo$E.HSP) ## this should be OK
  ##   result@mean_UP <- sum(snpg@locinfo$E.UP) ## this should be OK too
  result@mean_POP <- sum(snpg@locinfo$E.POP)
  result@mean_FSP <- sum(snpg@locinfo$E.FSP)
  attributes( result) <- c( attributes( result), returnList( bins, binprobs, eta, keep_thresh))

  result@call <- sys.call()

return( result)
}


"find_POPs" <-
function( snpg, subset1=1 %upto% nrow( snpg), subset2=subset1,
    limit_pairs=0.5*nrow(snpg),
    keep_thresh,
    eta= NULL,
    nbins= 50,
    maxbin = NULL,
    WPSEX_UP_POP_balance=0.99) {
###################
  # Sanity...
define_genotypes()
stopifnot(
    is.numeric( subset1),
    is.numeric( subset2),
    all( !duplicated( subset1)),
    all( !duplicated( subset2)),
    my.all.equal( subset1, subset2) || !length( intersect( subset1, subset2)),
    snpg %is.a% 'snpgeno',
    my.all.equal( snpg@diplos, genotypes6) || my.all.equal( snpg@diplos, genotypes4_ambig)
  )

  ## # SB checks--- but NB rowSums( pbonzer==1) !
  ## hspPower_change <- snpg@hspPower_checksum != with(snpg@locinfo, sum(pbonzer, snerr, useN))
  ## PLODSPA_change <- snpg@PLODSPA_checksum != with(snpg@locinfo, sum(useN, LOD6, LOD4, LOD3,
  ##                                                 PUP6, PUP4, PUP3))
  ## if(hspPower_change | PLODSPA_change) {
  ##   warning("snpg$locinfo appears to have been modified after hsp_power and/or prepare_PLOD_SPA were last called. I sure hope you know what you're doing...")
  ## }

  # Decide based #apparent exclusions of AA/BB form, using 4way genos, though
  # ... it's really AAO/BBO so not a true exclu but
  # ... close among pop-loci
  # Sticking with 4way genos so that genotyping errors are low

  # Instead of Nexclu, uses a wted sum of "exclus" in 4way genos to max expected diff between POP and UP
  # This version doesn't allow for geno errors, but does realize that AO/BO could happen
  # Doesn't bother with OO-AB

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

  # Algo should work either with 4way or 6way genotypes (no C allele though)
  # C code gets told, via AAish and BBish, which codes are really AAO and BBO
  # For 6way, remap AO->AA, BO->BB, and tell C just to check for AA and BB
  temp_snpg <- snpg[ , pop_loci]
  if( my.all.equal( snpg@diplos, genotypes6)) {
    recode_6_as_pseudo4 <- function( x) { x[ x=='AO'] <- AA; x[ x=='BO'] <- BB; x}
    temp_snpg <- recode_6_as_pseudo4( temp_snpg)
    AAish <- match( 'AA', snpg@diplos)
    BBish <- match( 'BB', snpg@diplos)
  } else { # 4way
    AAish <- match( 'AAO', snpg@diplos)
    BBish <- match( 'BBO', snpg@diplos)
  } # else NYI, which would have been picked up in sanity checks

  # Remove extranea
  attributes( temp_snpg) <- attributes( temp_snpg)[ 'dim']
  temp_snpg <- t( temp_snpg)

  n_loci <- ncol( snpg)

  # Prepare for diagnostics of #excl
    ## Defaults to two-sigma above *UP* mean... which will cover almost all

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

    K <- compile_vecless( K(0))
    dK <- compile_vecless( dK(0)) ## dK(0) should be the mean of wpsex for true UPs
    ddK <- compile_vecless( ddK(0)) ## should be the var of ditto

    #  n_sim_check <- 1000
    #  Ktest <- function( tt) {
    #    x <- matrix( runif( n_sim_check * n_loci) < pex_up, n_loci, n_sim_check)
    #    ewx <- exp( x*ww*tt)
    #    colSums( log( ewx))
    #  }
    #

    # Could do now predict binprobs *for UPs) via renorm_SPA_cumul
    # Let's not
    # What we need, is...
  if( is.null( maxbin)) {
      maxbin  <-  dK(0)+2*sqrt(ddK(0))
  }

  if(is.null(eta) ) {
      eta <- dK(0)-3*sqrt(ddK(0))
  }

  if( nbins<2) {
    warning( 'nbins<2 is senseless; gonna use 2')
    nbins <- 2
  }
  bins <- seq( from= 0, to= maxbin, length=nbins)
  ## binterval <- maxbin / (nbins-1) # bin *starting* at maxbin counts all bigger
                                        # binprobs not set
  binterval <- bins[2] - bins[1]

  # Trying special-cases here to minimize copying
  symmo <- my.all.equal( subset1, subset2)
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
        keep_n = limit_pairs,
        AAO= AAish,
        BBO= BBish,
        nbins = nbins,
        binterval = binterval
      )
  } else { # different subsets
    result <- POP_wt_paircomps_lots(
        geno1= temp_snpg[ ,subset1],
        geno2= temp_snpg[ ,subset2],
        w= ww,
        symmo= FALSE,
        eta= eta,
        max_keep_wpsex= keep_thresh,
        keep_n = limit_pairs,
        AAO= AAish,
        BBO= BBish,
        nbins = nbins,
        binterval = binterval
      )
  }

  # warning if we're running up against storage constraints
  if(length(result$big_wpsex) == limit_pairs){
    message("Returning just the ", limit_pairs, " pairs with the most POP-like wpsex scores, increase 'limit_pairs' if more are required")
  }

  n_wpsex_in_bin <- result$n_wpsex_in_bin  ## SB tweak; the 'below zero' bin must be empty
  # bins <- seq(0+binterval, 0+(nbins*binterval), binterval)
  ## SB addition. Bloody zero-base.

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
  attributes( result) <- c( attributes( result), returnList(
      bins, n_wpsex_in_bin, eta, keep_thresh))

  result@n_loci <- length( pop_loci)
  result@mean_UP <- dK( 0)
  result@var_UP <- ddK( 0)
  result@call <- sys.call()

return( result)
}


"find_POPs_lglk" <-
function(
    snpg,
    subset1= 1 %upto% nrow( snpg),
    subset2= subset1,
    gerr,
    limit_pairs= 0.5*nrow(snpg),
    keep_thresh,
    eta= NULL,
    nbins= 50,
    minbin= NULL,
    maxbin= NULL
    # WPSEX_UP_POP_balance=0.99
){
###################
  define_genotypes()

  # Sanity...
stopifnot(
    !missing( gerr),
    is.numeric( subset1), is.numeric( subset2),
    all( !duplicated( subset1)), all( !duplicated( subset2)),
    my.all.equal( subset1, subset2) || !length( intersect( subset1, subset2)),
    snpg %is.a% 'snpgeno',
    my.all.equal( snpg@diplos, genotypes6) || my.all.equal( snpg@diplos, genotypes4_ambig)
  )

  # A POP is "just like" an HSP where the ppn coin is 1.0 not 0.5
  # but we add a little bit of geno error for statistical lubrication
  # Call 'hsp_power' to get LODs manually for this case
  # then 'prepare_PLOD_SPA' to handle useN stuff
  # then use HSP machinery.

  snpg <- hsp_power( snpg, k= 1-gerr) # !! Captain Sneakypants !!
  snpg <- prepare_PLOD_SPA( snpg, n_pts_SPA_renorm=201)

  mc <- match.call()
  mc$gerr <- NULL
  mc[[1]] <- quote( find_HSPs)
  # Kinda "delayed assign" for snpg, to look here; avoids putting whole object into result@call
  mc$snpg <- substitute( sys.frame( n)$snpg, list( n=sys.nframe()))
  result <- eval.parent( mc)

  atts <- attributes( result) %without.name% 'call'
  atts <- atts[ names( atts) %that.dont.match% 'FSP']
  names( atts) <- sub( 'HSP', 'POP', names( atts))
  attributes( result) <- atts

  # assign extra info as attributes
  result@mean_UP <- snpg@Kenv$dK( 0)   # was called mean_theory
  result@var_UP <- snpg@Kenv$ddK( 0)
  result@mean_POP <- snpg@Kenv$dK( 0) + sum(snpg@locinfo$Ediff)
  result@call <- sys.call()

return( result)
}




#' Plot for splitting FSPs from POPs
#' 
#' Plots an absolute-frequency histogram for the output of
#' \code{find_FSPs_from_POPs_v2()}.
#' 
#' 
#' @param fsps2 the output of a call to \code{find_FSPs_from_POPs_v2()}
#' @param bin hist bin width. Used to define \code{breaks} (along with
#' \code{xlim}, if given), so you can't manually pass in \code{breaks}.
#' @param FSPmean plot the expected mean FP stat for FSPs? Default TRUE
#' @param POPmean plot the expected mean FP stat for POPs? Default TRUE
#' @param ... additional pars, passed to \code{hist()} (?should we force naming
#' them in a list?)
#' @seealso PLOD_loghisto
#' @keywords misc
#' @export FSP_POP_histo
"FSP_POP_histo" <-
function(fsps2, bin = 100, FSPmean = TRUE, POPmean = TRUE, main = "", ...) {

        palette(c("#0D0887FF", "#48039FFF", "#7401A8FF", "#9D189DFF", "#BF3984FF",
                  "#DA596AFF", "#EE7B51FF", "#FBA238FF", "#FCCE25FF"))

        if(exists("xlim")) {
            lb <- min(xlim)
            ub <- max(xlim)
        } else {
            lb <- min(fsps2$FPstat)
            ub <- max(fsps2$FPstat)+bin
        }

        hist.plod=hist(fsps2$FPstat, breaks=seq(lb, ub, bin),
                       col="lightgrey",xlab="FPstat", ...)
        if( POPmean) { abline(v = fsps2@E_FPstat[1], col = 1, lwd = 2) }
        if( FSPmean) { abline(v = fsps2@E_FPstat[2], col = 9, lwd = 2) }

        legendBits <- data.frame(allNames = c("POP","FSP"), allNumbers = c(1,9))
        if(!POPmean) {
            legendBits <- legendBits[legendBits$allNames != "POP",]
        }
        if(!FSPmean) {
            legendBits <- legendBits[legendBits$allNames != "FSP",]
        }

        legend("topright", legend = legendBits$allNames,
               lwd = 2, lty = 1, col = legendBits$allNumbers, bg = "white")
    }


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




#' QC checks on
#' 
#' This test looks at whether the allele frequencies in a given fish seem
#' right, or if there are discrepancies due to (i) degraded DNA or (ii) sample
#' contamination. Useful both for finding outlier samples, and for checking
#' whether the loci are collectively working as they should (and as is assumed
#' by all the calculations in \code{kinference}). The histogram should coincide
#' nicely with its predicted line.
#' 
#' 
#' @param snpg a \code{snpgeno} object
#' @param target which weighting should be used. \code{"rich"} is meant to be
#' more sensitive for detecting contaminated data (too many heterozygotes) and
#' \code{"poor"} for detecting DNA degredation (too few). In practice, there is
#' often little difference.
#' @param hist_pars list of parameters to pass to \code{hist}. If you are very
#' sneaky, you can pass in an \code{expression} to be evaluated inline instead
#' (ie fly-hacking).
#' @param showPlot show the plot? Default TRUE
#' @keywords misc
#' @export hetzminoo_fancy
"hetzminoo_fancy" <-
function( snpg, target=c( 'rich', 'poor'), hist_pars=list(), showPlot = TRUE) {
###################
  define_genotypes()
  extract.named( snpg@locinfo[ cq(  PUP4, pbonzer)]) ## used to also call use6; not used
  p0 <- pbonzer[,'O'] + pbonzer[,'C']
  pA <- pbonzer[,'A']
  pB <- pbonzer[,'B']

  v <- 2*pA*pB + sqr( p0) - sqr( 2*pA*pB-sqr( p0)) ## Pr(AB) + Pr(OO) - (Pr(AB) - Pr(OO))^2
  target <- match.arg( target)
  edash <- if( target=='rich') {
     (2*pA*p0+sqr(pA)) *   ## Pr(AA|AO)
     (1-sqr(1-pB)) +       ## original; Pr(BB|AB|BO)
     (2*pB*p0+sqr( pB)) *  ## Pr(BO|BB)
     (1-sqr(1-pA)) +       ## modified; Pr(AA|AB|AO)
     sqr( p0) *            ## Pr(OO)
     (1-sqr( p0))          ## Pr(!OO)
    } else {
      2*(pA*pB + pA*p0 + pB*p0) # Equals Pr(heterozygote) under 3-way HWE
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

  # MVB Dec 2020: I don't this recoding is required cos we only need AB and OO
  # ... anyway, version below won't work for genotypes4_ambig
  # recode4to6temp <- function( x) { x[ x=='AO'] <- AA; x[ x=='BO'] <- BB; x}
  # temp_snpg <- recode4to6temp( temp_snpg) # (AA,AO) -> AA; (BB,BO) -> BB

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

  dens_SPA <- renorm_SPA( K, dK, ddK, return_what='func', already_vectorized=TRUE)

  # optional graphics and/or user-specified outputs
  switch( mode( hist_pars),
    list = {
        hist_pars <- add_list_defaults( hist_pars,
            main=target,
            xlim= range( whmo), # so cutoff lines show
            xlab='', nclass=50)
        if (showPlot) {
            lv <- do.call( 'hist', c( list( x=whmo), hist_pars))
            with( lv, lines( mids, diff( breaks) * dens_SPA( mids) * sum( counts), col='green'))
        }
        # abline( v=ncuts, col='red')
      },
    expression = eval( hist_pars),
    NULL = NULL
  )

return( c( whmo))
}




#' PLOD histogram
#' 
#' Plots an absolute-frequency histogram for the output of
#' \code{\link{find_HSPs}}, with the lower bound set by the user. Lower bounds
#' should be set to exclude (as much as possible) the UP bump, as this will
#' otherwise swamp the signal from the HSP bump. Users must manually set a
#' lower bound for full-sibling PLODs (\code{fullsib_cut}) on order to exclude
#' full-siblings from the variance estimate for HSP PLODs.
#' 
#' 
#' @param hsps the output of a call to \code{find_HSPs()}
#' @param lb PLOD lower bound for plot extent. Should exclude the UP bump
#' @param ub PLOD upper bound for plot extent. Defaults to maximum PLOD score
#' plus a little padding
#' @param fullsib_cut PLOD score above which there are only full-sibs
#' @param bin hist bin width. Default 5. lb, ub, and bin together define
#' \code{breaks}, so you can't pass \code{breaks} via \code{...}
#' @param HSPmean plot the mean PLOD for HSPs? Default TRUE
#' @param HSPdist plot the distribution of PLOD for HSPs? Default TRUE
#' @param POPmean plot the mean PLOD for POPs? Default TRUE
#' @param FSPmean plot the mean PLOD for FSPs? Default TRUE
#' @param main graph title, passed straight to \code{hist()}
#' @param ... additional pars, passed to \code{hist()}
#' @seealso PLOD_loghisto
#' @keywords misc
#' @export HSP_histo
"HSP_histo" <-
function(hsps, lb, ub = max(hsps$PLOD)+bin, fullsib_cut, bin = 5,
             HSPmean = TRUE, HSPdist = TRUE, POPmean = TRUE, FSPmean = TRUE, main = "", ...) {

        palette(c("#0D0887FF", "#48039FFF", "#7401A8FF", "#9D189DFF", "#BF3984FF",
                  "#DA596AFF", "#EE7B51FF", "#FBA238FF", "#FCCE25FF"))

        hist.plod=hist(hsps$PLOD[hsps$PLOD > lb & hsps$PLOD < ub],breaks=seq(lb, ub, bin),
                       col="lightgrey",xlab="PLOD", main = main, ...)
        if( HSPmean) {
            E.hsp = hsps@mean_HSP
            abline(v=E.hsp,lwd=2,col=8)
        }
        if( POPmean) { abline(v = hsps@mean_POP, col = 1, lwd = 2) }
        if( FSPmean) { abline(v = hsps@mean_FSP, col = 9, lwd = 2) }
        if( HSPdist) {
            V.hsp=mean(sqr(hsps$PLOD[hsps$PLOD>E.hsp & hsps$PLOD < fullsib_cut]-E.hsp))
            obs.num <- hist.plod$counts
            exp.num <- 2*sum(hsps$PLOD>E.hsp & hsps$PLOD<fullsib_cut)*
                (pnorm(hist.plod$breaks[-1],E.hsp,sqrt(V.hsp))-
                 pnorm(hist.plod$breaks[-length(hist.plod$breaks)],E.hsp,sqrt(V.hsp)))
            points(hist.plod$mids,exp.num,pch=16,col=8,type='b')
        }
        legendBits <- data.frame(allNames = c("POP","HSP","FSP"), allNumbers = c(1,8,9))
        if(!POPmean) {
            legendBits <- legendBits[legendBits$allNames != "POP",]
        }
        if(!HSPmean) {
            legendBits <- legendBits[legendBits$allNames != "HSP",]
        }
        if(!FSPmean) {
            legendBits <- legendBits[legendBits$allNames != "FSP",]
        }

        legend("topright", legend = legendBits$allNames,
               lwd = 2, lty = 1, col = legendBits$allNumbers, bg = "white")
    }




#' Oddness metrics over an HSP histo
#' 
#' Plots an absolute-frequency histogram for the output of
#' \code{\link{find_HSPs}} (qv), with additional coloured regions showing
#' individuals with odd-looking CLOD scores (see \code{check_FPosity()} ),
#' ilglk stat (see \code{\link{ilglk_geno}}), or hetz stat (see
#' \code{\link{hetzminoo_fancy}}). You should choose a PLOD range for plotting
#' that excludes the UP bump, since UPs can swamp the signal from HSPs and
#' HTPs. \code{HSP_oddness_oneway} shows cases where at least one member of the
#' pair has an unusual score in the oddness metric, whereas
#' \code{HSP_oddness_twoway} shows cases where both members do.
#' 
#' 
#' @aliases HSP_oddness_oneway HSP_oddness_twoway
#' @param hsps the output of a call to \code{find_HSPs()}
#' @param snpg the \code{snpgeno} object from which \code{hsps} was built
#' @param lb PLOD lower bound for plot extent. Should exclude the UP bump
#' @param ub PLOD upper bound for plot extent. Defaults to the maximum PLOD
#' score plus a little padding
#' @param bin hist bin width. Default 5. lb, ub, and bin together define
#' \code{breaks}, so you can't pass \code{breaks} via \code{...}
#' @param CLOD_prop the quantile of CLOD below which animals are highlighted.
#' Default 0.001
#' @param ilglk_prop the quantile of ilglk stat below which animals are
#' highlighted. Default 0.001
#' @param hetz_prop the quantile of hetz stat below which animals are
#' highlighted. Default 0.001
#' @param ... additional pars, passed to \code{hist()}
#' @seealso PLOD_oddness_oneway
#' @keywords misc
#' @export HSP_oddness_oneway
"HSP_oddness_oneway" <-
function(hsps, snpg, lb, ub = max(hsps$PLOD)+10, bin = 5,
             CLOD_prop = 0.001, ilglk_prop = 0.001, hetz_prop = 0.001,
             ...) {

        palette("default")

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


"HSP_oddness_twoway" <-
function(hsps, snpg, lb, ub = max(hsps$PLOD)+10, bin = 5,
             CLOD_prop = 0.001, ilglk_prop = 0.001, hetz_prop = 0.001,
             ...) {

        palette("default")

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


"HSP_paircomps_lots" <-
function(pair_geno, LOD, geno1, geno2, symmo, eta, min_keep_PLOD, keep_n, minbin, binterval, nbins) {
    .Call(`_kinference_HSP_paircomps_lots`, pair_geno, LOD, geno1, geno2, symmo, eta, min_keep_PLOD, keep_n, minbin, binterval, nbins)
}




#' Locus selection for kin-finding
#' 
#' \code{hsp_power} can be used to predict how well a set of loci will work for
#' HSP-finding, and to prepare for some QC and kinference steps on serious
#' data. It returns the input \code{snpgeno} object with extra columns added to
#' the \code{locinfo} attribute, related to the per-locus mean and variance of
#' LOD (presumably an HSP/UP PLOD, though not inevitably) for different true
#' kin-types. It respects the per-locus decision about how precisely to
#' genotype (\code{useN=6/4/3}).
#' 
#' \describe{ \item{E.UP, V.UP}{mean & variance for UPs} \item{E.HSP,
#' E.POP,E.FSP}{as you would expect} \item{Ediff}{E.HSP - E.POP ie the
#' "absolute" power of that locus} \item{sdiff}{(E.HSP-E.POP)/sqrt(V.UP) which
#' is arguably better than \code{Ediff} for ranking loci} }
#' 
#' It also attaches \code{LOD}, \code{PUP}, and \code{ev01} elements (each a
#' matrix) to the \code{locinfo}. They have been made dull (see
#' \code{make_dull}) to improve your viewing experience, but they work fine for
#' all normal purposes (and you can always \code{unclass} them to remove the S3
#' class \code{dull}).
#' 
#' @param lociar \code{snpgeno} objects with the necessary ingredients
#' @param want_LOD_table can't think why you'd set this to FALSE
#' @param k target average kinship for LOD; 0.5 for HSPs, 0.25 for HTPs, etc.
#' @return \code{snpgeno} object with augmented columns in its "locinfo" attr.
#' @section Notes: \code{hsp_power} (and downstream) should get a refactor.
#' It's daft to store LODs for only one specific kin; it'd be better to always
#' calculate P1share and P2share as well as P0share (which is PUP), and then
#' compute whatever-is-needed later on-the-fly. As-is, we are re-computing P1
#' and P2 based on LOD and PUP OTF instead (which is also unsafe, because LOD
#' could have been calculated with k != 0.5).
#' @keywords misc
#' @examples
#' 
#' ## Need some examples!
#' # See vignette for now
#' 
#' @export hsp_power
"hsp_power" <-
function( lociar,
    want_LOD_table=TRUE, # T/F
    k # 0.5 for HSPs
){
############
  define_genotypes()
  li <- lociar$locinfo
  li1 <- li[1,]

`%without.names%` <- function( x, what) {
    new.names <- names( x) %except% what
    if( identical( new.names, names( x))) {
      return( x)   # also works if names(x) is NULL!
    }

    oatts <- attributes( x)
    # oatts must exist, since nameless-x returns earlier
    x <- x[ new.names]
    oatts$names <- new.names
    attributes( x) <- oatts
    return( x)
}

  temp0 <- with( li1, calc_g6probs_IBD0_scalar( pbonzer, snerr, record=TRUE))
  cg6p0 <- make_playback( calc_g6probs_IBD0_scalar, temp0)

  temp1 <- with( li1, calc_g6probs_IBD1_scalar( pbonzer, snerr, record=TRUE))
  cg6p1 <- make_playback( calc_g6probs_IBD1_scalar, temp1)

  temp2 <- with( li1, calc_g6probs_IBD2_scalar( pbonzer, snerr, record=TRUE))
  cg6p2 <- make_playback( calc_g6probs_IBD2_scalar, temp2)

  g6p0 <- with( li, cg6p0( pbonzer, snerr))
  g6p1 <- with( li, cg6p1( pbonzer, snerr))
  g6p2 <- with( li, cg6p2( pbonzer, snerr))

  s6 <- predict_hsp_util( g6p0, g6p1, g6p2, want_LOD_table, k=k)

  # For the 4-ways, must condense g6p's
  if( exists( 'genotypes4_ambig', inherits=FALSE)) { # TRUE unless overridden sneakily...
    extract.named( map6to4( g6p0, g6p1, g6p2))
    s4 <- predict_hsp_util( g4p0, g4p1, g4p2, want_LOD_table, k=k) %without.names% "matto"

    ### 3way too:
    extract.named( map6to3( g6p0, g6p1, g6p2))
    s3 <- predict_hsp_util( g3p0, g3p1, g3p2, want_LOD_table, k=k) %without.names% "matto"

    if( want_LOD_table) { # overrides predict_hsp_util's version of want_LOD_table
      # We want LOD6, PUP4, etc (matrices with cols "AB/AA" etc)
      # and ev01_6, ev01_4, etc (matrices with cols as per 'things' next)

      things <- cq( e0, e1, v0, v1)

      for( usy in cq( 3, 4, 6)){ # NB character!
        s <- get( 's' %&% usy) # s6 etc
        ev01 <- s[ things]
        names( ev01) <- things
        li[[ 'ev01_' %&% usy]] <- s$ev01 <- do.call( 'cbind', ev01)
        # eg li$ev01_4 will be a 4-col matrix

        li[[ 'LOD' %&% usy]] <- s@LOD # matrix
        li[[ 'PUP' %&% usy]] <- s@PUP # matrix
         # li$PUP6 <- s6@PUP, li$LOD4 <- s4@LOD, etc

        s@LOD <- s@PUP <- NULL
        s <- s %without.name% things
        assign( 's' %&% usy, s)
      }

      # s6@LOD <- s6@PUP <- s4@LOD <- s4@PUP <- s3@LOD <- s3@PUP <- NULL
      # s6$e0 <- s6$v0 <- s6$e1 <- s6$v1 <- NULL
    }

    # Now make "master" variables LOD, PUP, ev01 that correspond to useN
    li[ names( s6)] <- s6 # instead of cbind--- this will overwrite not add
    # Replace the ones that shouldn't be 6way:
    li[ li$useN == 4, names( s4)] <- s4[ li$useN == 4,] ## subs in 4-ways where useN == 4
    li[ li$useN == 3, names( s3)] <- s3[ li$useN == 3,] ## subs in 3-ways where useN == 3
  } else { # ... sneaky override, for non-ABCO systems
    # shouldn't really be called "...6" obvs
    s6$useN <- 6
    li[ names( s6)] <- s6 # instead of cbind--- this overwrites
  }

  li <- make_dull( li, names( li) %that.match% '^ev01')
  li <- make_dull( li, names( li) %that.match% '^(LOD|PUP)[0-9]') # you'll thank make for this :)

  lociar@locinfo <- li
  lociar@hspPower_checksum <- calc_hspPower_checksum( li)
return( lociar)
}




#' Kin-finding power for microhaplotyped loci
#' 
#' This is a short-term fudge for checking HSP-finding power of a bunch of loci
#' that (i) can have as many haplotypes as you like, but (ii) have no errors or
#' nulls. See \bold{Examples} for how you might use it.
#' 
#' If you want to explore the impact of missing genotypes (so that e.g. only
#' 90\% of loci are co-scored in some pairs) you'll have to do so manually. A
#' reasonable and very easy option is to multiply \code{Ediff} and \code{V.UP}
#' both by 0.9, then go thru the steps. If you choose the 0.9 conservatively-
#' ie it's highly likely that >0.9 of loci get co-scored- then the above calc
#' avoids any need to do much more complicated stuff (which I leave to you...).
#' 
#' At some point in future, \code{kinference} might be changed so that it can
#' handle >2 non-null alleles gracefully (ie microhaplotypes). But not yet. So
#' for now this version does some ghastly "live-hacking" of existing code for
#' \code{\link{hsp_power}} to implement no-errors no-nulls multi-allelic case.
#' It will be hard to follow, so use \code{mtrace} if you really want to see
#' what's going on. The guts of the code is in \code{\link{hsp_power}} and
#' \code{predict_hsp_util}.
#' 
#' @param lociar Usually, a matrix of allele frequencies (Locus * Alleles).
#' Locus names are set from the rownames, or "L1", "L2" etc if there are no
#' rownames. Allele names will be set to "A", "B", "C", etc, regardless of
#' colnames; you do not have a choice there. Will be renormalized so rows sum
#' to unity. NB \code{lociar} can also be a \code{snpgeno} object, as expected
#' for \code{hsp_power}. If so, then the allele freqs are assumed to live in
#' \code{lociar$locinfo$pbonzer}, and \bold{no} nulls or genotyping errors are
#' allowed for; hence, for a DartCap "ABCO"-style dataset, \code{hsp_power} and
#' \code{hsp_power2} will give \emph{different} answers.
#' @return If \code{lociar} is an allele-frequency matrix, then you get a
#' dataframe with one row per locus and columns "Ediff", "V.UP", and "sdiff".
#' "Ediff" is "E[LOD|HSP] - E[LOD|UP]"; "V.UP" is "V[LOD|UP]"; "sdiff" is
#' \code{sqrt(V.UP)/Ediff}, useful for ranking locus power. See \bold{Examples}
#' for use.
#' @seealso \code{\link{hsp_power}}
#' @keywords misc
#' @examples
#' 
#' ALF <- matrix( runif( 15), 3, 5) # 3 loci; 5 alleles
#' POW <- hsp_power2( ALF)
#' # look at the contents of each...
#' # Now do it for lots of loci. NB the allele freqs above are *insanely* good; you won't
#' # find anything like that in practice for lots'n'lots of loci
#' lots <- 500
#' ALF <- matrix( runif( lots*5), lots, 5)
#' POW <- hsp_power2( ALF)
#' # Now say we plan 1e6 pairwise comps, and might expect 100 HSPs
#' # Work relative to E[LOD|UP] which is not returned explicitly; treat that as "origin" ie 0
#' V <- sum( POW$V.UP)  # V[PLOD|UP]
#' E <- sum( POW$Ediff) # E[PLOD|HSP] - E[PLOD|UP]
#' E / sqrt( V) # 10.5 SDs--- pretty good.  Mean of HSPs is 10.5 UP-SDs above mean of UPs,
#' # ... so v. unlikely an UP will get as far as _typical_ HSP. But we need to be a bit
#' # ... more stringent than "typical"--- and, NB weaker kin
#' bigUP <- qnorm( 1e-6, mean=0, sd=sqrt( V), lower=FALSE) # most-kinlike UP
#' smallHSP <- qnorm( 1e-2, mean=E, sd=sqrt( 4*V))      # least-kinlike HSP
#' # ... so that's probably OK...
#' 
#' @export hsp_power2
"hsp_power2" <-
function( lociar,
    want_LOD_table=TRUE, # T/F
    k # 0.5 for HSPs
){
############
  # Uses the code of hsp_power, but sneakily redefining a few core functions...
  # ... hsp_power originally intended for 6way null-heavy ABCO-format Dartcap

  just_matrix <- lociar %is.not.a% 'snpgeno'
  if( just_matrix) {
    if( lociar %is.not.a% 'matrix') {
stop( "'lociar' must be 'snpgeno' or matrix of allele freqs (Locus X Alleles)")
    }
    # Turn lociar into a snpgeno with .$locinfo$pbonzer==lociar

    locnames <- rownames( lociar)
    if( is.null( locnames)) {
      locnames <- 'L' %&% seq.int( nrow( lociar))
    }

    lociar[ l, a] := lociar[ l, a] / (SUM_ %[a]% lociar[ l, a]) # normalize
    lociar <- as.array( lociar) # vecless xtensor legacy removal
    dimnames( lociar) <- NULL
    colnames( lociar) <- LETTERS[ seq.int( ncol( lociar))]

    data.frame <- data.frame
    formals( data.frame)$stringsAsFactors <- FALSE
    li <- data.frame( Locus=locnames, pbonzer=I( lociar)) # needs 2-step otherwise R messes up matrix

    genos <- expand.grid( colnames( lociar), colnames( lociar),
        stringsAsFactors=FALSE) %where% (Var1 <= Var2)
    diplos <- genos$Var1 %&% genos$Var2

    lociar <- snpgeno( 1, nrow( lociar), diplos, info=data.frame( Our_sample='THE_THING'),
        locinfo=li)
  }

  e <- new.env( parent=asNamespace( 'kinference'))
  e$lociar <- lociar

  result <- evalq( envir=e, {
    lociar$locinfo$snerr <- 0 # various things expect 'snerr' to exist

    define_genotypes <- eval( substitute( function( nlocal=sys.parent(), eg) mlocal({
        # Fakes various "genotypes_blah" categories that are expected by subsequent routines

        ALLeles <- THE_ALLELES # colnames( pbonzer) # eg A,B,C,D,E; use named() if we need to write A for "A"

        # need stringsAsFactors... FFS WTF were they (R) thinking ????
        eg <- expand.grid( ALLeles, ALLeles, stringsAsFactors=FALSE) %where% (Var1 <= Var2)
        genotypes_C <- eg$Var1 %&% eg$Var2
        genotypes6 <- genotypes_C
        ABCO <- named( ALLeles)
      }), list( THE_ALLELES=colnames( lociar$locinfo$pbonzer))))

    add_pairprob_error <- function( nlocal=sys.parent()) mlocal({
      pp_err <- pp_true                ? 0
      pp6_err <- pp_true               ? 0
    })

    # Make sure subsidiary functions pick up these versions
    environment( hsp_power) <-
        environment( calc_g6probs_IBD0_scalar) <-
        environment( calc_g6probs_IBD1_scalar) <-
        environment( calc_g6probs_IBD2_scalar) <-
        environment( add_pairprob_error) <-
        environment( predict_hsp_util) <-
        environment()
    # mtrace( hsp_power) # automate the debugging a bit
  return( hsp_power( lociar, want_LOD_table=TRUE, k=0.5))
  })

  result$locinfo <- result$locinfo %without.name% cq( snerr, use6)
  if( just_matrix){
    result <- result$locinfo[ cq( Ediff, V.UP, sdiff)]
  }

return( result)



  lociar@locinfo <- li
return( lociar)
}




#' Check individual multilocus genotypes for typicality
#' 
#' \code{ilglk_geno} computes log-likelihood of entire 4-way (not 6-way)
#' genotype of each individual, i.e., sum log Pr[ g(i,l)]; and compares the
#' distribution across individuals to theoretical distro given allele
#' frequencies. Significant mismatch is bad. Can also detect outlier
#' individuals, usually with lglks that are much too low (ie rather than too
#' high- I'm not sure what could generate "too typical a genome" at the
#' individual level).
#' 
#' You can use \code{locator(1)} to click the histogram to figure out where to
#' adjust the \code{xlim}/'ylim' values to change the range of the data to
#' inspect more closely- ie you then re-run the function with its
#' \code{...hist_par} argument set accordingly.
#' 
#' Currently, the SPA calcs are a wee bit slow because of heavy use of
#' \code{vecless} which in version 1.0 is sluggish. The lglks themselves are
#' computed in C and are blisteringly fast. If the SPA line (expected distro)
#' doesn't appear, let us know- needs fixing! There might e.g. be too many
#' loci, so that the calculation is falling over.
#' 
#' Haven't added any formal uh-oh criteria yet; that could be done via the SPA,
#' as in \code{dump_badhetz_fish}. But, reading off from the graph is probably
#' fine...
#' 
#' @param snpg a \code{snpgeno} (6-way genotype)
#' @param hist_pars \code{list()} passed to \code{hist} for controlling
#' histogram, eg \code{hist_pars=list(xlim=c(-12000, -6000))}, or use
#' \code{FALSE} to not plot.
#' @param showPlot show the histogram? Defaults to TRUE, but overrideen by
#' \code{hist_pars=FALSE}.
#' @return Vector of log-likelihood for each individual; also usually (but
#' optionally), a histogram of log-likelihood values across individuals.
#' @keywords misc
#' @examples
#' 
#' ## Need an example!
#' 
#' @export ilglk_geno
"ilglk_geno" <-
function(snpg, hist_pars=list(), showPlot = TRUE) {
  define_genotypes()
  extract.named( snpg@locinfo[ cq( pbonzer)])

  p0 <- pbonzer[,'O'] + pbonzer[,'C']
  pA <- pbonzer[,'A']
  pB <- pbonzer[,'B']

  n_samps <- nrow( snpg)
  n_loci <- ncol( snpg)
  minfo_fields <- cq( Our_plate, Our_sample) %that.are.in% names( snpg$info)
  snpg4 <- snpgeno(
      NULL,
      diplos=genotypes4_ambig,
      info= snpg@info[, minfo_fields],
      locinfo= snpg@locinfo[, 'Locus', drop=FALSE]
    )
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

  quick <- TRUE # used to be param
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

  dens_SPA <- renorm_SPA( K, dK, ddK, return='func', already_vectorized=TRUE)

  if( isFALSE( hist_pars)) {
    showPlot <- FALSE
  } else {
    hist_pars <- add_list_defaults( hist_pars,
        main   = 'Geno lglk by specimen',
        xlim   = range( ilglk),
        col    = "grey",
        border = NA,
        xlab   = '',
        nclass = 50)
  }

  if( showPlot) {
      lv <- do.call( 'hist', c( list( x=ilglk), hist_pars))

      # some dens_SPA can fail -- do lapply and then weed out baddies
      mids_SPA <- lapply(lv$mids, function(x) try(dens_SPA( x), silent=TRUE))
      good_ind <- unlist(lapply(mids_SPA, class) != "try-error")
      mids_SPA <- unlist(mids_SPA[good_ind])
      # plot predicted density. Slowish with vecless 1.0
      lines( lv$mids[good_ind], diff( lv$breaks)[good_ind] * mids_SPA * n_samps,
            col='green')
  }

return( ilglk)
}


"indiv_lglk_geno" <-
function(lpgeno, geno) {
    .Call(`_kinference_indiv_lglk_geno`, lpgeno, geno)
}


"K_indiv" <-
function(tt, geno, vec_LOD, Pg) {
    .Call(`_kinference_K_indiv`, tt, geno, vec_LOD, Pg)
}




#' add a kin-type legend with the default colour scheme
#' 
#' Package kinference uses a constant colour scheme for kin types, designed to
#' be colourblind-friendly and to allow clear visual distinction between kin
#' classes that share similar relatedness levels. We recommend that users
#' maintain this colour scheme in any custom plots they create. This utility
#' function adds a legend with kinship acronyms and their colours. By default,
#' it will display only the "top four" kin classes (UP, FSP, POP, and HSP)
#' 
#' 
#' @param position passed to \code{legend}. Must be one of "topleft", "top",
#' "topright", etc.
#' @param include a character vector of kin-classes to include. Limited to
#' "POP", "GGP", "HCP", "FCP", "UP", "HTP", "FTP", "HSP", and "FSP". Will
#' automatically include "UP", "POP", "FSP", and "HSP" unless these are
#' specifically \emph{excluded}. Will trigger a warning if any values are also
#' in \code{exclude}, as those values will be excluded.
#' @param exclude a character vector of kin-classes to exclude. Limited to
#' "POP", "GGP", "HCP", "FCP", "UP", "HTP", "FTP", "HSP", and "FSP". Will
#' trigger a warning if any values are also in \code{include}.
#' @param ... additional args, passed to legend()
#' @return adds a legend to the current plot
#' @keywords misc
#' @export kinlegend
"kinlegend" <-
function(position = "topright", include = character(), exclude = character(), ...) {

    palette(c("#0D0887FF", "#48039FFF", "#7401A8FF", "#9D189DFF", "#BF3984FF",
              "#DA596AFF", "#EE7B51FF", "#FBA238FF", "#FCCE25FF"))

    if(sum(include %in% exclude) > 0) {
        stop(paste("kin class(es) '",
                   paste(include[include %in% exclude], collapse = ", "),
                   "' is/are in both 'include' and 'exclude'", sep = "")
             )
    }

    if(!all(include %in% c("POP","GGP","HCP","FCP","UP","HTP","FTP","HSP","FSP"))) {
        warning("some unrecognised kin classes in 'include'")
    }
    if(!all(exclude %in% c("POP","GGP","HCP","FCP","UP","HTP","FTP","HSP","FSP"))) {
        warning("unrecognised kin classes in 'exclude'")
    }

    ## the default set
    allNames <- c("UP","POP","HSP","FSP")

    ## add any extras from 'include'
    allNames <- unique(c(include, allNames))

    ## remove those in 'exclude'
    allNames <- allNames[!allNames %in% exclude]

    ## find their palette numbers
    allNumbers <- match(allNames, c("POP","GGP","HCP","FCP","UP","HTP","FTP","HSP","FSP"))

    legend(position, legend = allNames,
               lwd = 2, lty = 1, col = allNumbers, bg = "white", ...)
}


"lglk_loci" <-
function( snpg) {
# Obsolete / not useful ?
# When loci are presumably misbehaving, in that ilglk_geno looks
# ... OK on a "trusted subset" but not on all--- this may help

#% 'lglk_loci' compares, for each locus, the average (across individuals) observed lglk with the theoretical mean and variance. The idea is to help figure out when some loci
#are going wrongish (e.g. you can get decent fits from a subset of loci). Of course, 'check6and4' pvals should be the main guide here; 'lglk_loci' can show an overall deviation
#, as well as any remaining locus-specific misbehaviour (but shouldn't be much locus-specific stuff thx2 'check6and4'). Overall too-good-to-be-true-ism (as seen for 'Glyphis
#garricki') _might_ come when ALF is estimated from very small datasets.

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
}


"make_pgeno" <-
function( pA, pB, pC, which_genotypes) {
###
# which_genotypes eg genotypes_ambig. If no C-genos requested then pC is forced to 0 so new O includes C
  if( !any( grepl( 'C', which_genotypes))) {
    pC <- 0
  }

  pO <- pmax( 0, 1 - pA - pB - pC) # rounding error guard
  pAO <- 2*pA*pO
  pAA <- sqr( pA)
  pAAO  <- pAA + pAO
  pAB  <- 2*pA*pB
  pBO <- 2*pB*pO
  pBB <- sqr( pB)
  pBBO  <- pBB + pBO
  pAC <- 2*pA*pC
  pCO <- 2*pC*pO
  pCC <- sqr( pC)
  pCCO <- pCC + pCO
  pBC <- 2*pB*pC
  pOO  <- pO^2
  pBBOO <- pBBO + pOO

  # Could be scalar or vector: c or cbind
  funco <- if( length( pO) > 1) cbind else c
  phat <- do.call( funco, FOR( which_genotypes, get( 'p' %&% .)))
  names( phat) <- sub( '[.].*', '', names( phat)) # loco R name-extrusion habit FFS
return( phat)
}


"map6to3" <-
function(g6p0, g6p1, g6p2){
  define_genotypes()

  map6to3 <- matrix( 0, 6, 3, dimnames=list( genotypes6, genotypes3))
  # AB & OO are OK; AAO should receive both AA and AO; etc
  mm <- match( genotypes6, substring( genotypes3, 1, 2), 0) # the "AA" bit of "AAO"...
  yup <- cbind( which(mm>0), mm[ mm>0])
  map6to3[ yup] <- 1
  mm <- match( genotypes6, substring( genotypes3, 2, 3), 0) # ... and the "AO" bit
  yup <- cbind( which(mm>0), mm[ mm>0])
  map6to3[ yup] <- 1
  mm <- match( genotypes6, substring( genotypes3, 3, 4), 0) # ... and the "OO" bit
  yup <- cbind( which(mm>0), mm[ mm>0])
  map6to3[ yup] <- 1

  # Really want g4p0[l,i,j] := map6to4[i,k6] %[k6]% g6p0[l,k6,m6] %[m6]% map6to4[m6,j]
  # ... but vecless can't presently handle multi-stages

  A[l,i,k] := g6p0[l,i,j] %[j]% map6to3[j,k]
  g3p0[l,i,j] := map6to3[ k,i] %[k]% A[l,k,j]

  A[l,i,k] := g6p1[l,i,j] %[j]% map6to3[j,k]
  g3p1[l,i,j] := map6to3[ k,i] %[k]% A[l,k,j]

  A[l,i,k] := g6p2[l,i,j] %[j]% map6to3[j,k]
  g3p2[l,i,j] := map6to3[ k,i] %[k]% A[l,k,j]


returnList( g3p0, g3p1, g3p2)
}


"map6to4" <-
function(g6p0, g6p1, g6p2){
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

  A[l,i,k] := g6p2[l,i,j] %[j]% map6to4[j,k]
  g4p2[l,i,j] := map6to4[ k,i] %[k]% A[l,k,j]

returnList( g4p0, g4p1, g4p2)
}


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
}


"OLD_find_POPs" <-
function( snpg, subset1=1 %upto% nrow( snpg), subset2=subset1,
    limit_pairs=0.5*nrow(snpg),
    keep_thresh,
    eta= NULL,
    nbins= 50,
    maxbin = NULL,
    WPSEX_UP_POP_balance=0.99) {
###################
## This is the non-lglk WPSEX/NABOO version

  # Sanity...
define_genotypes()
stopifnot( is.numeric( subset1) && is.numeric( subset2))
stopifnot( all( !duplicated( subset1)) && all( !duplicated( subset2)))
stopifnot( my.all.equal( subset1, subset2) || !length( intersect( subset1, subset2)))
stopifnot( (snpg %is.a% 'snpgeno') && (my.all.equal( snpg@diplos, genotypes6) || my.all.equal( snpg@diplos, genotypes4_ambig)))

  ## # SB checks--- but NB rowSums( pbonzer==1) !
  ## hspPower_change <- snpg@hspPower_checksum != with(snpg@locinfo, sum(pbonzer, snerr, useN))
  ## PLODSPA_change <- snpg@PLODSPA_checksum != with(snpg@locinfo, sum(useN, LOD6, LOD4, LOD3,
  ##                                                 PUP6, PUP4, PUP3))
  ## if(hspPower_change | PLODSPA_change) {
  ##   warning("snpg$locinfo appears to have been modified after hsp_power and/or prepare_PLOD_SPA were last called. I sure hope you know what you're doing...")
  ## }

  # Decide based #apparent exclusions of AA/BB form, using 4way genos, though
  # ... it's really AAO/BBO so not a true exclu but
  # ... close among pop-loci
  # Sticking with 4way genos so that genotyping errors are low

  # Instead of Nexclu, uses a wted sum of "exclus" in 4way genos to max expected diff between POP and UP
  # This version doesn't allow for geno errors, but does realize that AO/BO could happen
  # Doesn't bother with OO-AB

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

  # Algo should work either with 4way or 6way genotypes (no C allele though)
  # C code gets told, via AAish and BBish, which codes are really AAO and BBO
  # For 6way, remap AO->AA, BO->BB, and tell C just to check for AA and BB
  temp_snpg <- snpg[ , pop_loci]
  if( my.all.equal( snpg@diplos, genotypes6)) {
    recode_6_as_pseudo4 <- function( x) { x[ x=='AO'] <- AA; x[ x=='BO'] <- BB; x}
    temp_snpg <- recode_6_as_pseudo4( temp_snpg)
    AAish <- match( 'AA', snpg@diplos)
    BBish <- match( 'BB', snpg@diplos)
  } else { # 4way
    AAish <- match( 'AAO', snpg@diplos)
    BBish <- match( 'BBO', snpg@diplos)
  } # else NYI, which would have been picked up in sanity checks

  # Remove extranea
  attributes( temp_snpg) <- attributes( temp_snpg)[ 'dim']
##  attributes( temp_snpg) <- temp_snpg@dim # that's yer lot; just the scores
  temp_snpg <- t( temp_snpg)

  n_loci <- ncol( snpg)

  # Prepare for diagnostics of #excl
    ## Defaults to two-sigma above *UP* mean... which will cover almost all

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

    K <- compile_vecless( K(0))
    dK <- compile_vecless( dK(0)) ## dK(0) should be the mean of wpsex for true UPs
    ddK <- compile_vecless( ddK(0)) ## should be the var of ditto

    #  n_sim_check <- 1000
    #  Ktest <- function( tt) {
    #    x <- matrix( runif( n_sim_check * n_loci) < pex_up, n_loci, n_sim_check)
    #    ewx <- exp( x*ww*tt)
    #    colSums( log( ewx))
    #  }
    #

    # Could do now predict binprobs *for UPs) via renorm_SPA_cumul
    # Let's not
    # What we need, is...
  if( is.null( maxbin)) {
      maxbin  <-  dK(0)+2*sqrt(ddK(0))
  }

  if(is.null(eta) ) {
      eta <- dK(0)-3*sqrt(ddK(0))
  }

  if( nbins<2) {
    warning( 'nbins<2 is senseless; gonna use 2')
    nbins <- 2
  }
  bins <- seq( from= 0, to= maxbin, length=nbins)
  ## binterval <- maxbin / (nbins-1) # bin *starting* at maxbin counts all bigger
                                        # binprobs not set
  binterval <- bins[2] - bins[1]

  # Trying special-cases here to minimize copying
  symmo <- my.all.equal( subset1, subset2)
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
        keep_n = limit_pairs,
        AAO= AAish,
        BBO= BBish,
        nbins = nbins,
        binterval = binterval
      )
  } else { # different subsets
    result <- POP_wt_paircomps_lots(
        geno1= temp_snpg[ ,subset1],
        geno2= temp_snpg[ ,subset2],
        w= ww,
        symmo= FALSE,
        eta= eta,
        max_keep_wpsex= keep_thresh,
        keep_n = limit_pairs,
        AAO= AAish,
        BBO= BBish,
        nbins = nbins,
        binterval = binterval
      )
  }

  # warning if we're running up against storage constraints
  if(length(result$big_wpsex) == limit_pairs){
    message("Returning just the ", limit_pairs, " pairs with the most POP-like wpsex scores, increase 'limit_pairs' if more are required")
  }

  n_wpsex_in_bin <- result$n_wpsex_in_bin  ## SB tweak; the 'below zero' bin must be empty
  # bins <- seq(0+binterval, 0+(nbins*binterval), binterval)
  ## SB addition. Bloody zero-base.

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
  attributes( result) <- c( attributes( result), returnList(
      bins, n_wpsex_in_bin, eta, keep_thresh))

  result@n_loci <- length( pop_loci)
  result@mean_UP <- dK( 0)
  result@var_UP <- ddK( 0)
  result@call <- sys.call()

return( result)
}


"OLD_split_FSPs_from_POPs" <-
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


"paircomps" <-
function(pair_geno, LOD, geno1, geno2, symmo, granulum, granulum_loci) {
    .Call(`_kinference_paircomps`, pair_geno, LOD, geno1, geno2, symmo, granulum, granulum_loci)
}




#' PLOD histogram on log-scale
#' 
#' Plots a log-frequency histogram for the output of \code{find_HSPs()}, with
#' the expected mean PLOD for unrelated pairs, the expected distribution of
#' unrelated pairs, and the expected mean PLOD for HSPs. Expectations are
#' coloured according to the table below
#' 
#' Colour scheme for all kin-finding markers: Shamelessly cropped from package
#' \pkg{viridis}; defined as a one-off palette to avoid adding dependencies.
#' Apparently quite colourblind-friendly.
#' 
#' Kin class Hex Colour Number Colour UP #BF3984FF 5 Magenta (light) POP
#' #0D0887FF 1 Navy blue GGP #48039FFF 2 Violet HSP #FBA238FF 8 Orange FSP
#' #FCCE25FF 9 Yellow HCP #7401A8FF 3 Purple FCP #9D189DFF 4 Magenta (dark) HTP
#' #DA596AFF 6 Rose FTP #EE7B51FF 7 Coral
#' 
#' @param hsps the output of a call to \code{find_HSPs()}
#' @param UP,HSP,POP,FSP whether plot the expected (mean) PLOD for pairs of
#' that type? Defaults TRUE
#' @param showUP plot the expected density curve for unrelated pairs using the
#' SPA approximation (default TRUE), Normal approximation (default FALSE),
#' both, or neither. Either approximation will plot in colour 5, a light
#' magenta.
#' @param ... additional pars, passed to \code{plot}
#' @keywords misc
#' @export PLOD_loghisto
"PLOD_loghisto" <-
function(
  hsps,
  UP= TRUE,
  HSP= TRUE,
  POP= TRUE,
  FSP= TRUE,
  showUP= c(SPA= TRUE, Normal= FALSE),
  ...
){
  palette(c("#0D0887FF", "#48039FFF", "#7401A8FF", "#9D189DFF", "#BF3984FF",
            "#DA596AFF", "#EE7B51FF", "#FBA238FF", "#FCCE25FF"))

  binmids <- hsps@bins + (hsps@bins[2] - hsps@bins[1])/2

  ## c++ gives bins that are out in a different direction for bins and n_in_bin
  x <- hsps@bins
  y <- hsps@n_PLODs_in_bin

  l <- list( ...)
  ylim <- l$ylim
  ylim <- if( is.null( ylim)) c( 0, NA) else c( 0, ylim[ 2])

  plot( x[-1], log10( head( pmax( y, 0.1), -1)),
     ...,
     ylim= ylim,
     type= "S", xlab= "PLOD", ylab= "log10(Frequency)")
  # was: hsps@bins[-1], log10( hsps@n_PLODs_in_bin[1:(length(hsps@n_PLODs_in_bin)-1)]),

  # Pass thru _some_ graphical pars
  # linargs <- list( ...)[ cq( lty, lwd, col)]
  # linargs <- linargs %SUCH.THAT% length(.) # eliminate pure NULLs
  # ... which makes the lines() call look confusing...
  # do.call( 'lines', c( list(
  #     tail( x, -1), log10( y[-1]), type='s'),
  #    linargs))


  if( UP) { abline(v= hsps@mean_UP, col= 5, lwd= 2) }
  if( HSP) { abline(v= hsps@mean_HSP, col= 8, lwd= 2) }
  if( POP) { abline(v= hsps@mean_POP, col= 1, lwd= 2) }
  if( FSP) { abline(v= hsps@mean_FSP, col= 9, lwd= 2) }
  if( showUP["SPA"]) {
      lines( binmids,log10( diff( hsps@binprobs)*sum( y[ binmids<0])),
            lwd=2,col=5)
  }
  if( showUP["Normal"]) {
    lines( x, log( diff( c( 0,
        pnorm( binmids, mean= hsps@mean_UP, sd= sqrt( hsps@var_UP)) *
        sum( y[ binmids<0])))),
        lwd= 2, col= 5) ## Normal approx
  }

  # MVB pu-R-ist recode here ;) to avoid heavvvy data.frame
  legendBits <- c( UP=5, POP=1, HSP=8, FSP=9)
  legendBits <- legendBits[ unlist( mget( names( legendBits)))] # T or F for each
  legend("topright",
      legend= names( legendBits),
      col= legendBits,
      lwd= 2, lty= 1, bg= "white")

  if( FALSE){ # old code
    legendBits <- data.frame(allNames= c("UP","POP","HSP","FSP"), allNumbers= c(5,1,8,9))
    if(!UP) {
        legendBits <- legendBits[legendBits$allNames != "UP",]
    }
    if(!POP) {
        legendBits <- legendBits[legendBits$allNames != "POP",]
    }
    if(!HSP) {
        legendBits <- legendBits[legendBits$allNames != "HSP",]
    }
    if(!FSP) {
        legendBits <- legendBits[legendBits$allNames != "FSP",]
    }
    legend("topright", legend= legendBits$allNames,
           lwd= 2, lty= 1, col= legendBits$allNumbers, bg= "white")
  }
}




#' Oddness metrics
#' 
#' Plots the percentage of all pairs in each bin with an \code{unusually} low
#' CLOD score, ilglk stat, or hetz stat, across the range of PLOD.
#' \code{PLOD_oddness_oneway()} shows the percentage of cases where one member
#' has a low score, and \code{PLOD_oddness_twoway()} shows the percentage of
#' cases where both members of the pair have a low score.
#' 
#' 
#' @aliases PLOD_oddness_oneway PLOD_oddness_twoway
#' @param hsps the output of a call to \code{find_HSPs()}
#' @param snpg the \code{snpgeno} object from which \code{hsps} was built
#' @param lb PLOD lower bound for plot extent. Should exclude the UP bump
#' @param ub PLOD upper bound for plot extent. Defaults to the maximum PLOD
#' score plus a little padding
#' @param bin hist bin width. Default 5
#' @param CLOD_prop the quantile of CLOD below which animals are highlighted.
#' Default 0.001
#' @param ilglk_prop the quantile of ilglk stat below which animals are
#' highlighted. Default 0.001
#' @param hetz_prop the quantile of hetz stat below which animals are
#' highlighted. Default 0.001
#' @param ... additional pars, passed to \code{plot()}. \code{ylim} and
#' \code{breaks} are set internally, so you cannot pass them via \code{...}.
#' @seealso HSP_oddness_oneway
#' @keywords misc
#' @export PLOD_oddness_oneway
"PLOD_oddness_oneway" <-
function(hsps, snpg, lb = min(hsps$PLOD)-10, ub = max(hsps$PLOD)+10, bin = 5,
             CLOD_prop = 0.001, ilglk_prop = 0.001, hetz_prop = 0.001, ...) {

        palette("default")

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
             pch=(15), xlab="PLOD value", ylab="Percent of specimen pairs",
             ylim = c(0, yup), ...)
        points(hist1$breaks[-1], histB$counts/hist1$counts, type='b',col=(3), pch=(16))
        points(hist1$breaks[-1], histC$counts/hist1$counts, type='b',col=(4), pch=(17))
        legend('topright',c('low CLOD stat','low ilglk stat','low hetz stat'),
               title='ONE member of the pair has:', pch=15:17,col=2:4)
    }


"PLOD_oddness_twoway" <-
function(hsps, snpg, lb = min(hsps$PLOD)-10, ub = max(hsps$PLOD)+10, bin = 5,
             CLOD_prop = 0.001, ilglk_prop = 0.001, hetz_prop = 0.001, ...) {

        palette("default")

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
             pch=(15), xlab="PLOD value", ylab="Percent of specimen pairs",
             ylim = c(0, yup), ...)
        points(hist1$breaks[-1], histB$counts/hist1$counts, type='b',col=(3), pch=(16))
        points(hist1$breaks[-1], histC$counts/hist1$counts, type='b',col=(4), pch=(17))
        legend('topright',c('low CLOD stat','low ilglk stat','low hetz stat'),
               title='BOTH members of the pair have:', pch=15:17,col=2:4)
    }


"POP_paircomps_lots" <-
function(geno1, geno2, symmo, eta, max_keep_Nexclu, keep_n, bins, AAO, BBO) {
    .Call(`_kinference_POP_paircomps_lots`, geno1, geno2, symmo, eta, max_keep_Nexclu, keep_n, bins, AAO, BBO)
}


"POP_wt_paircomps_lots" <-
function(geno1, geno2, w, symmo, eta, max_keep_wpsex, keep_n, AAO, BBO, nbins, binterval) {
    .Call(`_kinference_POP_wt_paircomps_lots`, geno1, geno2, w, symmo, eta, max_keep_wpsex, keep_n, AAO, BBO, nbins, binterval)
}




#' Check POP/FSP splitter; private.
#' 
#' This doco shouldn't be visible in the package! (KEYWORDS internal)
#' 
#' \code{simcheck_FSP_POP} simulates POPs and FSPs based on a known set of
#' loci, ready for checking \code{find_FSPs_from_POPs} or variants. Quick.
#' 
#' \code{postprocess_simcheck_FSP_POP} can be used to graphically confirm that
#' the analytical probabilities of samenames and pseudoexclusion actually match
#' the simulated values- it is amazingly difficult to get the formulae right.
#' You wouldn't generally need \code{postprocess...} if you are working with a
#' real set of loci.
#' 
#' @aliases postprocess_simcheck_FSP_POP simcheck_FSP_POP
#' @param forp_sim result of previous call to \code{simcheck_FSP_POP}
#' @param findings result of previous call to \code{find_FSPs_from_POPs_v2(
#' forp_sim, ..., keep_indiv=TRUE)}- which generates the \code{is_same} and
#' \code{is_psex} matrices for (pairs * loci), and the by-locus vectors
#' \code{Pr_same_FSP}, \code{Pr_psex_FSP}, \code{Pr_same_POP},
#' \code{Pr_psex_POP}.
#' @param plot. whether to plot. Can't see why this would ever be false...
#' @param snpg a \code{snpgeno} object with known allele freqs and error rates
#' (i.e. \code{pbonzer} and \code{snerr} must exist in \code{snpg$locinfo})
#' @param chromos Controls the extent of linkage:how many "chromosomes" should
#' the loci be split between? There's no crossover within these "chromosomes",
#' but they inherit independently of each other. Loci are allocated to
#' chromosomes in sequence, so that if there are 3 chromosomes, the first locus
#' goes to "C1", the second to "C2", third to "C3", fourth to "C1", etc.
#' Returned object will have an extra field \code{$locinfo$chromosim}.
#' @param N number of pairs of each type (ie this many FSPs, and this many
#' POPs).
#' @return \code{simcheck_FSP_POP} returns a \code{snpgeno} object with
#' \code{4*N} rows, and an attribute "callo" recording the call. The
#' \code{$info} contains just one column, "Our_sample", with names like "F17_A"
#' and "P6_B"; those would respectively indicate the first (A) member of the
#' 17th FSP, and the second (B) member of the 6th POP. The \code{locinfo} field
#' is unchanged, except as noted under "chromos" above.
#' \code{postprocess_simcheck_FSP_POP} just plots.
#' @seealso \code{find_FSPs_from_POPs}
#' @keywords misc
#' @examples
#' 
#' ForP <- simcheck_FSP_POP( my_snpg, 10, 1000)
#' ForP$locinfo$useN[] <- 3 # minimal but most robust genotype
#' test3 <- find_FSPs_from_POPs_v2( ForP,
#'     cbind( seq( 1, by=2, length=1000),
#'       seq( 1, by=2, length=1000)),
#'     keep_indiv=TRUE) # need keep_indiv=T for postprocess...() to work
#' postprocess_simcheck_FSP_POP( ForP, test3)
#' hist( test3, nc=50)
#' abline( v=test$E_FPstat, col='blue') # POP and FSP means
#' \dontrun{
#' # Sim from real loci--- in this case, with 'useN==4' for all loci since 6 looked a bit iffy
#' library( atease)
#' s11nodup_all4 <- s11nodup
#' s11nodup_all4$locinfo$use6 <- NULL
#' s11nodup_all4$locinfo$useN <- 4L
#' simbo4 <- simcheck_FSP_POP( s11nodup_all4, N=1000, chromo=20)
#' simbo4_next <- find_FSPs_from_POPs_v2( simbo4, cbind( seq( 1, 3999, by=2), seq( 2, 4000, by=2)), keep=T)
#' kinference:::postprocess_simcheck_FSP_POP( simbo4, simbo4_next) # looks OK; not needed by "user"
#' hist( simbo4_next$FPstat, nc=50)
#' abline( v=simbo4_next@E_FPstat, col='green') # theory means
#' # 95% of all POPs should be Left of the line drawn next:
#' abline( v=simbo4_next@E_FPstat['POP']+2*sqrt( simbo4_next@V_FPstat), col='blue', lty=1)
#' # Can't say for FSPs, becoz linkage
#' # Real data:
#' testo4 <- find_FSPs_from_POPs( s11nodup_all4, pops_005)
#' abline( v=testo4$FPstat, col='red')
#' }
#' 
#' @export postprocess_simcheck_FSP_POP
"postprocess_simcheck_FSP_POP" <-
function( forp_sim, findings, plot.=TRUE) {
  # Are the names what we expect from simcheck_FSP_POP ..?
  stopifnot( all( grepl( '^[FP][0-9]+_[AB]$', forp_sim$info$Our_sample)))

  # ... doesn't proof against someone scrambling the order, of course!
  extract.named( attributes( findings)[ grep( '_same|_psex', atts( findings), value=TRUE)] )

  N <- nrow( is_same) / 2
  FSPs <- 1:N
  POPs <- N + FSPs

  if( plot.) {
    opar <- par( no.readonly=TRUE)
    on.exit( par( opar))

    # Code from 'handy::clear.graphs'
    par(mfrow = c(1, 1))
    plot(0, 0, type = "n", axes = F, xlab = "", ylab = "")
    par(mfrow = c(2,2), pty='s')
  }

  useN <- forp_sim$locinfo$useN

  ucols <- c( useN3='blue', useN4='orange', useN6='lightgreen')

  # Generate emp_Pr_psex_POP etc
  for( kinno in cq( FSP, POP)) {
    for( compo in cq( same, psex)) {
      emp <- colMeans( get( 'is_' %&% compo)[ get( kinno %&% 's'),])
      assign( sprintf( 'emp_Pr_%s_%s', compo, kinno), emp)
      if( plot.) {
        ana <- get( sprintf( 'Pr_%s_%s', compo, kinno))
        plot( ana, emp, xlim=range( c( ana, emp)), ylim=range( c( ana, emp)),
          xlab='', ylab='', main='', type='n', col='grey', asp=1)
        points( ana, emp, col=ucols[ 'useN' %&% useN], pch='.', cex=2)
        legend( 'topleft', legend=as.character( c( 3,4,6)),
            col=c( 'blue', 'orange', 'lightgreen'), pch=rep( '.', 3), pt.cex=10)
        abline( 0, 1)
        title( sprintf( '%s %s X: theory Y: data', compo, kinno), sub=1)
      }
    }
  }

return( NULL)
}


"predict_hsp_util" <-
function( pIBD0, pIBD1, pIBD2, want_LOD_table=FALSE, k=0.5) {
  # This version ignores the possibility of errors involving AB or OO...
  # ... which should be pretty rare

  define_genotypes()
  nl <- nrow( pIBD1)
  Phsp <- pIBD1 * k + pIBD0 * (1-k)
  Pup <- pIBD0
  Ppop <- pIBD1
  Pfsp <- pIBD0*0.25 + pIBD1*0.5 + pIBD2*0.25

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

  E.POP[l] := LOD[l,i,j] %[i,j]% Ppop[l,i,j]
  E.FSP[l] := LOD[l,i,j] %[i,j]% Pfsp[l,i,j]

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

    retval <- data.frame( Ediff, V.UP, sdiff, matto,
                         E.POP, E.FSP, E.UP, E.HSP)
    ## matto comes in as 4 named columns, not as matto.
    ## E.UP is not 100% efficient, but bonus points for readability
  if( want_LOD_table) {
    retval@LOD <- gpLOD
    retval@PUP <- gpPUP
  }
return( retval)
}




#' Prepare for kin-finding
#' 
#' \code{prepare_PLOD_SPA} is something you have to run before using some
#' kin-finding/QC tools, to set up your \code{snpgeno} object for fancy maths
#' woooo (saddlepoint approximations). There are no meaningful options, you
#' just have to run this. It can be \emph{slightly} slow which is why it's a
#' separate step.
#' 
#' 
#' @param geno6 a \code{snpgeno} object that has been thru \code{hsp_power}
#' @param n_pts_SPA_renorm how accurate to make the approximation. Default
#' should be fine.
#' @return Another \code{sngeno} object with an environment \code{Kenv}, which
#' contains functions (with their own preloaded data) allowing null
#' distributions (eg PLODs for true UPs) to be calculated. Various sanity
#' checks are incorporated to try to stop you from stuffing up with
#' out-of-synch loci etc later.
#' @keywords misc
#' @export prepare_PLOD_SPA
"prepare_PLOD_SPA" <-
function( geno6, n_pts_SPA_renorm=201) {
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
      e$sqr <- function( x) x*x
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

    oldClass( geno6) <- 'snpgeno' # no special SPA class now
    ## after subset operations.

    geno6@PLODSPA_checksum <- calc_PLODSPA_checksum( geno6$locinfo)

  # PUP and LOD are in Kenv now, so don't duplicate them in locinfo
  # They are really the "workhorse" versions, and are a bit cheaty, so don't want
  # them too public

return( geno6)
}




#' Re-estimate allele frequencies after read-in with load_whopper.
#' 
#' Returns locus-specific allele frequency estimates.
#' 
#' 
#' @param snpg an snpg object
#' @keywords misc
#' @export re_est_ALF
"re_est_ALF" <-
function( snpg) {
## check to be called after load_whopper loads entire dataset

  define_genotypes()
  n_samp <- nrow( snpg)
  n_loci <- ncol( snpg)
    gamb <- snpg # includes C but won't be used
    gamb@diplos  <- genotypes_ambig
  gamb[ snpg==AB] <- AB
  gamb[ snpg==OO] <- OO
  gamb[ snpg==AO] <- AAO
  gamb[ snpg==AA] <- AAO
  gamb[ snpg==BO] <- BBO
  gamb[ snpg==BB] <- BBO

  ## snpg@geno_amb <- gamb # required by...
  new_ALFs <- est_ALF_ABCO( snpg, geno_amb = gamb)
return( new_ALFs)
}


"set_thresholds" <-
function( keeping, nlocal=sys.parent()) mlocal({
## Used to exist to auto-set keep_thresh based on one_in_X_eta
## Obsolete AFAIK Dec 2020
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


"simcheck_FSP_POP" <-
function( snpg, chromos, N, check_genofreq=FALSE) {
## For internal checks only AFAICR
## For now, we should leave the function in the package, but not exported.
## "internal" keyword in doco _should_ enforce that, but...
## Anyway, Rd doc is now included at the end of the code here, for future ref

## Simulate N FSPs and N POPs using loci in snpg
## Check whether find_FSPs_from_POPs works OK

  extract.named( snpg$locinfo[ cq( pbonzer, perr, snerr)])
  n_loci <- nrow( perr)

  # Assign each locus to a chromo
  reppo <- n_loci %/% chromos
  mychro <- rep( 1:chromos, reppo+1)[ 1:n_loci]

  thrub <- snpg[ rep( 1, 4*N),]
  thrub$info <- thrub$info[, 'Our_sample', drop=FALSE]
  thrub$locinfo <- thrub$locinfo[ cq( Locus, pbonzer, snerr, perr, useN, use6, PUP6) %that.are.in%
    names( thrub$locinfo)]
  thrub$locinfo$chromosim <- 'C' %&% mychro

  # Sample names: F1_A, F1_B, F2_A, F2_B, ..., P1_A, P1_B, etc
  ForP <- c( rep( 'F', 2*N), rep( 'P', 2*N))
  pairnum <- rep( rep( 1:N, each=2), 2)
  AB <- rep( LETTERS[1:2], 2*N)
  thrub$info$Our_sample <- sprintf( '%s%i_%s', ForP, pairnum, AB) # ... not to be confused with alleles A&B !

  # Start by giving EVERYTHING a different chromo


#  pABO <- pbonzer[ , 1:3]
# pABO[,3] <- pbonzer[, 'C'] + pbonzer[,'O']
#  g1 <- g2 <- matrix( integer(0), 4*N, n_loci)
#  for( iloc in 1:n_loci) {
#    g1[ , iloc] <- rsample( 4*N, cq( A, B, O), prob=pABO[ iloc,], replace=TRUE)
#    g2[ , iloc] <- rsample( 4*N, cq( A, B, O), prob=pABO[ iloc,], replace=TRUE)
#  }

  get_one_copy <- function() {
      g <- matrix( 'O', 4*N, n_loci)
      g_A <- matrix( runif( 4*N*n_loci), 4*N, n_loci) < rep( pbonzer[,'A'], each=4*N)
      g_B_if_not_A <- matrix( runif( 4*N*n_loci), 4*N, n_loci) < rep( pbonzer[,'B']/(1-pbonzer[,'A']), each=4*N)
      g[ g_A] <- 'A'
      g[ !g_A & g_B_if_not_A] <- 'B'
    return( g)
    }

  g1 <- get_one_copy()
  g2 <- get_one_copy()

  # Merge genos that are IBD--- per chromo
  # FSPs
  FSPstart <- seq( 1, by=2, length=N)

  ibd_chromo1 <- matrix( runif( N*chromos) > 0.5, N, chromos)
  ibd1 <- ibd_chromo1[ , mychro]
  g1[ FSPstart+1, ][ ibd1] <- g1[ FSPstart,][ ibd1]
  ibd_chromo2 <- matrix( runif( N*chromos) > 0.5, N, chromos)
  ibd2 <- ibd_chromo2[ , mychro]
  g2[ FSPstart+1, ][ ibd2] <- g2[ FSPstart,][ ibd2]

  # POPs--- just make the g1's the same
  POPstart <- 2*N + seq( 1, by=2, length=N)
  g1[ POPstart+1,] <- g1[ POPstart,]

  # Merge g1&g2 into one genotype per fish/locus
  swappo <- g1 > g2
  g3 <- g1[ swappo]
  g1[ swappo] <- g2[ swappo]
  g2[ swappo] <- g3

  define_genotypes()
  thrub$locinfo@diplos <- genotypes6
  thrub[,] <- g1 %&% g2
  # thrub$locinfo$useN <- 6 # user can change post hoc

  # Apply XO/XX errors
  # Could use 'perr' directly (symm errors XO/XX), but safer/clearer to use snerr?
  miscall <- function( which_snerr) {
      matrix( runif( 4*N*n_loci), 4*N, n_loci) < rep( snerr[ , which_snerr], each=4*N)
    }

  isAA <- thrub==AA
  isAO <- thrub==AO
  thrub[ isAA & miscall( 'AA2AO')] <- AO
  thrub[ isAO & miscall( 'AO2AA')] <- AA

  isBB <- thrub==BB
  isBO <- thrub==BO
  thrub[ isBB & miscall( 'BB2BO')] <- BO
  thrub[ isBO & miscall( 'BO2BB')] <- BB

  if( check_genofreq) {
    # Deliberately indirect check, so that I'm not duplicating any ...
    # ... mistakes from the above code
    PUP6 <- snpg$locinfo$PUP6
    nam <- colnames( PUP6)

    # Stored with 21 columns compressing symmetric entries, so eg ...
    # ... OO/AB also covers AB/OO
    # Expand into full 36 entries...
    diffo <- substring( nam, 1, 2) != substring( nam, 4, 5)
    PUP6[ , diffo] <- 0.5 * PUP6[ ,diffo]
    revnam <- substring( nam, 4, 5) %&% '/' %&% substring( nam, 1, 2)
    PUP6_alt <- PUP6[,diffo]
    colnames( PUP6_alt) <- revnam[ diffo]
    PUP6 <- cbind( PUP6, PUP6_alt)
    nam <- colnames( PUP6)
    sum1 <- sum2 <- emp <-
        matrix( 0, n_loci, 6, dimnames=list( snpg$locinfo$Locus, genotypes6))

    for( ig in genotypes6) {
      sum1[ ,ig] <- rowSums( PUP6[ , substring( nam, 1, 2)==ig, drop=FALSE])
      sum2[ ,ig] <- rowSums( PUP6[ , substring( nam, 4, 5)==ig, drop=FALSE])
      emp[ ,ig] <- colMeans( thrub==ig)
    }

    thrub@ana <- sum1
    thrub@emp <- emp
  }

  callo <- sys.call()
  thrub@call <- callo
return( thrub)

Rd_doc <- r"---{
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kinference.R
\name{postprocess_simcheck_FSP_POP}
\alias{postprocess_simcheck_FSP_POP}
\alias{simcheck_FSP_POP}
\title{Check POP/FSP splitter; private.}
\usage{
postprocess_simcheck_FSP_POP(forp_sim, findings, plot. = TRUE)

simcheck_FSP_POP(snpg, chromos, N, check_genofreq = FALSE)
}
\arguments{
\item{forp_sim}{result of previous call to \code{simcheck_FSP_POP}}

\item{findings}{result of previous call to \code{find_FSPs_from_POPs_v2(
forp_sim, ..., keep_indiv=TRUE)}--- which generates the \code{is_same} and
\code{is_psex} matrices for (pairs * loci), and the by-locus vectors
\code{Pr_same_FSP}, \code{Pr_psex_FSP}, \code{Pr_same_POP},
\code{Pr_psex_POP}.}

\item{plot.}{whether to plot. Can't see why this would ever be false...}

\item{snpg}{a \code{snpgeno} object with known allele freqs and error rates
(i.e. \code{pbonzer} and \code{snerr} must exist in \code{snpg$locinfo})}

\item{chromos}{Controls the extent of linkage:how many "chromosomes" should
the loci be split between? There's no crossover within these "chromosomes",
but they inherit independently of each other. Loci are allocated to
chromosomes in sequence, so that if there are 3 chromosomes, the first locus
goes to "C1", the second to "C2", third to "C3", fourth to "C1", etc.
Returned object will have an extra field \code{$locinfo$chromosim}.}

\item{N}{number of pairs of each type (ie this many FSPs, and this many
POPs).}
}
\value{
\code{simcheck_FSP_POP} returns a \code{snpgeno} object with
\code{4*N} rows, and an attribute "callo" recording the call. The
\code{$info} contains just one column, "Our_sample", with names like "F17_A"
and "P6_B"; those would respectively indicate the first (A) member of the
17th FSP, and the second (B) member of the 6th POP. The \code{locinfo} field
is unchanged, except as noted under "chromos" above.
\code{postprocess_simcheck_FSP_POP} just plots.
}
\description{
This doco shouldn't be visible in the package! (KEYWORDS internal)
}
\details{
\code{simcheck_FSP_POP} simulates POPs and FSPs based on a known set of
loci, ready for checking \code{\link{find_FSPs_from_POPs}} or variants.
Quick.

\code{postprocess_simcheck_FSP_POP} can be used to graphically confirm that
the analytical probabilities of samenames and pseudoexclusion actually match
the simulated values--- it is amazingly difficult to get the formulae right.
You wouldn't generally need \code{postprocess...} if you are working with a
real set of loci.
}
\examples{

ForP <- simcheck_FSP_POP( my_snpg, 10, 1000)
ForP$locinfo$useN[] <- 3 # minimal but most robust genotype
test3 <- find_FSPs_from_POPs_v2( ForP,
    cbind( seq( 1, by=2, length=1000),
      seq( 1, by=2, length=1000)),
    keep_indiv=TRUE) # need keep_indiv=T for postprocess...() to work
postprocess_simcheck_FSP_POP( ForP, test3)
hist( test3, nc=50)
abline( v=test$E_FPstat, col='blue') # POP and FSP means
\dontrun{
# Sim from real loci--- in this case, with 'useN==4' for all loci since 6 looked a bit iffy
library( atease)
s11nodup_all4 <- s11nodup
s11nodup_all4$locinfo$use6 <- NULL
s11nodup_all4$locinfo$useN <- 4L
simbo4 <- simcheck_FSP_POP( s11nodup_all4, N=1000, chromo=20)
simbo4_next <- find_FSPs_from_POPs_v2( simbo4, cbind( seq( 1, 3999, by=2), seq( 2, 4000, by=2)), keep=T)
kinference:::postprocess_simcheck_FSP_POP( simbo4, simbo4_next) # looks OK; not needed by "user"
hist( simbo4_next$FPstat, nc=50)
abline( v=simbo4_next@E_FPstat, col='green') # theory means
# 95\% of all POPs should be Left of the line drawn next:
abline( v=simbo4_next@E_FPstat['POP']+2*sqrt( simbo4_next@V_FPstat), col='blue', lty=1)
# Can't say for FSPs, becoz linkage
# Real data:
testo4 <- find_FSPs_from_POPs( s11nodup_all4, pops_005)
abline( v=testo4$FPstat, col='red')
}

}
\seealso{
\code{\link{find_FSPs_from_POPs}}
}
\keyword{internal}


}---"

}


"simtest_Kstuff" <-
function( ck, n, nq=20) {
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




#' Find full-sib pairs among parent-offspring pairs or among half-sib pairs
#' 
#' For pairs already picked as likely parent-offspring pairs (POPs), i.e.,
#' those with a weighted pseudo-exclusion (WPSEX) statistic less than some
#' threshold, they might be full sibling pairs (FSPs). This function checks
#' potential POPs with very low WPSEX values for their potential to be FSPs.
#' 
#' The idea of \code{split_FSPs_from_POPs}--- though this is not the only
#' possible workflow--- is that pairs which are \emph{either} POPs \emph{or}
#' FSPs should stand out very clearly from everything else, via
#' \code{\link{find_POPs}}. Then the job is to pick between those
#' possibilities. The workflow is supposed to be:
#' 
#' \itemize{ \item nail POPs/FSPs first with \code{\link{find_POPs}} \item pick
#' between them with \code{split_FSPs_from_POPs} (update: this doesn't work
#' very well yet... use age info if you can) \item look for HSPs and filter out
#' already-known POPs and FSPs }
#' 
#' Hence the other function, \code{split_FSPs_from_HSPs}, is theoretically
#' unnecessary in that you have already run \code{\link{find_POPs}} and
#' \code{split_FSPs_from_POPs} so you should know which of your "HSPs" are
#' really something else. But, nevertheless it's handy to have.
#' 
#' However, an equally reasonable workflow might be:
#' 
#' \itemize{ \item nail HSPs and everything stronger with
#' \code{\link{find_HSPs}} \item split HSPs from POPs/FSPs with
#' \code{split_FSPs_from_HSPs} \item \code{split_FSPs_from_POPs} }
#' 
#' All \code{split_} functions return expected values under different possible
#' kin-types (not variances, since these cannot be predicted for all
#' kin-types). \subsection{Obsolete note The statistic for
#' \code{split_FSPs_from_POPs}- _which isn't as powerful as I'd hoped; I'm
#' going to redo it- is based on the weighted sum of the number of
#' exactly-matching 4-way genotypes, with weights chosen to have high power for
#' this particular discrimination. Weighting is optimized for the unlikely
#' scenario that POPs and FSPs are equally likely a priori, but in practice the
#' weights are not sensitive to this. The test is deliberately crude and
#' robust- e.g. it avoids exclusion-based checks- on the assumption that you
#' have enough loci to pick HSPs, so the more-related kin-types should be
#' slam-dunks. \bold{But} it doesn't seem powerful enough. More worked
#' needed...
#' 
#' \code{split_FSPs_from_HSPs} again uses 4-way genotypes only (to avoid having
#' to worry about errors) but in a properly optimal PLOD designed for FSP/HSP
#' discrimination- its expectation is positive for FSPs and negative for HSPs.
#' Theoretical means for those are returned as attributes (variances cannot be
#' predicted). Haven't added means for POPs or UPs since you're not "supposed"
#' to have those in the mix by the time you run \code{split_FSPs_from_HSPs},
#' but maybe I should fix that at some point. }
#' 
#' @aliases split_FSPs_from_HSPs split_FSPs_from_POPs split_FSPs_from_candiHSPs
#' @param snpg a \code{snpgeno} object
#' @param candiHSPs,candiPOPs candidate kin-pairs- normally, a dataframe with
#' rows being pairs and columns \emph{i} and \emph{j} (and possibly others)
#' e.g. from find_POPs() or find_HSPs(). Can also be a 2-column matrix (each
#' row again one pair).
#' @param gerr genotyping error rate:0.01 would mean 1\%. There's no default so
#' it's up to you; e.g., looking at replicate samples is a good way to estimate
#' it. Results of \code{split_FSP_from_POPs} should not be very sensitive to
#' the value, however.
#' @param use_obsolete_version the original code for POPs-vs-FSPs was based on
#' a different non-likelihood-based statistic; see \bold{Obsolete note}. It
#' turned out to have low statistical power, but did not require specifying a
#' \code{gerr}. For replicability purposes, you can still run it by setting
#' this parameter to \code{TRUE}.
#' @keywords misc
#' @examples
#' 
#' # pops_or_fsps <- find_POPs( mysnpg, ...)
#' ## do histograms etc to find likely ones
#' # discro <- split_FSPs_from_POPs( mysngp, pops_or_fsps %where% (wpsex < 0.042), gerr=0.01)
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
#' # discro2 <- split_FSPs_from_HSPs(  mysngp, h_or_f %where% (PLOD > 55))
#' # hist( discro2$PLOD_FH, nc=20, col='grey') # HSPs "should" be < 0, FSPs > 0
#' # abline( v=discro2@E_HSP, col='orange')
#' # text( discro2@E_HSP, par( 'usr')[4], 'POP', col='orange', pos=1) # below
#' # abline( v=discro2@E_FSP, col='lightblue')
#' # text( discro2@E_FSP, par( 'usr')[4], 'FSP', col='lightblue', pos=1) # below
#' 
#' @export split_FSPs_from_HSPs
"split_FSPs_from_HSPs" <-
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

  ret@E_FSP <- EPLOD_FH_F
  ret@E_HSP <- EPLOD_FH_H

  ret@call <- sys.call()

  return(ret)
}


"split_FSPs_from_POPs" <-
function( 
    snpg, 
    candiPOPs, 
    gerr, 
    use_obsolete_version=FALSE
){
## For pairs already picked as first-order kin (fsp or pops), eg via plod(hsp,up) > eta
## don't need full pairwise screening (do post hoc on a few hundred candidates), hence all in r.
## Based almost entirely on code from split_fsps_from_hsps(), just changed at the end...
## Do.in.envir() stuff is cos this belongs in 'kinference' pkg, but I've not put it there yet

  if( use_obsolete_version){
    # Non-lglk, based on ppn IBS *genotypes*
return( OLD_split_FSPs_from_POPs( snpg, candiPOPs))
  }

  if( missing( gerr) || 
    !is.numeric( gerr) || 
    length( gerr) != 1 || 
    !is.finite( gerr) || 
    gerr < 0 || 
    gerr >= 1
  ){
stop( "must define a geno error rate (gerr). Shouldn't need to be accurate...")
  }

  define_genotypes()

  # candidates normally from 'find_HSPs'; or can be M*2 matrix of rows in snpg that are poss HSPs
  # if former, make latter

  if( candiPOPs %is.a% 'data.frame') {
    candiPOPs <- as.matrix( candiPOPs[ cq( i, j)])
  }
  sibg <- just_sibg <- snpg[ c( candiPOPs),]

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


  nsib <- nrow( candiPOPs)
  nloci <- ncol(sibg)

  # if only R were zero-indexed...
  kappa_fsp <- c(1/4, 1/2, 1/4)
  # kappa_hsp <- c(1/2, 1/2, 0)
  kappa_pop <- c( gerr, 1-gerr, 0)

  p12fsp <- p12pop <- matrix(NA, nloci, nsib)


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
      # P[g_1l g_2l | POP]
      p12pop[l, i] <- kappa_pop[1] * P_k0[l, g1[i, l], g2[i, l]] +
                      kappa_pop[2] * P_k1[l, g1[i, l], g2[i, l]] +
                      kappa_pop[3] * P_k2[l, g1[i, l], g2[i, l]]
    }
  })
  OD_FP <- p12fsp/p12pop
  LOD_FP <- log(OD_FP)
  PLOD_FP <- colSums(LOD_FP)

  # Expectations
  # Need sum-to-1 here, so either rewrite when vecless2 appears, or work with compressed forms...
  OPHSP <- exp( OLOD) * OPUP # Pr[gg|HSP] <- 0.5 * PUP4 + 0.5 * Pr[gg|kappa=1]
  P_k0 <- OPUP
  P_k1 <- 2*OPHSP - OPUP
  P_k1[ P_k1 < 0] <- 0 # rounding error
  P_k2 <- 0 * P_k0 # get the shape right
  P_k2[ , diag( mg)] <- sqrt( P_k0[ , diag( mg)]) # only the cases where g1==g2

  p12fspa <- kappa_fsp[1] * P_k0 +  kappa_fsp[2] * P_k1 + kappa_fsp[3] * P_k2
  p12popa <- kappa_pop[1] * P_k0 + kappa_pop[2] * P_k1 + kappa_pop[3] * P_k2

  EPLOD_FP_F <- sum(log(p12fspa/p12popa) * p12fspa)
  EPLOD_FP_P <- sum(log(p12fspa/p12popa) * p12popa)

  # format a return object
  ret <- data.frame(PLOD_FP = PLOD_FP,
                    i       = candiPOPs[,1],
                    j       = candiPOPs[,2])

## Comment from FSP/HSP version was: Next 2 are COMPLETELY WRONG !!!
## ??? Is that true..??? They look OK to me!

  ret@E_FSP <- EPLOD_FP_F
  ret@E_POP <- EPLOD_FP_P

  ret@call <- mvb.sys.call(1) # don't ask, just works

return( ret)
}


"split_FSPs_from_POPs_6way" <-
function( snpg, candiPOPs, SDwt_POP=0.5) { # shouldn't have default--- hardwire or leave!

  # Don't need full pairwise screening for FSPs (do post hoc on a few hundred
  # candidate POPs), hence all in R.

  define_genotypes()

  # Transform to 4way genotypes
  # based on code in find_duplicates
  # careful, since "factor level" of AB and OO is different in 4way vs 6way
  snpg <- snpg[ c( candiPOPs),]

  # Code from hsp_power()
  li <- snpg@locinfo
  li1 <- li[1,]

  temp0 <- with( li1, calc_g6probs_IBD0_scalar( pbonzer, snerr, record=TRUE))
  cg6p0 <- make_playback( calc_g6probs_IBD0_scalar, temp0)

  temp1 <- with( li1, calc_g6probs_IBD1_scalar( pbonzer, snerr, record=TRUE))
  cg6p1 <- make_playback( calc_g6probs_IBD1_scalar, temp1)

  temp2 <- with( li1, calc_g6probs_IBD2_scalar( pbonzer, snerr, record=TRUE))
  cg6p2 <- make_playback( calc_g6probs_IBD2_scalar, temp2)

  g6p0 <- with( li, cg6p0( pbonzer, snerr))
  g6p1 <- with( li, cg6p1( pbonzer, snerr))
  g6p2 <- with( li, cg6p2( pbonzer, snerr))

  off <- 1 # until vecless has arbitrary-base arrays
  Pr_same_given_k <- array( 0, c( n_loci, 3))
  Pr_same_given_k[ l, {off+0}] := g6p0[l,g,g2=g] %[g]% g6p0[l,g,g2=g] # dunno if this will work...
  Pr_same_given_k[,{off+1}] <-
  Pr_same_given_k[,{off+2}] <- 1

  Pr_nsame_FSP <- c( 1/4, 1/2, 1/4)
  Pr_same_FSP[ l]:= Pr_nsame_FSP[ k] %[k]% Pr_same_given_k[ l, k]
  Pr_same_POP[ l]:= Pr_same_given_k[ l, {off+1}]

  SD_FSP <- sqrt( Pr_same_FSP * (1-Pr_same_FSP))
  SD_POP <- sqrt( Pr_same_POP * (1-Pr_same_POP))
  SD_denom <- SDwt_POP * SD_POP + (1-SDwt_POP) * SD_FSP
  wt <- (Pr_same_FSP - Pr_same_POP) / sqr( SD_denom)

}


"split_FSPs_from_POPs_v2" <-
function( snpg, candiPOPs, keep_indiv=FALSE) {
###### ADD ALIAS TO DOCO for find_FSPs_from_POPs and ensure this is exported
###### ... which this function should replace
###### ... keep the old one as 'old_find_FSPs_from_POPs' for now...
###### ... mainly so I don't have rewrite doco before today's mtg
    ## used to have unused var 'max_med_mul'. Removed.

    ## Don't need full pairwise screening for FSPs (do post hoc on a few hundred
## candidate POPs), hence all in R.

  define_genotypes()

  # 'candiPOPs' normally from 'find_POPs'; or can be M*2 matrix of rows in snpg that are poss POPs
  # if former, make latter

  if( candiPOPs %is.a% 'data.frame') {
    candiPOPs <- as.matrix( candiPOPs[ cq( i, j)])
  }
  snpg <- snpg[ c( candiPOPs),]

  n_pairs <- nrow( candiPOPs)

  useN <- snpg@locinfo$useN
  if( is.null( useN)) {
      useN <- ifelse( snpg@locinfo$use6, 6, 4)
    }

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
  # so that index {off} means 0 ibd, {off+1} means 1, {off+2} means 2

  Pr_same_given_k <- array( 0, c( n_loci, 3))

  Pr_same_given_k[ useN != 6,{off+2}] <- 1 # no errors between the "big 4" (or 3)

  Pr_same_given_k[ useN==3, {off+0}] <- (
      sqr( 2*pA*pB) +
      sqr( sqr( pA) + 2*pA*pO) +
      sqr( sqr( 1-pA))
    )[ useN==3]

  # Shared copy first, then other copies ie AB/O means ...
  # ... A shared, 2nd copy B and other 2nd copy O
  # ... (XuY) means X or Y (u for union)
  Pr_same_given_k[ useN==3, {off+1}] <- (
      pA * (sqr( pB) + sqr(1-pB)) +             # AB/B, A(AuO)/(AuO)
      pB * (sqr( pA) + sqr(1-pA)) +             # BA/A, B(BuO)/(BuO)
      pO * (sqr( pA) + sqr(1-pA))               # O(OuB)/(OuB)
    )[ useN==3]


  Pr_same_given_k[ useN==4, {off+0}] <- (sqr( 2*pA*pB) +
      sqr( sqr( pA) + 2*pA*pO) +
      sqr( sqr( pB) + 2*pB*pO) +
      sqr( sqr( pO))
    )[ useN==4]

  Pr_same_given_k[ useN==4,{off+1}] <- (
      pA * (sqr( pB) + sqr(1-pB)) +             # AB/B, A(AuO)/(AuO)
      pB * (sqr( pA) + sqr(1-pA)) +             # BA/A, B(BuO)/(BuO)
      pO * (sqr( pO) + sqr(pA) + sqr( pB))      # OX/X
    )[ useN==4]

  # Pretty ugly for 6way...
  # ... save a bit of work by re-using some PUP calcs
  P6 <- snpg@locinfo$PUP6
  snerr <- snpg$locinfo$snerr

  same6 <- do.on( strsplit( colnames( P6), '/'), .[1]==.[2])
  Pr_same_given_k[ useN==6, {off+0}] <- (rowSums( P6[ , same6])
    )[ useN==6]
  Pr_same_given_k[ useN==6, {off+2}] <- (
      2*pA*pB +
      sqr( pA) * ( 1- 2*( 1 - snerr[ , 'AA2AO']) * snerr[, 'AA2AO']) +
      sqr( pB) * ( 1- 2*( 1 - snerr[ , 'BB2BO']) * snerr[, 'BB2BO']) +
      2*pA*pO * ( 1- 2*( 1 - snerr[ , 'AO2AA']) * snerr[, 'AO2AA']) +
      2*pB*pO * ( 1- 2*( 1 - snerr[ , 'BB2BO']) * snerr[, 'BB2BO']) +
      sqr( pO)
    )[ useN==6]
  Pr_same_given_k[ useN==6, {off+1}] <- (
      pA * (
        sqr( pB) +
        sqr( pA) * ( sqr( snerr[ , 'AA2AO']) + sqr( 1-snerr[, 'AA2AO'])) +
        sqr( pO) * ( sqr( snerr[ , 'AO2AA']) + sqr( 1-snerr[, 'AO2AA'])) +
        2*pA*pO*( (1-snerr[ , 'AA2AO']) * snerr[, 'AO2AA'] + snerr[ , 'AA2AO'] * (1-snerr[ , 'AO2AA']))) +
      pB * (
        sqr( pA) +
        sqr( pB) * ( sqr( snerr[ , 'BB2BO']) + sqr( 1-snerr[, 'BB2BO'])) +
        sqr( pO) * ( sqr( snerr[ , 'BO2BB']) + sqr( 1-snerr[, 'BO2BB'])) +
        2*pB*pO*( (1-snerr[ , 'BB2BO']) * snerr[, 'BO2BB'] + snerr[ , 'BB2BO'] * (1-snerr[ , 'BO2BB']))) +
      pO * (
        sqr( pA) * ( sqr( snerr[ , 'AA2AO']) + sqr( 1-snerr[, 'AA2AO'])) +
        sqr( pB) * ( sqr( snerr[ , 'BB2BO']) + sqr( 1-snerr[, 'BB2BO'])) +
        sqr( pO) )
    )[ useN==6]

  Pr_nibd_FSP <- c( 1/4, 1/2, 1/4)
  Pr_nibd_POP <- c( 0, 1, 0)

  Pr_same_FSP[ l]:= Pr_nibd_FSP[ k] %[k]% Pr_same_given_k[ l, k]
  Pr_same_POP[ l]:= Pr_nibd_POP[ k] %[k]% Pr_same_given_k[ l, k] # or Pr_same_given_k[ l, {off+1}]

  # From find_POPs: UP case then POP case
  Pr_psex_given_k <- Pr_same_given_k
  Pr_psex_given_k[] <- NA

  Pr_psex_given_k[ useN==3,{off+0}] <- 2 * (
      (2*pA*pO + pA*pA) * (2*pB*pO + pB*pB + pO*pO)  # AAO/BBOO
    )[ useN==3]
  Pr_psex_given_k[ useN==3,{off+1}] <- (pO * (2*pA*(pB+pO))
    )[ useN==3]
  Pr_psex_given_k[ useN==3,{off+2}] <- 0

  Pr_psex_given_k[ useN==4,{off+0}] <- 2 * (
      (2*pA*pO + pA*pA) * (2*pB*pO + pB*pB) +  # AAO/BBO
      2*pA*pB * sqr( pO) # AB/OO
    )[ useN==4]
  Pr_psex_given_k[ useN==4,{off+1}] <- (pO * (2*pA*pB)
    )[ useN==4]
  Pr_psex_given_k[ useN==4,{off+2}] <- 0

  # Cols of P6 are alphabetical, but by *2nd* geno before 1st !
  Pr_psex_given_k[ useN==6,{off}] <-
    rowSums( P6[,c( 'OO/AB', 'BB/AA', 'OO/AA', 'OO/BB', 'BO/AA', 'BB/AO')])[ useN==6]
  Pr_psex_given_k[ useN==6,{off+1}] <- 2*(
      0 +  # AB/OO
      pO * (pA*pB*snerr[,'AO2AA']*snerr[,'BO2BB']) +  # AA/BB O shared; 1 A, 1 B, both misclassed
      pO * (pO*pA*snerr[,'AO2AA']) +  # AA/OO O shared, 1 O & 1 A, A misclassed
      pO * (pO*pB*snerr[,'BO2BB']) +  # BB/OO O shared, 1 O & 1 B, B misclassed
      pO * (pB*pA*snerr[,'AO2AA']*(1-snerr[,'BO2BB'])) + # AA/BO, O shared, A misclassed, B not
      pO * (pA*pB*snerr[,'BO2BB']*(1-snerr[,'AO2AA'])) # BB/AO, O shared, B misclassed, A not
    )[useN==6]
  Pr_psex_given_k[ useN==6,{off+2}] <- 0 # Not possible AFAICS

  Pr_psex_FSP[ l]:= Pr_nibd_FSP[ k] %[k]% Pr_psex_given_k[ l, k]
  Pr_psex_POP[ l]:= Pr_nibd_POP[ k] %[k]% Pr_psex_given_k[ l, k] # Pr_psex_given_k[ l, {off+1}]

  # For brevity:
  px <- Pr_psex_FSP
  ps <- Pr_same_FSP

  # 2x2 eqn to solve:
  # [ px(1-px), -ps*px   ] = [ Dxx] * lambda
  # [ -ps*px,   ps(1-ps) ]   [ Dss]

  inv_DET <- 1 / (px * (1-px) * ps * (1-ps) - sqr( ps*px))
  inv_DET <- 1 / (ps*px*(1-ps-px))
  Dx <- Pr_psex_FSP - Pr_psex_POP
  Ds <- Pr_same_FSP - Pr_same_POP

  # 2x2 matrix inverse...
  ws <- (px*(1-px)*Ds + ps*px*Dx) * inv_DET
  wx <- (ps*px*Ds + ps*(1-ps)*Dx) * inv_DET

  # Check:
  # rbind( ps*(1-ps)*ws - ps*px*wx, Ds)
  # rbind( -ps*px*ws + px*(1-px)*wx, Dx)

  # Paranoia... only happens if ps==0 or px==0
  wx[ !is.finite( inv_DET)] <- 0
  ws[ !is.finite( inv_DET)] <- 0

  # Not-very-optimal cap on ws...
  #ws <- pmin( ws, max_med_mul * median( ws))

  # Scaling could go here, eg to ensure ws[stat|FSP]=n_loci ignoring linkage
  V_FSP <- px*sqr(wx) + ps*sqr(ws) - sqr( px*wx + ps*ws)
  rescalor <- sqrt( n_loci / sum( V_FSP))
  wx <- wx * rescalor
  ws <- ws * rescalor
  V_FSP_nolink <- px*sqr(wx) + ps*sqr(ws) - sqr( px*wx + ps*ws)
  V_POP <- (px-Dx)*sqr(wx) + (ps-Ds)*sqr(ws) - sqr( (px-Dx)*wx + (ps-Ds)*ws)

  E_FSP <- wx * px + ws * ps
  E_POP <- wx * Pr_psex_POP + ws * Pr_same_POP

  # g1 & g2 are 6way. Produce 4way and 3way equivs...
  # Yet to write 'recode_geno'...
  snpg6 <- snpg

  snpg4 <- snpg
  snpg4@diplos <- genotypes4_ambig
  snpg4[ snpg6==AO] <- AAO
  snpg4[ snpg6==AA] <- AAO
  snpg4[ snpg6==BO] <- BBO
  snpg4[ snpg6==BB] <- BBO
  snpg4[ snpg6==OO] <- OO # need to do OO & AB too, since codes are different in 4way vs 6way
  snpg4[ snpg6==AB] <- AB

  snpg3 <- snpg4
  snpg3@diplos <- genotypes3
  snpg3[ snpg4==AB] <- AB
  snpg3[ snpg4==AAO] <- AAO
  snpg3[ snpg4==BBO] <- BBOO
  snpg3[ snpg4==OO] <- BBOO

  g1_6 <- snpg6[ 1 %upto% n_pairs,]
  g2_6 <- snpg6[ n_pairs + (1 %upto% n_pairs),]

  g1_4 <- snpg4[ 1 %upto% n_pairs,]
  g2_4 <- snpg4[ n_pairs + (1 %upto% n_pairs),]

  g1_3 <- snpg3[ 1 %upto% n_pairs,]
  g2_3 <- snpg3[ n_pairs + (1 %upto% n_pairs),]

  is_same <- is_psex <- matrix( FALSE, n_pairs, n_loci)
  is_same[, useN==3] <- (g1_3==g2_3)[, useN==3]
  is_same[, useN==4] <- (g1_4==g2_4)[, useN==4]
  is_same[, useN==6] <- (g1_6==g2_6)[, useN==6]

  # With psex, the g1/g2 order can matter--- so, do it one way, then swap g1 & g2
  is_psex[, useN==3] <- (
      (g1_3==AAO & g2_3==BBOO) |
      (g2_3==AAO & g1_3==BBOO)
    )[ , useN==3]
  is_psex[, useN==4] <- (
      (g1_4==AB & g2_4==OO) |
      (g1_4==AAO & g2_4==BBO)
    )[ , useN==4]
  is_psex[, useN==4] <- is_psex[, useN==4] | (
      (g2_4==AB & g1_4==OO) |
      (g2_4==AAO & g1_4==BBO)
    )[ , useN==4]
  is_psex[, useN==6] <- (
      (g1_6==AB & g2_6==OO) |
      (g1_6==AA & g2_6==BB) |
      (g1_6==AA & g2_6==OO) |
      (g1_6==BB & g2_6==OO) |
      (g1_6==AA & g2_6==BO) |
      (g1_6==BB & g2_6==AO)
    )[, useN==6]
  is_psex[, useN==6] <- is_psex[, useN==6] |       (
      (g2_6==AB & g1_6==OO) |
      (g2_6==AA & g1_6==BB) |
      (g2_6==AA & g1_6==OO) |
      (g2_6==BB & g1_6==OO) |
      (g2_6==AA & g1_6==BO) |
      (g2_6==BB & g1_6==AO)
    )[ ,useN==6]

  stat[i]:= wx[l] %[l]% is_psex[ i, l] + ws[ l] %[l]% is_same[ i, l]

  ret <- data.frame(
      FPstat= stat,
      i = candiPOPs[,1],
      j = candiPOPs[,2]
    )

  ret@E_FPstat <- c( POP=sum( E_POP), FSP=sum( E_FSP))
  # No point in returning V_FSP, since the variance will depend on linkage
  ret@V_FPstat <- c( POP=sum( V_POP))


  if( keep_indiv) {
    attributes( ret) <- c( attributes( ret), returnList(
        is_same, is_psex,
        Pr_same_FSP, Pr_psex_FSP,
        Pr_same_POP, Pr_psex_POP))
  }
  ret@call <- sys.call()
return( ret)

#  Set up for real data--- in this case, with 'useN==4' for all loci since 6 looked a bit iffy
#  s11nodup_all4 <- s11nodup
#  s11nodup_all4$locinfo$use6 <- NULL
#  s11nodup_all4$locinfo$useN <- 4L
#  simbo4 <- simcheck_FSP_POP( s11nodup_all4, N=1000, chromo=20)
#  simbo4_next <- find_FSPs_from_POPs_v2( simbo4, cbind( seq( 1, 3999, by=2), seq( 2, 4000, by=2)), keep=T)
#  kinference:::postprocess_simcheck_FSP_POP( simbo4, simbo4_next) # looks OK; not needed by "user"
#  hist( simbo4_next$FPstat, nc=50)
#  abline( v=simbo4_next@E_FPstat, col='green') # theory means
#  95% of all POPs should be Left of the line drawn next:
#  abline( v=simbo4_next@E_FPstat['POP']+2*sqrt( simbo4_next@V_FPstat), col='blue', lty=1)
#  Can't say for FSPs, becoz linkage
#
#  Real data:
#  testo4 <- find_FSPs_from_POPs( s11nodup_all4, pops_005)
#  abline( v=testo4$FPstat, col='red')
#

}


"split_HSPs_from_HTPs" <-
function( snpg, candiHTPs) {
  # For pairs already picked as possible HSPs, they might be HTPs

  # Don't need full pairwise screening for FSPs (do post hoc on a few hundred
  # HSPs), hence all in R.

  define_genotypes()

  # HSPs normally from 'find_HSPs'; or can be M*2 matrix of rows in snpg that are poss HSPs
  # if former, make latter

  if( candiHTPs %is.a% 'data.frame') {
    candiHTPs <- as.matrix( candiHTPs[ cq( i, j)])
  }
  sibg <- just_sibg <- snpg[ c( candiHTPs),]

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

  nsib <- nrow( candiHTPs)
  nloci <- ncol(sibg)

  # if only R were zero-indexed...
##   kappa_fsp <- c(1/4, 1/2, 1/4)
    kappa_hsp <- c(1/2, 1/2, 0)
    kappa_htp <- c(3/4, 1/4, 0)


  p12hsp <- p12htp <- matrix(NA, nloci, nsib)

  # split sibg into g1 and g2 for the two parts of the pairs
  g1 <- sibg[1:nsib, ]
  g2 <- sibg[(nsib+1):(2*nsib), ]

  # for loop version of the code
  g1 <- as.character(g1)
  g2 <- as.character(g2)
  evalq(for(i in 1:nsib){
    for(l in 1:nloci){
      # P[g_1l g_2l | HSP]
      p12hsp[l, i] <- kappa_hsp[1] * P_k0[l, g1[i, l], g2[i, l]] +
                      kappa_hsp[2] * P_k1[l, g1[i, l], g2[i, l]] +
                      kappa_hsp[3] * P_k2[l, g1[i, l], g2[i, l]]
      # P[g_1l g_2l | HTP]
      p12htp[l, i] <- kappa_htp[1] * P_k0[l, g1[i, l], g2[i, l]] +
                      kappa_htp[2] * P_k1[l, g1[i, l], g2[i, l]] +
                      kappa_htp[3] * P_k2[l, g1[i, l], g2[i, l]]
    }
  })
  OD_ST <- p12hsp/p12htp
  LOD_ST <- log(OD_ST)
  PLOD_ST <- colSums(LOD_ST)

  # Expectations
  # Need sum-to-1 here, so either rewrite when vecless2 appears, or work with compressed forms...
  OPHTP <- exp( OLOD) * OPUP # Pr[gg|HTP] <- 0.5 * PUP4 + 0.5 * Pr[gg|kappa=1]
  P_k0 <- OPUP
  P_k1 <- 2*OPHTP - OPUP
  P_k1[ P_k1 < 0] <- 0 # rounding error
  P_k2 <- 0 * P_k0 # get the shape right
  P_k2[ , diag( mg)] <- sqrt( P_k0[ , diag( mg)]) # only the cases where g1==g2

  p12hspa <- kappa_hsp[1] * P_k0 +  kappa_hsp[2] * P_k1 + kappa_hsp[3] * P_k2
  p12htpa <- kappa_htp[1] * P_k0 + kappa_htp[2] * P_k1 + kappa_htp[3] * P_k2

  EPLOD_ST_HS <- sum(log(p12hspa/p12htpa) * p12hspa)
  EPLOD_ST_HT <- sum(log(p12hspa/p12htpa) * p12htpa)

  # format a return object
  ret <- data.frame(PLOD_ST = PLOD_ST,
                    i       = candiHTPs[,1],
                    j       = candiHTPs[,2])

  # Next 2 are COMPLETELY WRONG !!!
##  ret@E_HSP <- EPLOD_ST_F
##  ret@E_HTP <- EPLOD_ST_H

  ret@E_HSP <- EPLOD_ST_HS
  ret@E_HTP <- EPLOD_ST_HT

  ret@call <- sys.call()

  return(ret)
}




#' Predict variance of PLOD for HCPs and HTPs
#' 
#' Aim is to work out how much your putative half-sibling pairs (HSPs) might be
#' contaminated by half-thiatic pairs (HTPs) or half-cousin pairs (HCPs) (or,
#' theoretically, by more remote kin). HSP-selection is presumably based on the
#' pairwise PLODs for HSP:UP, taking all pairs where that PLOD exceeds some
#' threshold. Given the allele freqs, the \emph{mean} PLOD is predictable when
#' the truth is UP, HCP, HTP, or HSP. The variance is only predictable for UPs,
#' though, because linkage makes loci non-independent for kin. However, an
#' empirical variance can be estimated for HSPs based on the observed PLODs
#' above some safe threshold (to exclude weaker kin), typically the mean PLOD
#' when truth is HSP. Based on the empirical variance for HSPs and the
#' analytical variance for UPs, we basically know how much linkage there might
#' be, so we can predict the PLOD variances for the other kin-pair types. The
#' wrinkle is that those more-remote variances also depend somewhat on the
#' finer-scale organization of the genome, i_e. whether it's lots of
#' chromosomes with no crossover, or a few chromosomes with lots of crossover.
#' \code{var_PLOD_kin} therefore calculates two versions, one assuming the
#' genome is entirely made up of equal-sized chromosome with zero crossover,
#' and the other assuming the genome is a single chromosomes with crossover
#' according to a memoryless random process. The output (basically, two
#' variance estimates which ought to bound the true variance for the
#' "contaminating" kin-type of interest--- subject to statistical noise) can be
#' fed into \code{autopick_HSP_threshold} (qv) to do what it says.
#' 
#' The "per-locus LOD" (whose properties are stored in the columns \code{e0},
#' \code{e1}, \code{v0}, \code{v1} in \code{linfo}) is created by calling
#' \code{\link{hsp_power}} (qv). The normal use-case would be that you've done
#' so with \code{k=0.5}, so that the (P)LOD pertains to HSP::UP comparisons.
#' However, if you called it with \code{k=0.25} then the (P)LOD would be
#' designed for HTP::UP comparisons, and so on. In fact, you could even
#' hand-tweak the calculations to contain LODs for HTP::HSP comparisons, which
#' \emph{might} in principle improve the resolution (but you'd have to fiddle
#' manually; you could actually do it based on two calls to
#' \code{\link{hsp_power}}, one with \code{k=0.5} and one with \code{k=0.25},
#' and manipulating the results). The other calculations in this function are
#' "agnostic WRTO", ie not intrinsically dependendent on, the values of
#' \code{e0/e1/v0/v1}, so the rest of the calcs should just work.
#' 
#' It's assumed that lots of loci are being used, so that the mix of loci on
#' each "chromo", or the splatter of loci along the single "megachromo", always
#' matches the overall population, on law-of-large-numbers grounds.
#' 
#' Stuff like uncertainly in allele frequencies, and in the PLOD variance for
#' HSPs, needs to be accounted for externally, by repeatedly drawing from the
#' posteriors and re-calculating the PLODs and re-running this function.
#' 
#' If the variance estimates show really good separation between the kin-pair
#' types, then one could refine the "preliminary variance" step by reducing the
#' super-high threshold (and assuming a truncated-Normal distribution). This
#' might be worthwhile if the preliminary variance otherwise has to be based on
#' a very small number of no-brainer HSPs. The "logical conclusion" of That
#' Kind Of Thing is some kind of MLE involving estimating the population of
#' different types of kin, and we really don't want to go there for now (since
#' that should include the population dynamics shebang). In other words, we'd
#' end up linking the genetic kin-finding model to the population dynamics
#' model, which makes life statistically harder. And god knows it's hard
#' enough. Anyway, if we were taking that approach, it might well be better to
#' avoid PLODs altogether and instead go for inferences about the actual ppn of
#' co-inherited loci, from which estimates-of-co-inherited-variance and
#' inferences about kin-ppns can be made. \subsection{Subtypes of kinGiven the
#' loci and the crossover rates, the PLOD variance for different kin-types is
#' mainly determined by the number of meioses. However, at least for the
#' with-crossover version, there is also \emph{some} effect of the \emph{type}
#' of kin within a given order: GGPs and HSPs would have \emph{slightly}
#' different variances. For CKMR purposes, the commonest type of kin of given
#' order are those born closest in time, so the algorithm always uses the type
#' with \emph{single} shared ancestor and minimax number of generations since
#' shared ancestor. This means HSPs for \code{n_meio=2} (FTPs have _two shared
#' ancestors; GGPs entail 2 generations of gap whereas HSPs have only 1), HTPs
#' for \code{n_meio=3}, HC1Ps for \code{n_meio=4}, etc. If you really wanted to
#' look at that, you could use the \code{V_allX} code inside this function,
#' which takes two arguments \code{short} and \code{long} for the length of
#' chains since the shared ancestor: for HSP, these are both 1, but for GGPs,
#' one is 0 and the other is 2. But since the single-chromo
#' equal-linkage-distance model is highly approximate anyway, do you really
#' care? }
#' 
#' @param linfo either a \code{snpgeno} object, or its "locinfo" attribute (or
#' a fake one). The "locinfo" should be a dataframe with columns \code{e0},
#' \code{e1}, \code{v0}, \code{v1}, \code{count}. Each row is one "type" of
#' locus, i_e., with roughly the same values of e/v 0/1, and \code{count} says
#' how many such loci there are. e/v 0/1 are means and variances of the
#' per-locus LOD (note no P) when the locus is or isn't co-inherited. See
#' \bold{Details}.
#' @param emp_V_HSP empirical variance of PLOD for deffo HSPs. You're supposed
#' to be running this on real data, so that \code{emp_V_HSP} is an actual
#' number; however, for testing purposes, you can set up an artificial version
#' via \code{C_equiv} below.
#' @param n_meio Target number of meioses:2 for 2nd-order kin (e.g. HSPs; also
#' GGPs and FTPs), 3 for 3rd-order (e.g. HTPs), etc. This is by far the main
#' driver of variance, but technically the only one; see \bold{Subtypes of
#' kin}.
#' @param debug Logical flag. Defaults to FALSE.
#' @param C_equiv for artificial test, with \code{emp_V_HSP} set to the
#' no-crossover variance from \code{C_equiv} chromos (need not be integer).
#' Ignored if \code{emp_V_HSP} is set.
#' @return Matrix with two rows \code{V0} and \code{Vx}, and one column for
#' each element of \code{n_meio} (which is always augmented to include 2),
#' named "M2" etc. The two rows pertain respectively to the no-crossover
#' multiple-chromosome scenario, and the single-chromosome multiple-crossover
#' scenario. The matrix also has an attribute \code{info}, which is a numeric
#' vector of elements named \code{V_UP}, \code{V_HSP}, \code{C_hat}, and
#' \code{rho_hat} (note that \code{V_HSP} should duplicate the first column of
#' the matrix).XXX>'. \code{C_hat} is estimated equivalent number of
#' chromosomes for the no-crossover scenario, and \code{rho_hat} is per-locus
#' crossover rate for the all-crossover scenario.
#' @keywords misc
#' @examples
#' 
#' # COMPLETELY MADE-UP e/v values! Nothing to do with genetics :)
#' var_PLOD_kin( data.frame( count=45, ev01= I( cbind( e0=-1, e1=2, v0=0.03, v1=0.02))), C_equiv=22, n_meio=3:4)
#' #        M2    M3    M4
#' #  V0 208.2 156.6  91.9
#' #  Vx 208.2 186.2 101.8
#' #  attr(,"info")
#' #      V_UP    V_HSP    C_hat  rho_hat
#' #    1.3500 208.2273  22.0000   0.2616
#' 
#' @export var_PLOD_kin
"var_PLOD_kin" <-
function(
    linfo,
    emp_V_HSP= V_noX( C_equiv, 2),
    n_meio,
    debug=FALSE,
    C_equiv=NULL
){
#########
stopifnot( # user != 'bozo',
    n_meio==floor( n_meio),
    all( n_meio >= 2)
  )

  n_meio <- as.integer( sort( unique( c( 2, n_meio))))

  if( linfo %is.a% 'snpgeno') {
    linfo <- linfo@locinfo
  }

  # linfo should be a DF with element ev01 having cols e0, e1, v0, v1
  # as produced by hsp_power
  # Becos of make_dull() [as of mvbutils 2.8.437] 'x %is.a% "matrix"' does NOT work
  # after make_dull( x); the implicit c( 'matrix', 'array') for class( unclass( x)) disappears
  # under make_dull() . Sigh. This is kind-of an R "feature" re implicit S3 classes...
  # But, is.matrix() does work

  THINGS <- cq(e0, e1, v0, v1)
stopifnot(
    is.matrix( linfo$ev01),
    ncol( linfo$ev01)==4,
    all( THINGS %in% colnames( linfo$ev01))
  )

  count_l <- linfo$count # only present if linfo is simulated eg 100 loci like this one, 100 like the next...
  if( is.null( count_l)){
    count_l <- rep( 1, nrow( linfo))
  }

  L <- sum( count_l)
  pi <- count_l / L # sum(pi)==1

  # Extract columns of ev01. Slightly odd code--- maybe don't need to bother
  ev01 <- unclass( linfo$ev01) # get rid of "dull" attr
  for( thing in THINGS ) {
    assign(thing %&% "_l", linfo$ev01[,thing]) # the "_l" is a pseudo subscript...
  }

  # What are the "typical" properties of a locus?
  e0 <- pi %**% e0_l
  e1 <- pi %**% e1_l
  v0 <- pi %**% v0_l
  v1 <- pi %**% v1_l

  # noX = no crossover, i_e. all in separate equal chromos
  V_noX <- function( C, meioses) {
      p <- 2 ^ (1-meioses) # Marginal prob of coinheritance = 1/2 for HSPs, 1/8 for HCPs
      EofV <- L * v0 + L * p *(v1-v0)
      VofE <- sqr( L * (e1-e0)) * p * (1-p) / C
    return( EofV + VofE)
  }

  # MoM for HSPs:
  if( !is.null( C_equiv)) {
    C_hat <- C_equiv <- min( L, C_equiv)
  } else {
    C_hat <- if( V_noX( 1, 2) < emp_V_HSP)
        1
      else if( V_noX( L, 2) > emp_V_HSP)
        L
      else
        find.root( V_noX, target=emp_V_HSP, start=1, step=1, fdirection='decreasing',
            min.x=1, max.x=L, meioses=2)
  }
  V0 <- do.on( n_meio, V_noX( C_hat, meioses=.))

  # Entirely Xover on one single chromo
  # Should these be done separately by locus type, then averaged??
  dee <- 1 %upto% (L-1)
  e2_1 <- pi %**% (sqr( e1_l) + v1_l)
  e2_0 <- pi %**% (sqr( e0_l) + v0_l)
  e1_0 <- e0 # pi %*% e0_l
  e1_1 <- e1

  V_allX <- function( rho, meioses) {
    # Assumes we are HSP or descendent ie stepladder; not extension ladder a la GGP

    pHSP_11 <- (1/4) * (1+exp( -4*rho*dee))

    # Extrapolate to number of meioses, m rather than HSP
    pm_11 <- ((0.25 * (1 + exp( -2*rho*dee))) ^ (meioses-2L)) * pHSP_11
    pm <- 2^(1L-meioses) # AKA pm_1: marginal prob of coin at a single locus

    pm_10 <- pm - pm_11
    # pm_01 <- pm_10 is just symmetry
    pm_00 <- 1 - pm_11 - 2*pm_10

    EL <- L * ( e1_1*pm + e1_0*(1-pm))
    EL2 <- L * ( e2_1*pm + e2_0*(1-pm)) +
        2 * (L-dee) %**% ( sqr( e1_1)*pm_11 + 2*e1_1*e1_0*pm_10 +  sqr( e1_0)*pm_00)
    VL <- EL2 - sqr( EL)
  }

  if( debug) {
    mtrace( V_allX) # surely wanna
    0 # help debugging...
  }

  # If no Odis (e_g. with very few loci!) then no point in going further
  rho_hat <- if( C_hat > L-1) 100 else
      find.root( V_allX, target= emp_V_HSP, start= 1/L, step= 0.2/L,
          fdirection= 'decreasing', min.x= 0,
          meioses= 2L) # for HSP case
  Vx <- do.on( n_meio, V_allX( rho_hat, meioses= .)) # for target kinship

  stuff <- rbind( V0, Vx)
  colnames( stuff) <- sprintf( 'M%i', n_meio)
  
  stuff@info <- unlist( returnList( V_UP=L*v0, V_HSP=emp_V_HSP, C_hat, rho_hat)) 
  # ... though V_HSP will be in V0 & Vx anyway

return( stuff)
}

