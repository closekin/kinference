# This is package kinference

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

  # As of v1.2, next should prolly be TRUE...
  assign( 'are_we_deprecating_yet', FALSE, envir=asNamespace( 'kinference'))

  # Zap attached plain-text docus, if any
  dedoc_namespace( pkgname)
}


".onUnload" <-
function( libpath){
  library.dynam.unload( 'kinference', libpath)
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


"autopick_threshold" <-
structure( function(
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

  E_HSP <- sum( li$E_HSP)
  E_UP <- sum( li$E_UP)

  E2 <- E_HSP
  E3 <- (E_UP + E_HSP) / 2 # 0 is good-enuf approx *IF* stat...
  # ...is actually HSP::UP PLOD, but it _might_ be something else eg HSP::HTP
  E4 <- (3*E_UP + 1*E_HSP) / 4

  if( E3 > min_PLOD){
    warning( sprintf(
        'fitrange_PLOD goes below 3rd-order mean, which is %5.2f; ' %&%
        'probably a bad idea', E3))
  }
  if( E2 > max_PLOD){
stop( sprintf( 'fitrange_PLOD exceeds 2nd-order mean, which is %5.2f;' %&%
      'noooo!', E2))
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

    # We observe Nobs kin that are either 2nd or 3rd (or 4th) order
    # but only if they're within fitrange_PLOD. So the total number of 2nd+3rd(+4th) would be...
    Nall <- Nobs / Pr_in_range
    N3 <- Nall * P3 # just 3rds

    # Threshold is where "FPtol_pairs" of 3rds would be expected
    # above it; 4ths irrel
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
    newvals3 <- returnList( SD2, SD3, SD4=0.1, P2=1-P3, P3, P4=0,
        Nall, lglk, thresh, Pr_FNeg)
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
      newvals4 <- returnList( SD2, SD3, SD4, P2, P3, P4,
          N3, Nall, lglk, thresh, Pr_FNeg)
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
    old_palette <- palette()
    on.exit( palette( old_palette)) # rude to just force it!
    kincols <- kinPalette( setPalette=TRUE)

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
      abline( v=thresh, lty=2, lwd = 2, col='black')
      probblies <- Nall * diff( CDF( histo$breaks,
          c( P2, P3, P4), c( SD2, SD3, SD4)))
      lines( histo$mids, probblies, col='black', lty=1, lwd = 3)

      # Might as well see the components...
      kinship_order <- c( POP=1, HSP=8, HTP=6, HCP=3)
      for( ord in c( 2, 3, if( use4th) 4 else NULL)){
        Eo <- get( 'E' %&% ord)
        SDo <- get( 'SD' %&% ord)
        Po <- get( 'P' %&% ord)

        probbly <- diff( pnorm( histo$breaks, mean=Eo, sd=SDo)) * Nall * Po
        lines( histo$mids, probbly, col=kincols[ kinship_order[ ord]],
            lty=2, lwd = 2)
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
    rect( par( 'usr')[1], 0, min_PLOD, par( 'usr')[4],
        density=DENSITY, border=NA, col= COL)

    # Means: show them last, so not covered by shading
    abline( v=c( E2, E3, E4), col= kincols[ cq( HSP, HTP, HCP)],
        lty = 2, lwd = 2)
    # ... E4 *shouldn't* appear; should be off LHS!

    legend("topright", legend = c( 'Overall', 'Threshold',
        'HSP', 'HTP', if(use4th) 'HCP' else NULL),
      lwd= c( 3, 2, 2, 2, if(use4th) 2 else NULL ),
      lty= c( 1, 2, 2, 2, if(use4th) 2 else NULL ),
      col= c( 'black', 'black', kincols[ 'HSP'], kincols[ 'HTP'],
          if( use4th) kincols[ 'HCP'] else NULL ),
      bg='white')
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
, doc =  mvbutils::docattr( r"{
autopick_threshold      package:kinference

PLOD threshold for HSPs


DESCRIPTION

This function proposes a PLOD threshold for excluding almost all 3KPs and 4KPs from true HSPs, and computes the associated false-negative probability (i.e., that a true HSP will have a PLOD below that threshold).

The modus operandi is to fit a mixture distribution to observed PLODs within some limited range that is expected to contain only 2KPs, 3KPs, and _perhaps_ a few 4KPs (and specifically excluding 1KPs). It is implicitly assumed that there _is_ a clear 2KP bump, and that it peaks close to the theoretical mean; if not, there is no point proceeding until you find more 2KPs (and perhaps tackle any QC issues). The threshold is then "chosen" (or proposed; the actual choice is really up to you) so that the expected number of false-positives from 3KPs plus 4KPs (i.e., with PLODs above the threshold) matches whatever you decide. A histogram with expected values is plotted (unless you tell it not to).

The mixture distribution treats the expected PLOD values for true 2KPs, 3KPs, and (if used) 4KPs as known (since those expectations can be calculated just from allele frequencies). However, because of linkage, the corresponding variances cannot be predicted a priori. Instead, the variance of the 2KP bump is estimated empirically from those PLODs which are above the 2KP mean, on the assumption that there will be few if any 3KPs with PLODs that large (hence the need for an upper limit on the included PLODs). Then, 'var_PLOD_kin' (qv) is used to infer the extent of linkage and to predict what the 3KP and 4KP variances might be, under two extreme assumptions about the detailed nature of linkage. The two assumptions lead to upper and lower bounds on the variances of 3KPs and 4KPs. Next, a Normal mixture distribution is fitted by ML, to estimate the height and variance of the 3KP bump, subject to those bounds; it can optionally allow for the presence of 4KPs as well. Once those parameters have been estimated, 'autopick_threshold' calculates the threshold PLOD above which 'FPtol_pairs' of 3KPs (and perhaps 4KPs) would be expected. Finally, it calculates the proportion of 2KP PLODs that would be expected to fall below that threshold.

Since it is better to make the threshold too high rather than too low, you can use the 'selecto' option to choose how the 3KP variance is estimated, subject to its bounds. The options are (i) strict ML estimation, to give the best fit to the overall PLODs, or (ii) choosing whichever 3KP variance leads to the _highest_ threshold (but still fitting the 3KP bump height by ML conditional on that variance).

The calculations in 'var_PLOD_kin' are in fact geared specifically to HSPs, HTPs, and HC1Ps, as opposed to the other subtypes of 2/3/4KPs. Accordingly, it is a good idea to restrict the pairs to those which, if the are indeed 2KPs, are a priori likely to be HSP rather than GGP or HTP. (If they are not 2KPs, then it probably matters less exactly what subtype they are.) Having said that, the sky will _probably_ not fall in if some GGPs are also included; in some circumstances it might even be worth doing so deliberately, in order to increase the 2KP sample size. But the general principle of kin-finding applies: it is best to concentrate on the comparisons that yield the outcomes you're interested in, rather than including lots of extraneous comparisons and then trying to filter the results.

Despite the name 'auto...', _you_ still have to supply sensible values for a couple of parameters, based on looking at your data and understanding what you are trying to do. So it's not _completely_ automated process, and it never will be! Choosing a threshold is _not_ an "optimization process" with explicit bias/variance tradeoffs; rather, it's about ensuring that you have adequate "engineering tolerance" in the next stage of CKMR. The False-Negative Probability will, to a great extent, compensate for the choice of threshold (i.e. removing any bias in the fitted CKMR model) _unless_ you set the threshold too low, and thus end up with some 3KPs in your set of "definite 2KPs".


.FITRANGE.AND.USE4TH

Assuming that you are using a 2KP-oriented PLOD, then the range of PLODs you fit to should extend from somewhere above 0 (which is very close to the expected PLOD for 3KPs), up to the RHS of the 2KP bump, but not so large as to include any FSPs and POPs. You can use the 'fitrange_PLOD' parameter to control this, instead of subsetting the 'kin' argument manually. If you set the lower limit high enough, then you don't need to worry about 4KPs intruding (their contribution would be negligible), so you can get away with fitting a 2-component mixture (simpler, less to go wrong...) by setting 'use4th=FALSE'. But if you push the lower range closer to 0 (which does give you a larger sample size for fitting), then you might need to set 'use4th=TRUE'. The (substantial) downside of doing that, is that there are often more PLODs close to 0 than near the 2KP mean, so the mixture-fit (which has to make more assumptions when also using 4KPs, and those assumptions may not be perfect) will concentrate its efforts on getting a good fit near 0, rather than near the 2KP mean which is what really matters. It is worth experimenting.


USAGE

autopick_threshold(
  x,
  kin,
  fitrange_PLOD,
  FPtol_pairs,
  use4th,
  selecto = c("ML", "paranoid"),
  NVAR = 10,
  plot_bins = NULL,
  shading_density = 10,
  want_all_results = FALSE,
  ...
)


ARGUMENTS

  x: a 'snpgeno' or its 'locinfo' attribute. Must already have been prepared by running 'kin_power'.

  kin: a dataframe of "close-ish" kin-pairs and their PLODs, presumably from running 'find_HSPs'; must have a column "PLOD".

  fitrange_PLOD: two numbers, specifying the range of PLODs from 'kin' to use in _fitting_ (though all are _plotted_, by default)

  FPtol_pairs: how many expected False-Positive 3KPs (and 4KPs, if 'use4th=TRUE') should the threshold exclude?

  use4th: whether to allow for 4KPs when fitting.

  selecto: whether to choose the threshold based on the best mixture fit ("ML"), or the most conservative ("paranoid").

  NVAR: how many variances to try, between the limits set by 'var_PLOD_kin'

  plot_bins: bin-width for histogram plotting. Default NULL means no plot.

  shading_density: By default, all PLODs in 'kin' will be included in the histogram, even though only a subset are used in fitting. The histogram bars _not_ used in fitting (i.e., below 'fitrange_PLOD[1]') will be lightened in colour, according to this parameter. Results are graphics-device-dependent, so you may need to experiment away from the default; larger numbers usually mean lighter shading. Setting 'shading_density=NA' should result in a light transparent rectangle covering the entire LHS of the graph, which you might prefer. You can also set 'xlim' as usual, to remove those left-hand bars altogether.

  want_all_results: if TRUE, return dataframe(s) containing results for each variance explored. This lets you examine "sensitivity".

  ...: other parameters passed to 'hist', eg 'xlim', 'ylim', 'col'. Many others will be ignored, and some will cause problems.


VALUE

The proposed threshold, with lots of attributes. You can use those to calculate false-neg probabilities for _other_ possible thresholds, as per *Examples*; you _don't_ have to accept the one that is proposed here! Threshold choice is _up to you_ (and badly-chosen thresholds are not the fault of 'kinference')!


SEE.ALSO

'kin_power', 'var_PLOD_kin'


EXAMPLES

library( atease) # for convenient 'thresh@info' below
dropbears1 <- kin_power( dropbears, k = 0.5)
hsps <- find_HSPs( dropbears1, keep_thresh = -10)
histoPLOD( hsps, log=TRUE) ## observe HSP - POP gap centred ~ PLOD = 120.
## Use that air-gap to set fitrange_PLOD.
thresh <- autopick_threshold( x= dropbears1, kin= hsps,
  fitrange_PLOD= c(0, 120), FPtol_pairs= 1, use4th= TRUE, plot_bins= 5)
thresh ## the 2nd order / 3rd order threshold value
thresh@info["Pr_FNeg"] ## ...
## ... the estimated 2KP false-neg rate, given that threshold
}")

)

"basic_sanity_checks_pairfinding" <-
function( check_geno_encoding=TRUE,
    nlocal=sys.parent()) mlocal({
## Called by all find_<blah>() functions
stopifnot(
    snpg %is.a% 'snpgeno',
    is.numeric( subset1),
    is.numeric( subset2),
    all( !duplicated( subset1)),
    all( !duplicated( subset2)),
    my.all.equal( subset1, subset2) ||
        !length( intersect( subset1, subset2)),
    all( c( subset1, subset2) %in.range% (1:nrow( snpg))),
    ij_numeric || !is.null( rowid_field( snpg)),
    # Most find_<blah> requires g6 or g4ambig
    # but that might change
    !check_geno_encoding || my.all.equal( snpg@diplos, genotypes6) ||
        my.all.equal( snpg@diplos, genotypes4_ambig)
  )
})


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


"calc_kinPower_checksum" <-
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
# Not sure what this was meant for. It's similar to code in other places...

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

return( returnList(
  pIBD0, pIBD1, pIBD2
))
}


"chain_pairwise" <-
structure( function( thing) {
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
, doc =  mvbutils::docattr( r"{
chain_pairwise      package:kinference
get_chain

Sib-groups within HSPs


DESCRIPTION

For checking veracity of _potential_ half-sibs or other kin-pairs. 'chain_pairwise' organizes them into chains within which each sample can be linked to another by a succession of direct pairwise links. The general idea is that real HSPs will be in clusters of 2 or 3; a spurious sample with a "lucky" genotype that wants to be everybody's mate will appear in a big but incomplete chain of mostly false-positives, where the direct pairwise links between the other chain-members are weak. You'd only run 'chain_pairwise' for pairs with a PLOD (or whatever statistic is being used) within a particular suspect range, so each chain may have false-negatives (i.e. missing direct links), but the general idea should be clear.

'chain_pairwise' and 'get_chain' are diagnostic tools for when 'kinference' goes wrong. For a comprehensive approach to complicated links between _genuine_ sib-groups (full and half), e.g. within larval samples such as for Atlantic Bluefin Tuna, see the 'sibgrouper' package (if not on R-universe, then contact the 'kinference' authors).


USAGE

chain_pairwise( thing)
get_chain( thing, seed)


ARGUMENTS

  thing: output from 'find_HSPs' or 'find_POPs' etc, or some subset thereof

  seed: one sample ID, interpreted as a row-number in 'thing'. To do:also allow names, via 'info' attr.


DETAILS

'get_chain' finds the chain for one specific sample.


VALUE

'chain_pairwise' returns a list of matrices, each for one chain; the rows and columns of each matrix are the samples in that chain. A "+" in the matrix indicates that those two samples have a direct pairwise link (i.e., they appear together in one row of 'thing'); a "." means not. The rows and columns of each matrix are sorted so that the linkiest samples are on the bottom and right. 'get_chain' returns the row-subset of 'thing' that is chained to 'seed'.
}")

)

"check6and4" <-
structure( function(
  geno6,
  thresh_pchisq_6and4,
  return_what=c( 'just_pvals', 'all'),
  extra_title = "",
  show6 = FALSE
){
##########
  n_fish <- nrow( geno6)
  n_loci <- ncol( geno6)
  define_genotypes()
  diplos <- geno6@diplos
    stopifnot(
        my.all.equal( sort( genotypes6), sort( diplos)) |
        my.all.equal( sort( genotypes4_ambig), sort( diplos))
    )
    if( my.all.equal( sort( genotypes4_ambig), sort( diplos))) {
      message("genotypes are stored 4-way. Skipping 6-way GOF tests...")
      do6 <- FALSE
      if( is.null(geno6@locinfo$snerr)) {
        snerrmat <- matrix(0, ncol(geno6), 4,
            dimnames = list(NULL, c("AA2AO", "AO2AA", "BB2BO", "BO2BB")))
        geno6@locinfo$snerr <- snerrmat
      }
    } else {
      message("genotypes are stored 6-way. Doing both 6-way and 4-way GOF tests...")
      do6 <- TRUE
    }
    if( do6 & !show6) {
      message( "calculating both 6-way and 4-way pvals, plotting 4-way...")
    }

  p6 <- with( geno6@locinfo,
      calc_g6probs( pbonzer[,'A'], pbonzer[,'B'], pbonzer[,'C'],
          snerr=snerr))[,genotypes6]

  p4 <- matrix( NaN, nrow(p6), 4)
  colnames(p4) <- genotypes4_ambig
  p4[,"OO"] <- p6[,"OO"]
  p4[,"AB"] <- p6[,"AB"]
  p4[,"AAO"] <- p6[,"AA"] + p6[,"AO"]
  p4[,"BBO"] <- p6[,"BB"] + p6[,"BO"]

  # At some point, will also need p6_IBD... later...

  # Chi-sq check: requires gpred & gobs attr on geno6
  g6pred <- p6 * n_fish
  g6obs <- matrix( 0, n_loci, 6, dimnames=list( NULL, genotypes6))

  if(do6) {
    for( ig in genotypes6) {
      g6obs[,ig] <- colSums( geno6==match( ig, diplos))
    }

    if(show6) {
      pval6 <- chisq_genofreq_check( geno6, gobs=g6obs, gpred=g6pred,
          test='G', showPlot = TRUE,
          thresh_pchisq_loci=thresh_pchisq_6and4, trim=FALSE,
          extra_title = extra_title)@locinfo$pval
    } else {
      pval6 <- chisq_genofreq_check( geno6, gobs=g6obs, gpred=g6pred,
          test='G', showPlot = FALSE,
          thresh_pchisq_loci=thresh_pchisq_6and4, trim=FALSE,
          extra_title = extra_title)@locinfo$pval
    }
    g4obs <- gtab6to4( g6obs)
    g4pred <- gtab6to4( g6pred)
  }

  if( !do6) {
    g4pred <- p4 * n_fish
    g4obs <- matrix( 0, n_loci, 4, dimnames=list( NULL, genotypes4_ambig))

    for( ig in genotypes4_ambig) { ## need to fix this for diplos != genotypes6
      g4obs[,ig] <- colSums( geno6==match( ig, diplos))
    }
  }

  pval4 <- chisq_genofreq_check( geno6, gobs=g4obs, gpred=g4pred,
      test='G', thresh_pchisq_loci=thresh_pchisq_6and4, trim=FALSE,
      extra_title = extra_title)@locinfo$pval

  return_what <- match.arg( return_what)
  if( return_what=='all') {
    if( do6) {
      geno6@locinfo$pval6 <- pval6
    }
    geno6@locinfo$pval4 <- pval4
return( geno6)
  } else {
    if( do6) {
return( returnList( pval6, pval4))
    } else {
return( returnList( pval6 = NULL, pval4))
    }
  }
}
, doc =  mvbutils::docattr( r"{
check6and4      package:kinference

Locus QC check


DESCRIPTION

Compares 6-way (if available) and 4-way genotype counts to HWE expectations, based on already-estimated allele frequencies. Plots a histogram of p-values across all loci, and (more importantly) plots observed / expected counts for each genotype by locus.

An overall p-value is calculated for each locus, based on chisq(1) of the G-statistic; it is not obvious how many DoF should be used. Anyway, the p-value itself is just indicative; few if any CKMR datasets actually "fit properly" because sample sizes are so large that even trivial misfits lead to significant departures from HWE null hypothesis, even though the final outcomes of kin-finding (ie PLODdograms) can look perfectly good. This is usually obviuos from the initial histogram of p-values, which in theory should be uniform for the "good" loci; in practice, it usually has heavy left-skew. Thus, the p-values are not to be taken literally, but just as a relative guide; the worst loci have the smallest p-values, and you can use the p-value as a criterion for dropping some loci. Since QC is iterative (samples, loci, samples, loci, ...) you can always revisit the decision later.

The p-value threshold ('thresh_pchisq_6and4') simply determines the colour used to plot each locus. You can manually changing the threshold and re-run until you have a visually satisfactory pattern of orange vs green. Green means good (or good enough); orange means bad, ie below the lower threshold. (You are actually allowed to supply two numbers for the threshold, in which case loci with in-between values will display pink; that is probably a design flaw, because it's confusing to have more than one threshold). Once you have identified a threshold that _looks good_, then you can use the 'pval4' or 'pval6' return-value to keep or discard loci, depending whether they are above or below that threshold.

The plots never look absolutely perfect, and there is no absolute criterion for "how bad is too bad". So, judgement and experience are required. But remember that keep-or-drop decisions aren't final; the whole QC process can be somewhat iterative, and the ultimate test is the PLODdogram(s) at the end. If it looks bad, ie bumps in the wrong places, then you may not have been restrictive enough (ie your threshold might be too high); whereas if the bumps are in the right places but are too wide for kin-finding, then your threshold might be too low.

With 6-way genotyping (eg SBTuna), calculations are done separately for the 6-way and 4-way versions. It is quite possible for a locus to look bad in the 6-way version, but good in the 4-way version; if so, don't just throw it out entirely, but try setting 'useN=4' for that locus, or 'useN=3' if null frequency is dangerously low.

For examples, see the 'kinference-vignette'.


USAGE

check6and4(
  geno6,
  thresh_pchisq_6and4,
  return_what= c("just_pvals", "all"),
  extra_title= "",
  show6= FALSE
)


ARGUMENTS

  geno6: a 'snpgeno' object with 4-way (or optionally 6-way) genotypes

  thresh_pchisq_6and4: a pair of thresholds for "bad" and "really bad" p-values. These determine the color in which each locus appears in all subplots.

  return_what: one of 'just_pvals' or 'all'; see value

  extra_title: a character string to be added to the bottom-right corner of all plots. Best if < 25 characters.

  show6: show the plots for 6-way goodness-of-fit? Defaults to TRUE. If diplos is anything other than genotypes6, should be FALSE.


VALUE

Creates per-locus vectors 'pval6' and 'pval4' for 6-way and 4-way genotypes respectively. If 'return_what="just_pvals"', these are returned in a list; if 'return_what="all"', they are added as columns to 'geno6$locinfo'.
}")

)

"chisq_genofreq_check" <-
structure( function( lociar,
    gpred= lociar@gpred,
    gobs= lociar@gobs,
    thresh_pchisq_loci,  # NULL to not worry; 1 value a threshold; 2 vals to inspect "iffy" ones
    showPlot = TRUE,
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

  if(showPlot) {
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

   ##   plot.new() ## should help with knitr, which currently pushes 6and4 plots off the page

    # Locussy fits
    opar <- par(mfrow=c(2, ncol( gpred) %/% 2),
        mar=c( 3, 3, 0, 0)+0.1, oma=c( 2, 2, 3, 1))
    # omi=c(0,0,.6,0),mai=c(.8,.8,.2,.2)) Paige uses absolute margins
    on.exit( par( opar))

    gtypes <- colnames( gobs)
    for( g in 1:ncol( gpred)) {
      # all of 'em
      # pch='.' is too small
      plot(gpred[,g], gobs[,g], pch=16, cex=0.4, col='lightblue',
          xlab='', ylab='', main='',
          xlim=c( 0, max( c( gpred[,g], gobs[,g]))),
          ylim=c( 0, max( c( gpred[,g], gobs[,g]))))
      abline(0,1,col=8,lwd=2)
      points( gpred[ keep_loci,g], gobs[ keep_loci,g], pch=16, cex=0.4, col='green') # overplot to make visible against the line
      mtext( side=3, gtypes[ g], line=-1)
      # bad ones with X
      points( gpred[ !keep_loci, g], gobs[ !keep_loci, g], col='magenta',
          pch=16, cex=0.6)
      if( any( iffy_loci)) {
          points( gpred[ iffy_loci, g], gobs[ iffy_loci, g], col='orange',
              pch=4, cex=1)
          points( gpred[ iffy_loci, g], gobs[ iffy_loci, g], col='orange',
              pch=16, cex=0.6)
          ## because otherwise we have a visible magenta shadow in the middle of each X
      }
    }
    mtext( 'Green = "good"', side=3, cex=1.5, outer=TRUE, line=1.5)
    mtext( 'Expected', side=1, cex=1.5, outer=TRUE)
    mtext( 'Observed', side=2, cex=1.5, outer=TRUE)
          mtext(extra_title, side = 1, adj = 1, padj = 4) ## SB

    if( liffies) {
      # MVB: the "ifs" below look dodgy. Surely there should be an 'else' instead of the comma..?
      legend( 'bottomright', pch=c( 16, rep( 4, liffies)), col=c( 'magenta', if( liffies>1) 'orange', 'green'), pt.cex=c( 0.6, if( liffies>1) 1, 0.6),
          legend=c( sprintf( '%4.1e < Pr ...', thresh_pchisq_loci[ 2]),
          if( liffies>1) sprintf( '... < %4.1e', thresh_pchisq_loci[1]), '... < 1') )
    }
  }

  if( trim) {
    lociar <- lociar[ , keep_loci, ,drop=FALSE]
  }
return( lociar)
}
, secret_doc =  mvbutils::docattr( r"{
chisq_genofreq_check      package:kinference

Check observed genotypes against HWE expectations


DESCRIPTION

Checks observed genotype frequencies against expected frequencies, presumably with expectation defined by HWE.


USAGE

chisq_genofreq_check(
  lociar,
  gpred = lociar@gpred,
  gobs = lociar@gobs,
  thresh_pchisq_loci,
  showPlot = TRUE,
  test,
  trim,
  seq_paxis = 0.025,
  extra_title = ""
)


ARGUMENTS

  lociar: a snpgeno object

  gpred: predicted allele frequencies

  gobs: observed allele frequencies

  thresh_pchisq_loci: a param. Presumably, a threshold p-val for flagging loci with suspicious-looking allele frequencies.

  test: a character string, either "Pearson" or "G"

  trim: TRUE or FALSE. TRUE will keep only above max thresh_pchisq_loci. Arguably better as to be done post-hoc.

  seq_paxis: numeric#'

  extra_title: a character string to be added to the bottom-right corner of all plots. Best if < 25 characters.


SEE.ALSO

kinference::check6and4
}")

)

"drop_dups_pairwise_equiv" <-
structure( function( ij, want_groups=FALSE) {
  if(ncol(ij) == 3) {
    ij <- ij[,2:3]
  }  ## so that users don't have to specify it every time

  if( !nrow( ij)){
return( integer())
  }

  ij <- as.matrix( ij) # in case it was a data.frame
  uij <- unique( c( ij))
  m <- nrow( ij)
  ij <- matrix( match( ij, uij), m, 2) # make sure it's numeric
  n <- max( ij)

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
}
, doc =  mvbutils::docattr( r"{
drop_dups_pairwise_equiv      package:kinference

Grouping duplicate samples


DESCRIPTION

'find_duplicates' only does pairwise comparisons. However, tissue from the same animal may turn up in multiple samples, so that one sample may turn up in many duplicate-pairs, and the pairs are linked. This function constructs equivalence classes- each corresponding notionally to one _animal_- showing which samples belong in each class. It can either return the entire set of classes, or it can pick just one sample from each class and then return the "surplus" duplicate samples. With the latter, if you then drop those elements as per 'vignette("kinference")', only one sample from each animal will be retained.

The algorithm merges two classes whenever any sample in the first class is flagged as a duplicate of any sample in the second. Thus, it will be sensitive to false-positive duplicates (though less so to false-negative ones). It's up to you to make sure that the input really contains true duplicates!


USAGE

drop_dups_pairwise_equiv(ij, want_groups = FALSE)


ARGUMENTS

  ij: 2-column matrix or data.frame; possibly row numbers in a dataset, or strings as per 'rowid_field' (qv).

  want_groups: if 'TRUE', also return the equivalence-classes themselves, as attribute 'groups'.


DETAILS

Input should be row numbers in a 'snpgeno' objects of duplicates, as a two-column data.frame or matrix with each row being a pair of duplicates, or the output from 'find_duplicates' (a 3-col matrix). Identifies 'groups' of equivalent observations (e.g., if i and j are duplicates, and j and k are duplicates, then i, j, and k are all equivalent). Outputs a vector of the row numbers for all-but-one of each group.


VALUE

Surplus elements in 'ij', perhaps plus attributes 'groups' if 'want_groups=TRUE'. You can look at that to figure out which elements are being retained (one "representative" from each equiv class). If 'ij' has no rows, an empty integer vector is returned, without any attributes.


SEE.ALSO

chain.pairwise


EXAMPLES

pairs <- matrix( c(
294, 289,
328, 294,
904, 857,
905, 904),
    ncol=2, byrow=TRUE)
drop_dups_pairwise_equiv( pairs, TRUE)
#[1] 289 328 857 905
#attr(,"groups")
#attr(,"groups")$`5`
#[1] 294 328 289
#
#attr(,"groups")$`6`
#[1] 904 905 857
}")

)

"DUP_paircomps_incomplete_lots" <-
function(geno1, geno2, symmo, max_diff_ppn, limit) {
    .Call(`_kinference_DUP_paircomps_incomplete_lots`, geno1, geno2, symmo, max_diff_ppn, limit)
}


"DUP_paircomps_lots" <-
function(geno1, geno2, symmo, max_diff_loci, keep_n, nbins, binterval, maxbin) {
    .Call(`_kinference_DUP_paircomps_lots`, geno1, geno2, symmo, max_diff_loci, keep_n, nbins, binterval, maxbin)
}


"est_ALF_6way" <-
structure( function( snpg, control=list()) {
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
, doc =  mvbutils::docattr( r"{
est_ALF_6way      package:kinference

Estimate ALFs from 6-way genotypes and snerr


DESCRIPTION

Performs 6-way re-estimation of ALFs, given 6-way genotypes, starting 4-way estimates of ALFs (from est_ALF_ABO_quick), and estimates of "snerr" (single-to-null error rates for apparent homozgyotes, per locus). Used as a second-pass estimate after ALFs have already been estimated roughly from 4-way genotypes (hence the phrase "re-estimation"). Used in CSIRO pipelines, e.g. for SBTuna, where the genotyping and estimation of snerr has already been done by routines in the 'genocalldart' package. If this makes no sense to you, then just stay away from 6-way genotyping!


USAGE

est_ALF_6way(snpg, control = list())


ARGUMENTS

  snpg: a 'snpgeno' object with 6-way genotypes (i.e., 'diplos(snpg)' matches 'get_genotype_enconding()$genotypes6'), with 'snerr' and 'pbonzer' included

  control: as per 'nlminb'


SEE.ALSO

'est_ALF_ABCO', 're_est_ALF', and 'est_ALF_ABO_quick'


EXAMPLES

head( bluefin$locinfo$snerr) ## has to exist for 6-way genotypes
bluefin$locinfo$pbonzer <- NULL ## remove pre-existing allele freq
## estimates!
bluefin <- est_ALF_ABO_quick( bluefin)
head( bluefin$locinfo$pbonzer)
bluefin <- est_ALF_6way( bluefin)
head( bluefin$locinfo$pbonzer)
}")

)

"est_ALF_ABCO" <-
structure( function( lociar, geno_amb = attr( lociar, 'geno_amb')) {
########## Taken largely from "pipeline_for_SBT_baits.r"
########## MVB: I'd like to clean this up
########## Careful "parallel Newton-Raphson" could allow vectorization and whoosh-factor, but NFN I guess

  define_genotypes() # AAO etc
  ## geno_amb <- lociar@geno_amb

stopifnot(
    !is.null( geno_amb),
    is.character( geno_amb),
    my.all.equal( geno_amb@diplos, genotypes_ambig))

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
  gobs <- gpred <- matrix( 0, n_loci, length( genotypes_ambig),
      dimnames=list( NULL, genotypes_ambig))
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
    pC <- if( !has_C) 0 else
        1 - sqrt( 1- sum( gobs[ll, cq( CCO, AC, BC)]) / (1 + n_fish))

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

    start_par <- c( logit( pA), logit( pB / (1-pA)),
        if( has_C) logit( pC / (1-pA-pB)))

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
, doc =  mvbutils::docattr( r"{
est_ALF_ABCO      package:kinference
re_est_ALF

Estimate allele frequencies olde-style, including nulls


DESCRIPTION

Most users should avoid these! 're_est_ALF' calculates maximum-likelhood estimates of minor, null, and 3rd-allele freqs from a 'snpgeno' object. Its workhorse is 'est_ALF_ABCO' which has a weirder syntax unless you are using an old CSIRO "pipeline", and is documented here for that reason only.

're_est_ALF' will accept genotypes 4-way genotypes including "AAO" and "BBO", 6-way genotypes (but there are better options in that case; see below), or triallelic "ABCO" genotypes. The latter corresponds to 'genotypes_ambig' as seen in the code of 'define_genotypes'; they allow for an optional 3rd allele and nulls, but do not distinguish between single-nulls and homozygotes. Historically, for CSIRO users, "ABCO"-type genotypes are produced by 'genocalldart::geno_deambig_ABC'.

Null-allele frequency has to be estimated from HWE deviations, so good estimates require a decent sample size.

The original use-case for these functions was datasets where bona fide 3rd alleles are common; even though they are not used in any 'kinference' step (because they get recoded to nulls, i.e. "neither A nor B" en route), it's useful to have them around for ALF estimation. There's a bit of history behind this, which I shan't go into here.

If there are no 3rd alleles, then 're_est_ALF' is not the best choice. If you really do have 6-way genotypes (i.e. differentiating single-nulls from homozygotes, at least approximately) then you _could_ use 're_est_ALF' but the problems are:

  're_est_ALF' does not use the extra statistical information on single-null vs homozygote that is available with 6-way genotypes, whereas 'est_ALF_6way' does;

  and in the absence of that information, it's much slower than 'est_ALF_ABO_quick' would be, so why not use that instead?!


USAGE

re_est_ALF( snpg)
est_ALF_ABCO( lociar, geno_amb = attr( lociar, 'geno_amb'))


ARGUMENTS

  snpg: a 'snpgeno' object.

  lociar: a 'loc.ar' or 'snpgeno' object, normally with a 'geno_amb' attribute from 4-way genotyping (see next).

  geno_amb: a set of 4-way genotypes. In CSIRO's 6-way genotyping pipeline (in package 'genocalldart', not for general use), 4-way genotypes get called first and stored in a 'geno_amb' attribute, before making initial allele frequency estimates with 'est_ALF_ABCO'. those are then used as starting values in 6-way genotype-calling and allele-frequency estimation, and the 'geno_amb' attribute gets discarded..


VALUE

're_est_ALF' returns the input, adding a 4-column matrix 'pbonzer' to the '$locinfo' attribute, plus attributes 'gobs' and 'gpred' showing observed and expected counts of each genotype per locus. 'pbonzer' is what you want for subsequent calculations. 'est_ALF_ABCO' is similar but creates a new attribute '$locinfo$pambig', instead of creating/modifying 'pbonzer'.


SEE.ALSO

'est_ALF_ABO_quick' for normal users, 'est_ALF_6way' for special 6-way people.


EXAMPLES

head( bluefin$pbonzer) # C alleles sometimes present
bluefin$locinfo$pbonzer <- NULL ## remove pre-existing ALFs
bluefin$locinfo$snerr <- NULL ## remove pre-existing snerr
bluefin <- re_est_ALF( bluefin)
head( bluefin$locinfo$pbonzer) # C-alleles have gone to 0
}")

)

"est_ALF_ABO_quick" <-
structure( function(
    x=NULL,
    AB,
    AAO,
    BBO,
    OO,
    tol= 1e-7,
    EMtol= 1e-3,
    quietly= FALSE,
    MAX_AITKEN= 40,
    return_unconverged= FALSE
){
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

stopifnot( my.all.equal( x@diplos, genotypes4_ambig) |
    my.all.equal(x@diplos, genotypes6))

    output6way <- FALSE
    if(my.all.equal( x@diplos, genotypes6)) {
      oldx <- x
      x@diplos <- genotypes4_ambig
      x[ oldx == "AA"] <- "AAO"
      x[ oldx == "AO"] <- "AAO"
      x[ oldx == "BB"] <- "BBO"
      x[ oldx == "BO"] <- "BBO"
      x[ oldx == "AB"] <- "AB"
      x[ oldx == "OO"] <- "OO"
      output6way <- TRUE
    }

    AB <- colSums( x=='AB')
    AAO <- colSums( x=='AAO')
    BBO <- colSums( x=='BBO')
    OO <- colSums( x=='OO')
  } # if !null x

  n <- AB + AAO + BBO + OO
  LLL <- length( AB)

stopifnot( all( AB >= 0), all( AAO >= 0), all( BBO >= 0), all( OO >= 0),
    n > 0, length( AAO)==LLL, length( BBO)==LLL, length( OO)==LLL)

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

  # Corner cases: no A's and/or no B's
  patho <- which( (AB==0) & ((AAO==0) | (BBO==0)))
  if( length( patho)){
    alpha[ patho] <- beta[ patho] <- omega[ patho] <- 0
    has_OO <- OO[ patho] > 0
    omega[ patho[ has_OO]] <- (AAO + BBO + 2*OO)[ patho[ has_OO]]
    has_A <- AAO[ patho] > 0
    alpha[ patho[ has_A]] <- 1-omega[ patho[ has_A]]
    has_B <- BBO[ patho] > 0
    beta[ patho[ has_B]] <- 1-omega[ patho[ has_B]]
  }

  iAit <- seq_along( AB) %except% patho
  nits <- 0
  n_super_its <- 0 # for Aitkening
  # MAX_AITKEN <- 20 # few will take more than 10
  while( length( iAit) && (n_super_its < MAX_AITKEN)){
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

    # Next will not check OK during update (but that's OK :) cos then it relies entirely on EM...
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
warning( sprintf(
    "Still %i not-fully-converged loci after maximum [%i] Aitken superloops",
      length( iAit), MAX_AITKEN))
    if( return_unconverged){
return( iAit)
    }
  }

  if( missing( x)){
    mat <- cbind( alpha, beta, omega)
    mat@nits <- nits
return( mat)
  } else {
      if(output6way) { x <- oldx }
      x@locinfo$pbonzer <- matrix( c( alpha, beta, 0*beta, omega), ncol=4,
        dimnames=list( NULL, cq( A, B, C, O)))
return( x)
  }
}
, doc =  mvbutils::docattr( r"{
est_ALF_ABO_quick      package:kinference
est_ALF_nonulls

Estimate allele frequencies with nulls


DESCRIPTION

'est_ALF_ABO_quick' is the recommended way to estimate allele frequencies from called genotypes, _provided that_ you are using biallelic SNPs wtih 4-way genotyping and you believe that null alleles are a real, repeatable thing in your data. If you have biallelic SNPs but you don't believe you have repeatable, heritable null alleles, then use 'est_ALF_nonulls' instead. If you are using genuine 6-way genotyping with nulls, see *Beyond 4 way genotypes*.

So... 'est_ALF_ABO_quick' is for fast maximum-likelihood estimation of A (nominally major), B (nominally minor), and O (bona fide null) allele frequencies for a set of loci. It uses the EM algorithm plus Aitken acceleration; this means the whole calculation can be vectorized across loci, which (compared to direct maximization of the log-likelihood) more than compensates for the notorious inefficiency of EM (and also Aitken helps _a lot_). See the MS or the code for more details. Missing data not allowed.

'est_ALF_nonulls' is very simple; it just looks at the ratio of A to B alleles across all samples at each locus, and sets the null frequency to 0. This _can_ cope with missing data, and any sample-loci recorded as "OO" (double null) is treated as such. Thus, you could do 'est_ALF_nonulls' followed by a random-imputation step to fill in the missings, if you really can't unmissingize them (eg by using a more confident genotype-calling algorithm to your raw data). See *Examples*.


USAGE

est_ALF_ABO_quick(x = NULL, AB, AAO, BBO, OO,
    tol = 0.0000001, EMtol = 0.001,  quietly = FALSE,
    MAX_AITKEN= 40, return_unconverged= FALSE)
est_ALF_nonulls( x=NULL, AB, AA, BB,
  pbonzer_format= FALSE)


ARGUMENTS

  x: a 'snpgeno' object, or NULL to use the next 4 or 3 args explicitly. 'diplos(x)' should either be 'genotypes6' or 'genotypes4_ambig'. For 'est_ALF_nonulls', the latter is treated as if AAO always truly means AA, and BBO means BB; OO is treated as missing.

  AB, AAO, BBO, OO, AA, BB: vectors (per locus) of counts of these genotypes. Can't mix with non-null 'x'.

  tol: final convergence tolerance (in Aitken steps)

  EMtol: tolerance within the EM steps; after this is achieved, try an Aitken step

  quietly: if TRUE, then at the end print information on the number of iterations required

  MAX_AITKEN: Maximum number of Aitken-accelerations to allow. Most loci converge within 10. We have seen a few where over 30 Aitkens are needed; however, the results after just 10 were still pretty good. You can try pushing it higher, but do check your weird loci for extreme weirdness (see next argument).

  return_unconverged: if TRUE and some loci still haven't converged after 'MAX_AITKEN' iterations, then just return the indices of those loci. By default (ie if this is FALSE), you'll just get the estimates with a warning, and that's probably fine.

  pbonzer_format: if 'est_ALF_nonulls' is called directly on 'AB' etc rather than on 'x', then 'pbonzer_format' determines the format of the returned allele frequency estimates:TRUE means you get a 4-column matrix suitable for 'kinference' (to go in '$locinfo$pbonzer') and FALSE means you just get a vector of the major (A) allele frequencies.


VALUE

If 'x' is supplied, then its 'locinfo' attribute will be augmented with the 'pbonzer' (allele frequency) matrix required by most 'kinference' functions. Note that 'pbonzer' has 4 columns always, so here the 3rd column ("C") is set to zero. If 'x' is not supplied, then a 3-column matrix matrix is returned. Rowsums of the matrix are always 1 in either case.


BEYOND.4.WAY.GENOTYPES

'est_ALF_ABO_quick' does actually accept 6-way genotypes, but it would be a bit weird to use it, because the first thing it then does is to re-encode the genotypes as 4-way, thus sacrificing some statistical information. It starts by merging the single-nulls with true homozygotes, and then tries to "unmerge" them statistically! It would be better to use 'est_ALF_6way', though you might need to first run 'est_ALF_ABCO' or 'est_ALF_ABO_quick' to get starting values.

'est_ALF_ABCO' (qv) is a much slower version that can handle triallelic SNPs (ie potentially with a C allele). Like 'est_ALF_ABO', it does not distinguish between single-nulls and homozygotes. It should give the same results as 'est_ALF_ABO' for loci without a C allele. Syntax is a bit different because of its role in the legacy 6-way pipeline, so you might prefer to use 're_est_ALF' (qv) which does the same thing but with different (easier?) formatting of the input.

'est_ALF_6way' (qv) uses 6-way genotypes, where single nulls are called separately from homozygotes, but potentially with error. It requires error estimates in '$locinfo$snerr'. Unless you are dealing with a legacy 6-way dataset, you don't want to go near this!


EXAMPLES

#dropbears$locinfo$pbonzer <- NULL ## no population allele frequency estimates!
#dropbears <- est_ALF_ABO_quick( dropbears)
#head( dropbears$locinfo$pbonzer)
#dropbears <- est_ALF_nonulls( dropbears)
#head( dropbears$locinfo$pbonzer)
## Randomly make some values missing...
#nb <- nrow( dropbears)
#nl <- ncol( dropbears)
#nmiss <- round( 0.1 * nb * nl)
#missij <- cbind( rsample( nmiss, 1:nb, replace=TRUE),
#    rsample( nmiss, 1:nl, replace=TRUE))
#missbears <- dropbears
#missbears[ missij] <- 'OO'
#missbears <- est_ALF_nonulls( missbears)
#plot( missbears$locinfo$pbonzer[,1], dropbears$locinfo$pbonzer[,1]) # very similar
#abline( 0, 1)
#imputor <- function( x){
#  # I don't *recommend* this, but you *could* randomly impute "missing" (OO) genotypes like this
#  # if you are sure there's no nulls
#  # I haven't tested the code!
#  # Caveat emptor... and read it carefully to work out what it's
#  # (hopefully)
#  # doing...
#  pA <- x$locinfo$pbonzer[,1]
#  pB <- 1-pA
#  pAB <- 2*pA*pB
#  pAA <- sqr( pA)
#  pBB <- sqr( pB)
#  misso <- which( x=='OO', arr.ind=TRUE)
#  r <- runif( nrow( misso))
#  x[ misso] <- 'BB' # default
#  x[ misso[ r < (pAB+pAA)[ misso[,2]],] <- 'AA'
#  x[ misso[ r < pAB[ misso[,2]],] <- 'AB'
#return( x)
#}
}")

)

"est_ALF_nonulls" <-
function(
  x=NULL,
  AB,
  AA,
  BB,
  pbonzer_format= FALSE
){
  if( !is.null( x)){
stopifnot( missing( AB), missing( AA), missing( BB))

    # Either 6way (preferable) or 4way is OK
    gbasics::define_genotypes()
stopifnot( my.all.equal( x@diplos, genotypes4_ambig) |
    my.all.equal(x@diplos, genotypes6))

    AB <- colSums( x=='AB')
    if( my.all.equal( diplos( x), genotypes6)){
stopifnot( sum( x=='AO')==0, sum( x=='BO')==0)
      AA <- colSums( x=='AA')
      BB <- colSums( x=='BB')
    } else if( my.all.equal( diplos( x), genotypes4_ambig)){
      AA <- colSums( x=='AAO') # let's hope you mean it!
      BB <- colSums( x=='BBO')
    } else {
stop( "diplos(x) must be 'genotypes6' or 'genotypes4_ambig'")
    }
  } # if !null x

  n <- AB + AA + BB
  LLL <- length( AB)

stopifnot( all( AB >= 0), all( AA >= 0), all( BB >= 0),
    n > 0, length( AA)==LLL, length( BB)==LLL)

    pA <- (2*AA + AB) / (2*n)

    pbonzer <- cbind( A=pA, B=1-pA, C=0, O=0)
    if( is.null( x)){
return( if( pbonzer_format) pbonzer else pA)
    } else {
      x$locinfo$pbonzer <- pbonzer
return( x)
   }
}


"find.root" <-
function(
  f,
  start, step,
  fdirection = ifelse(f(start + step, ...) > f0,
      "increasing", "decreasing"),
  target = 0,
  min.x =  -Inf, max.x = Inf,
  args.to.uniroot = list(),
  ...
){
## Copied directly into 'kinference' from package 'handy2' to avoid faffing
## No longer so necessary thx2 changes in stats::uniroot, but don't wanna risk breaking anything in 'kinference' with last-minute changes
## Prolly should move to 'mvbutils' or something

  # 2024: added some sanity checks and to avoid array warnings
  target <- as.vector( target)
  f0 <- as.vector( f(start, ...))
stopifnot(
    length( target) == 1,
    length( f0) == 1
  )

  if( f0==target)
return( start)

  step <- abs(step) * ifelse(xor(f0 < target, fdirection == "decreasing"), 1, -1)
  bound <- ifelse(xor.thing <- (step < 0), min.x, max.x)
  for (bracket_iteration in seq_len(1024L)) {
    new <- start + step
    if(xor(new > bound, xor.thing))
      new <- (start + bound)/2
    if (!is.finite(new) || new == start)
      stop("Cannot bracket target within the supplied bounds")
    f1 <- as.vector( f(new, ...))
    if (length(f1) != 1L || !is.finite(f1))
      stop("Non-finite or non-scalar function value while bracketing")
    if (f1 == target) return(new)
    if(xor(f0 < target, f1 < target))
      break
    if (bracket_iteration == 1024L)
      stop("Cannot bracket target after 1024 iterations")
    start <- new
    f0 <- f1
    step <- step * 2
  }
  o <- order(x <- c(start, new))
  fvals <- c(f0, f1) - target
  dotnames <- paste(c(names(f)[1], names(list(...))), collapse = ",")
  ff <- function(...) as.vector( f(...)) - target

#  ff <- c(f[ - length(f)], list(FFFF = 1, target = 1), parse(text = "FFFF(" %&% dotnames %&% ") - target"))
#  mode(ff) <- "function"

  do.call("uniroot", c(list(ff, x[o]),
      # f.lower = fvals[o[1]], f.upper = fvals[o[2]]),
      args.to.uniroot, list(...)))$root
}


"find_duplicates" <-
structure( function(
  snpg,
  subset1= 1 %upto% nrow( snpg),
  subset2= subset1,
  max_diff_loci,
  limit_pairs= 0.5*nrow(snpg),
  nbins= 50,
  maxbin= ncol( snpg)/2,
  show_plot = TRUE,
  ij_numeric= is.null( rowid_field( snpg))
){ #################
  define_genotypes() # must be first
  basic_sanity_checks_pairfinding()
stopifnot( # ... and....
    max_diff_loci <= maxbin, # will bork on missing
    nbins >= 0 # 0 is OK; just return actual pairs below max_diff_loci
  )

  # Count #loci with different 4way genos. Errors in 4ways should be low.
  temp_snpg <- make_6way_into_4way( snpg)

  # Sort loci to get most informative/random ones first
  # use snpg1 for this, arbitrarily
  gtab <- matrix( 0, 4, ncol( snpg),
      dimnames=list( genotypes4_ambig, NULL))
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

  # Gotta ensure number of bins matches nbins even with rounding error...
  # nbins=0 is OK cos it saves heaps of time,
  # but obvs you don't get binned results
  binterval <- maxbin / max( nbins-1, 1)
  # bin *starting* at maxbin counts all bigger; avoid woe when nbins=0
  bins <- seq( from= 0, to= maxbin, length= max( 1, nbins-1))[-1]
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
  result <- with( result, data.frame(
      ndiff=big_similar, i= subset1[ big_i], j= subset2[ big_j]))
  if( !ij_numeric){
    result <- make_ij_character( result, temp_snpg)
  }

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
, doc =  mvbutils::docattr( r"{
find_duplicates      package:kinference
find_HSPs
find_POPs
find_POPs_lglk
%upto%

Kin-finders for loads-of-SNPs datasets


DESCRIPTION

These take a 'snpgeno' dataset that has been processed as far as 'check6and4' (and for HSPs, 'kin_power') and find various relations between the samples. Relationships include duplicates (DUPs/dupes/dups), parent-offspring pairs (POPs) and half-sibling pairs (HSPs) or other 2nd-order kin, plus of course unrelated pairs (UPs). You can specify the same or different subsets of the 'snpgeno' for comparison: e.g., first subset for the adults, second for the juveniles.

There are two versions aimed at POPs currently called 'find_POPs' and 'find_POPs_lglk'. The former uses a "weighted pseudo-exclusion" ("wpsex") statistic that allows for null alleles and is robust to genotyping errors. The latter uses a likelihood-based statistic (again allowing for nulls), but you do have to provide a guesstimate of genotyping error rate (to robustify the calculation- otherwise, a single genotyping error in a true POP could give a log-likelihood of -Inf). 'find_POPs_lglk' is newer, easier to explain, and perhaps less arbitrary, but we have used the "wpsex" version on all our real CKMR datasets (>10). Time will tell whether one is better/easier than the other; finding POPs ought to be pretty easy, so the results really should be the same.

'find_HSPs' should really be called 'find_2KPs' because it cannot discriminate amongst second-order kinships; there is no way to distinguish genetically between HSPs, Grandparent-Grandchild Pairs, and Full-Thiatic Pairs (eg aunt/nephew) with 'snpgeno' data alone. But, for historical reasons, it's still called 'find_HSPs'. Note that 'find_HSPs' can also be tricked into targeting [some] other types of kin, such as 3KPs; see *Details*, but watch out.


USAGE

find_duplicates(
  snpg,
  subset1 = 1 %upto% nrow(snpg),
  subset2 = subset1,
  max_diff_loci,
  limit_pairs = 0.5 * nrow(snpg),
  nbins = 50,
  maxbin = ncol(snpg)/2,
  show_plot = TRUE,
  ij_numeric= is.null( rowid_field( snpg))
)
find_HSPs(
  snpg,
  subset1 = 1 %upto% nrow(snpg),
  subset2 = subset1,
  limit_pairs = 0.5 * nrow(snpg),
  keep_thresh,
  eta = NULL,
  nbins = 50,
  minbin = NULL,
  maxbin = NULL,
  ij_numeric= is.null( rowid_field( snpg))
)
find_POPs(
  snpg,
  subset1 = 1 %upto% nrow(snpg),
  subset2 = subset1,
  limit_pairs = 0.5 * nrow(snpg),
  keep_thresh,
  eta = NULL,
  nbins = 50,
  maxbin = NULL,
  WPSEX_UP_POP_balance = 0.99,
  ij_numeric= is.null( rowid_field( snpg))
)
find_POPs_lglk(
  snpg,
  subset1 = 1 %upto% nrow(snpg),
  subset2 = subset1,
  gerr,
  limit_pairs = 0.5 * nrow(snpg),
  keep_thresh,
  eta = NULL,
  nbins = 50,
  minbin = NULL,
  maxbin = NULL,
  ij_numeric= is.null( rowid_field( snpg))
)


ARGUMENTS

  snpg: a 'snpgeno' object

  subset1, subset2: numeric vectors of which samples to use (not logical, not negative). Defaults to all of them. Iff 'subset1' and 'subset2' are identical, only half the comparisons are done (i.e., not _i_ with _j_ then _j_ with _i_). Some sanity checks are made.

  max_diff_loci: ('find_duplicates') max number of discrepant 4-way genotypes to tolerate in "identical" fish. Only the pairs with fewer than 'max_diff_loci' discrepancies will be retained. Try increasing this from say 10 upwards, and hopefully nothing much will change (though at some point things will change a lot, as you get into the non-duplicate bit of the distribution). Once you have identified duplicates, see 'drop_dups_pairwise_equiv' for how to remove them.

  limit_pairs: Integer. Defines the _maximum_ number of candidate pairs to keep. Will provide a warning if the number of identified pairs equals 'limit_pairs'.

  nbins, minbin, maxbin: 'find_XXX' functions summarise their pairwise comparison statistics into bins (in the part of the range where exact values are uninteresting), as well as returning specific pairs that pass the "interesting" threshold. 'nbins' sets the number of bins, 'minbin' sets the top value of the lowest bin (so that bin stretches from -Inf to 'minbin' for HSPs); 'maxbin' sets the highest. For HSPs, the minimum is 3 bins (-Inf:minbin),[minbin:maxbin),[maxbin:Inf). 'minbin' is not used for duplicates, nor for the non-likelihood 'find_POPs', since the statistics there are defined so that the lowest possible value is 0. The defaults for 'minbin' and/or 'maxbin' may not be what you need in all cases, so be prepared to select manually and then re-run. For duplicates, where calculations can be slow for big datasets, you can set 'nbins=0' to disable binning and focus instead on just finding the pairs with fewer than 'max_diff_loci' discrepancies. Each pairwise calculation normally loops over all the loci, but is aborted when the running total of discrepant loci reaches 'maxbin' (or, if 'nbins=0', when it reaches 'max_diff_loci'), thus saving considerable time. It is therefore not sensible to have 'maxbin<max_diff_loci' (think about it!).

  show_plot: whether to plot log histogram. Regardless, plot will not be shown if other arguments would lead to stupid result (e.g. no bins...).

  ij_numeric: if FALSE, use the rowid field (see 'with_rowid_field') to label the pair-members, rather than their row numbers.

  keep_thresh: ('find_HSPs' and 'find_POPs') is the analog of 'max_diff_loci' for 'find_duplicates'. It determines which pairs to retain for individual inspection. For 'find_HSPs' and 'find_POPs_lglk', this is the lowest retained PLOD; for wpsex-based 'find_POPs', it's the highest retained 'wpsex'. Set it with the aim of including anything interesting (ie not _missing_ any interesting pairs) and do expect false positives; that is, be willing to have some weaker kin in there, and to subsequently filter those out yourself, as per vignette. For HSPs, and for POPs with 'find_POPs_lglk', values like 0 (near the HTP mean) or -5 are a good start. For POPs with wpsex-based 'find_POPs', experiment (perhaps starting with 0.1). You may have to re-run the function a few times if you have been too brutal or too generous here - though "too generous" can be fixed post hoc just by filtering the result, as long as you haven't generated tooooo many pairs (see parameter 'limit_pairs').

  eta: ('find_HSPs', 'find_POPs_lglk', and 'find_POPs') Not essential; limit for calculating empirical mean and var PLOD, to compare with theoretical 'mean_UP' and 'var_UP'. If you care about this (and you might not, since for 'find_HSPs' the observed/expected binwise comparison is perhaps clearest), then set it to somewhere above 0 that should include almost all UPs and exclude most strong kin; that's an _upper_ limit for HSPs and lglk-based POPs, but a _lower_ limit for wpsex-based POPs in 'find_POPs'. To use 'eta', be prepared to look at the histograms and think. The general idea is that the number of UPs should dominate any other kin-type in large sparsely-sampled datasets, so there shouldn't be much problem even if you do accidentally "contaminate" the empirical UP statistics with a few weakish kin at the top end.

  WPSEX_UP_POP_balance: ('find_POPs') loci receive a weight which is proportional to (difference in probability of pseudo-exclusion between UP and POP) / (variance of indicator of pseudo-exclusion). But, should this be variance assuming UP or POP? This parameter sets the balance; bigger values make it more UPpity, so placing more emphasis on avoiding false-positives, which is probably the Right Thing To Do in most cases. 0.99 could be completely fine... but hopefully the value won't affect the result much anyway.

  gerr: ('find_POPs_lglk') Genotyping error rate (apart from any AA/AO-type errors)- which had better be a small number. You have to pick it yourself, but it is only used to "robustify" the lglk-based (PO)PLOD for testing POPs vs UPs, and thus can be a rough guesstimate. FWIW we have used 0.01 (i.e. 1%), which is considerably higher than suggested by an analysis of our replicate samples, but it is "safe" while still being small enough not to muck up overall statistical performance. You should really do the same thing yourself, and if you are very paranoid then try sensitivity analyses; but in practice, the results of 'find_POPs_lglk' are liable to be so clear-cut that you may not feel it necessary to try more than one small value.


DETAILS

Some categories will "catch" others (e.g. 'find_HSPs' will certainly include any POPs too), so you may need the splitter routines such as 'split_POPs_from_HSPs' afterwards. The safest general-purpose strategy - but often _not_ the most sensible, if your data is nicely organized and you know what you want - would be:

'find_duplicates' and then get rid of them

'find_HSPs' to get _all_ kin (though you will usually have to sacrifice some HSPs to false-neg because you'll need a threshold)

'split_POPs_from_HSPs' to split HSPs from POPs/FSPs

'split_POPs_from_FSPs' to split the latter.

The non-splitter functions, i.e. 'find_XXX', might be run on huge numbers of samples, entailing a 'choose(huge, 2)' number of comparisons. You don't want all those individual comparison results, and your computer certainly wouldn't enjoy trying to keep them! So the general idea is to set a threshold for what constitutes "maybe worth keeping individually" (that you expect will be generous enough to contain everything you _do_ want, plus some dross), and then to retain just binned counts of the relevant comp statistic for all comps (usually, the vast majority) which don't make your threshold.

In addition, the 'limit_pairs' argument is there to prevent your computer locking out with bazillions of unwanted pairs (in case you guess the bin limit inapproriately); the comparisons will be stopped if 'limit_pairs' is hit, with a warning. In that case, you probably need to change a threshold, or re-run with larger 'limit_pairs'. The default isn't meant to correspond to any biomathematical logic, it's just to stop blue smoke coming out your USB ports.

For 'find_duplicates', there are at least two different use-cases. First, you might want an initial run on a non-too-large subset of your data, to check that dups _can_ be clearly distinguished and to look at typical extent of genotyping errors (based on clear duplicates that don't match at every locus). For that, you can set 'nbins' and choose some reasonable guess as to 'max_diff_loci' (say, 50 loci). Because you set 'nbins>0', _every_ pair (almost...) gets checked at _all_ loci, so it can be slow. Thus, if you have done this before and have a good sense of "how bad can a real duplicate be?", then set 'nbins=0' (and 'max_diff_geno' to a small but safe value that won't miss any realistic duplicate-with-genotyping-error) so it will abort a comparison early as soon as it reaches 'max_diff_geno' differing loci. That saves a _lot_ of time on big datasets! You won't get a histo of number-of-diffs, but you don't need one for that use-case. The "almost" is that 'find_duplicates' uses "transitivity" (if A is a dup of B and of C, then we don't need to check B vs C), so it only counts differences for not-yet-known duplicates _based on_ 'max_diff_loci'. To discard duplicates and to find entire equivalence-classes of duplicates, e.g. from a control specimen included in numerous plates, see 'drop_dup_pairwise_equiv'.

'find_HSPs' relies on pre-computed values of "LOD" and "PUP" that have been set by 'kin_power'. Normally you would call the latter with 'k=0.5', since that's what HSPs are. However, the devious user can try _different_ values of 'k'- which is how 'find_POPs_lglk' works- and then the target of 'find_HSPs' will become "kin with that value of 'k'". Be very careful!!!


VALUE

A 'data.frame' with extra attributes (see below) and at least 3 columns: statistic 'PLOD' or 'wpsex' or 'ndiff' (number of mismatching genotypes), then 'i' and 'j' which are the indices/rows in 'snpg' of the two members of each pair. (In ancient times 'i' and 'j' misleadingly referred to the subsets instead, but we have moved on!) The attributes in all cases include 'bins' (upper boundaries), some kind of count statistic for number of comparisons in each bin (names vary), 'binprobs' (theoretical CDF for UPs in 'find_HSPs'; should also exist for _POPs_ in the old "wpsex"-based 'find_POPs', but currently doesn't), some of the input parameters, and the 'call' that invoked the function. 'find_POPs' adds a column named 'nABOO', showing the number of AB/OO exclusions for that potential POP. This is a useful additional diagnostic; it should be close to 0 for true POPs (it can only result from genotyping error or mutation, whereas AAO/BBO can result from nulls). For UPs, we have seen values typically in the low 20s, which is pretty good separation. A 2D scatter-plot of 'wpsex' and 'nABOO' can be more informative than either statistic on its own. 'find_HSPs' and 'find_POPs' have a bunch of extra attributes which should be reasonably clear. For 'find_HSPs', 'mean_sub_PLOD' and 'var_sub_PLOD' are the empirical means & var below 'eta', andy they should be close to 'mean_UP' and 'var_UP' _iff_ 'eta' has been chosen sensibly. For 'find_POPs', the same goes for 'mean_wpsex_hi' and 'var_wpsex_hi'. For duplicates, not _all_ pairwise duplicates are recorded, unless the subsets are different - otherwise you could have quadratic horror of enormous numbers of pairs arising from a cluster of say 100 identical controls! Since duplication is transitive (ie if _i_ & _j_ are the same, and _i_ & _k_ are the same, then _j_ & _k_ must also be the same), only the necessary ones are recorded to allow you to filter out yourself afterwards. For example, if samples 1, 3, 5, and 6 are all duplicates, you'll get this:

  i j:

  [1] 3 1:

  [2] 4 3:

  [3] 6 4: but you won't see the pairings for 1/4, 1/6, 3/6. If you just want to strip out all duplicates bar one in each group (and you don't care which one is kept), then you can use the function 'drop_dups_pairwise_equiv' - see *Examples*. For POPs and HSPs, the items below are also returned as attributes (which can be more conveniently accessed by '@' if 'atease' is loaded, as per *Examples*). The main point is that the "boring" below-threshold pairs get put into bins and are not kept individually. You often don't need to access these directly, since 'histoPLOD' may be enough. The names sometimes change depending on which statistic is being used.

  - eta, keep_thresh: see above.

  - <mean/var>_sub_PLOD, <mean/var>_hi_wpsex: empirical values for the statistic when it is on the boring (UP) side of 'eta' (ie nearly always).

  - <mean/var>_<kin>: (where 'kin' is 'UP', 'POP', or sometimes 'HSP', 'FSP') predictions from theory, to compare with previous and show graphically.

  - n_<stat>_in_bin: (where 'stat' is 'PLOD', 'wpsex', or 'ndiff'):number of pairs whose statistic fell into each bin.

  - bins: Upper cutpoints for the bins.


EXAMPLES

## find_duplicates
library( atease) # @ used below for clarity
library( mvbutils) # IDNK if this is needed explicitly
define_genotypes()    ## creating 'genotypes4_ambig' etc
x <- matrix(sample(c("AAO", "AB", "BBO", "OO"), 10000, TRUE,
    prob = c(0.2, 0.45, 0.2, 0.15)), nrow = 100, ncol = 100)
minisnpg <- snpgeno(x = x, diplos = genotypes4_ambig)
## seed some duplicates in. Sample 2 will be a copy of sample 1,
## and samples 21, 22, and 23 will be an exactly-matching group.
minisnpg[2,] <- minisnpg[1,]
minisnpg[ 22,] <- minisnpg[ 23,] <- minisnpg[ 21,]
## seed some inexact duplicates in. Make sample 3 nearly match
## sample 4, and the same for samples 13 and 14
minisnpg[4,1:80] <- minisnpg[3,1:80]
minisnpg[14,1:80] <- minisnpg[13,1:80]
## find exact duplicates
exact <- find_duplicates( minisnpg, max_diff_loci = 0, show_plot=FALSE)
exact ## finds ij pairs {1, 2} and the group {21, 22, 23}
## as pairs {21, 22} and {22, 23}
## find (in-)exact duplicates
## first a plot, as a guide to where to look...
find_duplicates( minisnpg,
    maxbin= 100, max_diff_loci= 1, show_plot=TRUE)
## looks like there's a few inexact ones
inexact <- find_duplicates( minisnpg,
    maxbin= 100, max_diff_loci= 25, show_plot=FALSE)
inexact ## finds the exactly-matching pairs as before, plus
## the inexactly-matching row pairs {3, 4} and {13, 14} with
## >0 differences
## to remove duplicates, keeping only one member of each
## group, use drop_dups_pairwise_equiv
droppers <- drop_dups_pairwise_equiv( inexact[,2:3])
droppers     ## note that _all but one_ of each _group_ of
## (near-)duplicates is included in 'droppers'
## Drop all-but-one of each set of duplicates:
minisnpg_nodups <- minisnpg[ - droppers,]
## find_HSPs (PLOD_HU)
## bluefin data
## stripped-down data-cleaning for example - see the
## vignette for approach for real data!
pvals <- check6and4( bluefin, thresh_pchisq_6and4 = c( 0.001, 0.0001))
bluefin_1 <- bluefin[ , pvals$pval4 > 0.01]
ilglks <- ilglk_geno( bluefin_1)
bluefin_2 <- bluefin_1[ ilglks > -1050,]
bluefin_3 <- est_ALF_ABO_quick( bluefin_2)
bluefin_4 <- bluefin_3[ , bluefin_3$locinfo$pbonzer[,"B"] > 0.02]
bluefin_5 <- kin_power( bluefin_4, k = 0.5)
dups <- find_duplicates( bluefin_5, max_diff_loci = 20)
bluefin_6 <- bluefin_5[ -c(drop_dups_pairwise_equiv( dups)) ,]
hsps <- find_HSPs( bluefin_6, keep_thresh = 0)
histoPLOD( hsps, log=FALSE, lb = 0, fullsib_cut = 75)
head( hsps)     ## PLODs and row numbers for each pair member in bluefin
library( atease)  ## and then it's sooo much easier to just write...
hsps@mean_HSP     ## mean expected PLOD for true HSPs. mean_UP (unrelated),
## mean_POP (parent-offspring), mean_FSP (full-sibling) follow the
## same format
## subset comparisons
## limit comparisons to those between animals on plate 3 and animals on
## the other two plates. Useful when, e.g., looking for kinships between
## adults and juveniles, but not kinships between adults and other
## adults. For demonstration, treat plate 3 as adults and plates 1 and 2
## as juveniles (they're not really).
adults <- which( bluefin_6$info$Our_plate == "plate3")
juvs <- which( bluefin_6$info$Our_plate %in% c( "plate1", "plate2"))
hsps_subset <- find_HSPs( bluefin_6,
    subset1 = adults, subset2 = juvs, keep_thresh = 0)
hsps_subset$iplate <- bluefin_6$info$Our_plate[ hsps_subset$i]
hsps_subset$jplate <- bluefin_6$info$Our_plate[ hsps_subset$j]
## all pairs are between a plate 3 fish and a fish on plate 1 or plate 2
head(hsps_subset)
## find_POPs and find_POPs_lglk
pops <- find_POPs(bluefin_6, keep_thresh = 0.1)
## keep_thresh = 0.1 is not magic, but often OK. If you
## don't see a big spike out to the right of your plot,
## you should generally set keep_thresh higher
pops@mean_UP     ## expected mean 'wpsex' for unrelated pairs
pops@var_UP     ## expected variance 'wpsex' for unrelated pairs
hist( pops$wpsex)     ## true POPs should have a wpsex near zero, and
## a nABOO of exactly zero, unless sequencing error has occurred
with( pops, plot.default( wpsex ~ nABOO))
## plot the full distribution from binned records
plot( pops@bins, pops@n_wpsex_in_bin, type = "s")
## add expected mean wpsex for unrelated pairs
abline( v = pops@mean_UP, col = kinPalette("UP"), lwd = 2)
pops_lglk <- find_POPs_lglk(bluefin_6, keep_thresh = 0, gerr = 0.01)
pops_lglk@mean_UP     ## as before
pops_lglk@mean_POP
## plot the full distribution from binned records
histoPLOD( pops_lglk, log=TRUE)
}")

)

"find_dups_with_missing" <-
structure( function(
  snpg,
  subset1= 1 %upto% nrow( snpg),
  subset2= subset1,
  max_diff_ppn,
  limit= 10000,
  ij_numeric= is.null( rowid_field( snpg))
){ ##############
  define_genotypes()
  basic_sanity_checks_pairfinding()
stopifnot( # and...
    !missing( max_diff_ppn)
  )

  og <- options( vecless.print=FALSE)
  on.exit( options( og))

  temp_snpg <- make_6way_into_4way( snpg)

  # Sort loci to improve speed (highest chance of mismatch first)
  # so that exit on non-dups is quicker
  gtab <- matrix( 0, 3, ncol( snpg),
      dimnames=list( genotypes4_ambig %except% OO, NULL))
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
  attributes( temp_snpg) <- attributes( temp_snpg)[ 'dim']
  # ... like unclass but more drastic
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
  result <- with( result, data.frame(
      ppn_diff= big_ndiff / big_ncomp,
      i= subset1[ big_i], j= subset2[ big_j],
      ndiff=big_ndiff, ncomp=big_ncomp))
  if( !ij_numeric){
    result <- make_ij_character( result, temp_snpg)
  }


  result@call <- sys.call()

return( result)
}
, doc =  mvbutils::docattr( r"{
find_dups_with_missing      package:kinference
%upto%

Duplicate-finding with some missing genotypes


DESCRIPTION

This function is not for CKMR datasets (where missingness is not allowed) but rather for "Gene-Tagging" (individual mark-recapture using genotypes as tags). For GT, the genotyping method has to be cheap rather than high-quality, so that some genotypes just end up missing (ie the genotyping pipeline has decided they are unscorable). Rather than "imputing" those missing genotypes (which you'd have to do for finding kin-pairs), for the specific case of finding duplicate samples it is probably better to just compare across the loci that _are_ both called in each pair of samples.

'find_dups_with_missing' looks at all pairwise comparisons and "accepts" any where there are not too many inconsistencies (which can arise within true duplicates, from allelic dropout when read-depth is low). This means that a very-poor-quality-DNA sample might actually match with lots of others, simply because it has so few genotypes called! So, subset your data beforehand to remove Bad Eggs.


USAGE

find_dups_with_missing(
  snpg,
  subset1 = 1 %upto% nrow(snpg),
  subset2 = subset1,
  max_diff_ppn,
  limit = 10000,
  ij_numeric= is.null( rowid_field( snpg))
)


ARGUMENTS

  snpg: a 'snpgeno' object

  subset1, subset2: numeric vectors of which samples to use (not logical, not negative). Defaults to all of them. Iff 'subset1' and 'subset2' are identical, only half the comparisons are done (i.e., not _i_ with _j_ then _j_ with _i_). Some sanity checks are made.

  max_diff_ppn: What _proportion_ of non-missing (ie scored-in-both) loci to treat as the threshold for duplicity?

  limit: if you hit this many "duplicates", it will stop, to avoid blowing out memory. It means you set 'max_diff_ppn' too high. For consistency, we should probably have called this 'keep_n' as per other 'find_XXX' functions.

  ij_numeric: if FALSE, use the rowid field (see 'with_rowid_field') to label the pair-members, rather than their row numbers.


VALUE

A dataframe with columns 'ppn_diff', then 'i' and 'j' which show pairs of samples that are within 'max_diff_ppn'. (Note that 'i' and 'j' refer to rows in 'snpg' itself, not to whatever was passed in the subsets. This is what any normal person would expect, but long ago the software worked differently...)
}")

)

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
  maxbin= NULL,
  ij_numeric= is.null( rowid_field( snpg))
){ ##################
    define_genotypes()

  basic_sanity_checks_pairfinding()
  # Also, snpg should have been thru 'prepare_PLOD_SPA' so it has @PPS
stopifnot(
    'Kenv' %in% names( attributes( snpg)),
    !missing( keep_thresh)
  )

  kinPower_change <-
      snpg@kinPower_checksum != calc_kinPower_checksum( snpg$locinfo)
  PLODSPA_change <-
      snpg@PLODSPA_checksum != calc_PLODSPA_checksum( snpg$locinfo)

  if(kinPower_change | PLODSPA_change) {
      warning("snpg$locinfo appears to have been modified after kin_power() was last called. You generally want to run find_HSPs with a locinfo that is up-to-date with itself. Consider re-running kin_power() on your snpg and trying again")
  }

  ## if snpg@diplos is genotypes4_ambig, convert whole lot to pseudo 6-way
  if( my.all.equal( diplos( snpg), genotypes4_ambig)) {
    tempsnpg <- snpg
    diplos( snpg) <- genotypes6
    snpg[ tempsnpg == OO] <- OO
    snpg[ tempsnpg == AB] <- AB
    snpg[ tempsnpg == AAO] <- AA
    snpg[ tempsnpg == BBO] <- BB
    rm(tempsnpg)
    snpg <- kin_power(snpg, k = 0.5) ## HACK!
  }

  og <- options( vecless.print=FALSE)
  on.exit( options( og))

  # For 4way loci, temporarily treat XO as XX...
  # ... have already adjusted the LOD entries so that
  # ... new_LOD6( XX/..) <- LOD4( XXO/..)
  # ... use the LOD that's in Kenv, where SPA is calculated

  ### Do we need LOD6/4/3 ?
  # extract.named( snpg@locinfo[ cq( useN, LOD6, LOD4, LOD3)])
  useN <- snpg@locinfo$useN
  temp_snpg <- snpg
  recode4to6temp <- function( x){
      x[ x=='AO'] <- AA; x[ x=='BO'] <- BB; x
    }
  recode3to6temp <- function( x){
      x[ x=='AO'] <- AA; x[ x=='BO'] <- BB; x[ x=='OO'] <- BB; x
    }

  temp_snpg[ , useN == 4] <- recode4to6temp( snpg[, useN == 4])
  temp_snpg[ , useN == 3] <- recode3to6temp( snpg[, useN == 3])
  temp_LOD <- snpg@Kenv$LOD # already done in prepare_PLOD_SPA, and based on useN to switch 6/4/3

  # Remove extranea
  attributes( temp_snpg) <- attributes( temp_snpg)[ 'dim']
  temp_snpg <- t( temp_snpg)

  li <- snpg@locinfo # convenient

  # MVB: don't trust length of seqs with non-integer steps!
  # bins <- seq(minbin+binterval, (minbin + (nbins*binterval)), binterval)

  EHSP <- sum( li$E_HSP) / 3 # who knows?!
  if( is.null( eta)) {
    eta <- EHSP / 3 # who knows?!
  }

  if( is.null( minbin)) {
    VUP <- sum( li$V_UP)
    EUP <- sum( li$E_UP)

    ncomps <- length( subset1) * length( subset2)
    # ... might be out by 2 if s1==s2...
    # ... which is irrel on log scale

    # Bins start at 2SD below lowest expected PLOD
    # usually should be enuf to show weird stuff below min
    nSD <- abs( qnorm( 1/ncomps)) + 2
    minbin <- EUP - nSD * sqrt( VUP) # 0.003% of UPs below that
  }

  if( is.null( maxbin)) {
    EPOP <- sum( li$E_POP)
    maxbin <- EPOP + (EPOP - EHSP) * 0.5 # changed from 0.1 in 6/2025
  }
  bins <- seq( from=minbin, to=maxbin, length=nbins)
  binterval <- bins[2] - bins[1]

  # Distro of PLOD|UP via SPA. Can _slightly_ exceed 1 at far-right tail, and oscillate, so squish it with a big hammer...
  binprobs <- snpg@Kenv$CDF( bins)
  binprobs <- c( cummax( binprobs), max( c( binprobs, 1)))
  binprobs <- binprobs / tail( binprobs, 1)

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
  # ... prolly needed in find_POPs, find_duplicates, etc, too
  if( !missing(keep_thresh) && (length(xresult$big_PLOD) == limit_pairs)){
    message( "Returning just the ", limit_pairs,
        " pairs with the highest PLOD scores; increase 'limit_pairs' to get more"
      )
  }

  result <- with( xresult, data.frame(
      PLOD=big_PLOD, i= subset1[ big_i], j= subset2[ big_j]))
  if( !ij_numeric){
    result <- make_ij_character( result, temp_snpg)
  }

  attributes( result) <- c( attributes( result),
      xresult %without.name% cq( big_PLOD, big_i, big_j))

  # assign extra info as attributes
  ## result@mean_UP <- snpg@Kenv$dK( 0)   # was called mean_theory
  result@var_UP <- snpg@Kenv$ddK( 0)
  ## result@mean_HSP <- snpg@Kenv$dK( 0) + sum(snpg@locinfo$Ediff) ##
  result@mean_HSP <- sum(snpg@locinfo$E_HSP)
  result@mean_UP <- sum(snpg@locinfo$E_UP)
  result@mean_POP <- sum(snpg@locinfo$E_POP)
  result@mean_FSP <- sum(snpg@locinfo$E_FSP)
  attributes( result) <- c( attributes( result),
      returnList( bins, binprobs, eta, keep_thresh))

  result@call <- sys.call()

return( result)
}


"find_POPs" <-
function(
  snpg,
  subset1=1 %upto% nrow( snpg), subset2=subset1,
  limit_pairs=0.5*nrow(snpg),
  keep_thresh,
  eta= NULL,
  nbins= 50,
  maxbin = NULL,
  WPSEX_UP_POP_balance=0.99,
  ij_numeric= is.null( rowid_field( snpg))
){
###################
  define_genotypes()
  basic_sanity_checks_pairfinding()

  # Decide based #apparent exclusions of AA/BB form, using 4way genos, though
  # ... it's really AAO/BBO so not a true exclu but
  # ... close among pop-loci
  # Sticking with 4way genos so that genotyping errors are low

  # Instead of Nexclu, uses a wted sum of "exclus" in 4way genos
  # to max expected diff between POP and UP
  # This version doesn't allow for geno errors, but
  # does realize that AO/BO could happen
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

  # Want POP and UP means as far apart as possible on the scale of SDs...
  # but, which SD? Make it WPSEX_UP_POP_balance * SD[UP] +
  # ... (1-WPSEX_UP_POP_balance) * SD[POP]

  # Optimal wt would depend on p0 and to some extent on pA
  # wt should be 1 if p0==0 and 0 if p0==1
  delta <- pex[,'UP'] - pex[,'POP'] # mathematically I think this *can't* be -ve
  SD <- sqrt( pex * (1-pex))
  SD_combo <- WPSEX_UP_POP_balance * SD[,'UP'] +
      (1-WPSEX_UP_POP_balance) * SD[,'POP']
      # ... %*% c( WPSEX_UP_POP_balance, 1-WPSEX_UP_POP_balance)
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
    recode_6_as_pseudo4 <- function( x) {
        x[ x=='AO'] <- AA; x[ x=='BO'] <- BB; x
      }
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
        ddKK[ l, j] := sqr( ww[ l]) * retw[ l, j] / sqr( 1+retw[ l, j])
        # ... guessing that is more accurate tahn 1-1/1+x
        ddK[ j] := ddKK[ +., j]
     return( c( ddK))
     }

    K <- compile_vecless( K(0))
    dK <- compile_vecless( dK(0))
    ## ... dK(0) should be the mean of wpsex for true UPs
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
    message("Returning just the ", limit_pairs,
    " pairs with most POP-like wpsex; increase 'limit_pairs' to get more")
  }

  n_wpsex_in_bin <- result$n_wpsex_in_bin
  # SB tweak; the 'below zero' bin must be empty
  # bins <- seq(0+binterval, 0+(nbins*binterval), binterval)
  # SB addition. Bloody zero-base.

  # construct the result
  result <- with( result, data.frame(
      wpsex=big_wpsex, i= subset1[ big_i], j= subset2[ big_j]))
  if( !ij_numeric){
    result <- make_ij_character( result, temp_snpg)
  }


  # calculate nABOO, only for interesting pairs
  snpg_i <- snpg[ result$i, pop_loci]
  snpg_j <- snpg[ result$j, pop_loci]
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
  maxbin= NULL,
  # WPSEX_UP_POP_balance=0.99
  ij_numeric= is.null( rowid_field( snpg))
){ ###################
  define_genotypes()
  basic_sanity_checks_pairfinding()

  # More sanity...
stopifnot(
    !missing( gerr)
  )

  # A POP is "just like" an HSP, except ppn coin is 1.0 not 0.5
  # but we add a little bit of geno error for statistical lubrication
  # Call 'kin_power' to get LODs manually for this case
  # then 'prepare_PLOD_SPA' to handle useN stuff
  # then use HSP machinery.

  snpg <- kin_power( snpg, k= 1-gerr) # !! Captain Sneakypants !!
  # snpg <- prepare_PLOD_SPA( snpg, n_pts_SPA_renorm=201) # now automatic

  mc <- match.call()
  mc$gerr <- NULL
  mc[[1]] <- quote( find_HSPs)
  # Kinda "delayed assign" for snpg, to look here;
  # ... avoids putting whole object into result@call
  mc$snpg <- substitute( sys.frame( n)$snpg, list( n=sys.nframe()))
  result <- eval.parent( mc)

  atts <- attributes( result) %without.name% 'call'
  atts <- atts[ names( atts) %that.dont.match% 'FSP']
  attributes( result) <- atts

  # Not sure if the FSP & HSP means are OK (might be fine; HSP looks OK...)
  # Confusing cozza the sneaky use of kin_power(), too lazy to figure it out
  # ... so, don't plot them subsequently
  result@trustplot_FSP_HSP_means <- FALSE

  # ... which is slightly lower than sum( snpg@locinfo$E_POP). Cozza gerr?
  result@call <- sys.call()

return( result)
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


"hetzminoo_fancy" <-
structure( function(
  snpg,
  target=c( 'rich', 'poor'),
  hist_pars=list(),
  showPlot = TRUE
){
###################
  define_genotypes()

  # Be kind to the user with informative error messages (for once...)
stopifnot( snpg %is.a% 'snpgeno',
    all( cq( AB, OO) %in% diplos( snpg)),
    !any( grepl( 'C', diplos( snpg))))
  pbonzer <- snpg@locinfo$pbonzer %||% stop(
      "No allele freq estimates!")
  PUP4 <- snpg@locinfo$PUP4 %||% stop(
      "No PUP4; did you forget 'kin_power'?")

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
, doc =  mvbutils::docattr( r"{
hetzminoo_fancy      package:kinference

QC checks on sample heterozygosity


DESCRIPTION

This test looks for samples with anomalous numbers of heterozygotes and/or double-nulls, which can result from (i) degraded DNA or (ii) sample contamination. Useful both for finding outlier samples, and for checking whether the loci are collectively working as they should (as is assumed by all the calculations in 'kinference'). The histogram should coincide nicely with its predicted line.

The basis for the "hetzminoo" statistic is explained in the accompanying MS. Note that the 'target' argument provides two variants, "rich" and "poor", designed to test for contamination and degradation respectively. With "rich", you should look for outlying samples on the _RHS_ of the histogram (too many heterozygotes): with "poor", on the _LHS_ (too few). The choice of target affects the weightings across loci, as explained in the MS. In practice, there is often little visual difference, and a bad sample looks bad bad in both. Nevertheless, you _should_ run both variants. See 'kinference-vignette' (qv) for examples.


USAGE

hetzminoo_fancy(
  snpg,
  target = c("rich", "poor"),
  hist_pars = list(),
  showPlot = TRUE
)


ARGUMENTS

  snpg: a 'snpgeno' object, with allele frequencies already estimated and 'kin_power' (qv) already run.

  target: which potential problem to focus on.

  hist_pars: list of parameters to pass to 'hist'. If you are very sneaky, you can pass in an 'expression' to be evaluated inline instead (ie fly-hacking). No, there's no example showing that!

  showPlot: show the plot? Default TRUE


VALUE

A vector of "hetzminoo" scores. You can then use it to subset your data by removing samples with unpleasant values.
}")

)

"histoPLOD" <-
structure( function(
  PLODs,
  log= c( FALSE, TRUE)[0],
  mean_show= c( 'POP', 'UP', if( !isFALSE( PLODs@trustplot_FSP_HSP_means)) c( 'HSP', 'FSP')),
  UP_distro_show= c( SPA=TRUE, Normal=FALSE),
  HSP_distro_show= !log && !is.null( PLODs$mean_HSP),
  lb, ub=NULL,
  fullsib_cut= NULL,
  bin= 5,
  main= deparse1( substitute( PLODs), width.cutoff=50),
  ...
){
## Wrapper for pre-existing PLOD_loghisto and HSP_histo,
## which have inconsistent names and arguments
## So don't expect pretty!

stopifnot(
    !missing( log),
    isTRUE( log) || isFALSE( log)
  )

  if( log){
    PLOD_loghisto(
        PLODs,
        UP= 'UP' %in% mean_show,
        HSP= 'HSP' %in% mean_show,
        POP= 'POP' %in% mean_show,
        FSP= 'FSP' %in% mean_show,
        showUP= UP_distro_show,
        main= main,
        ...,
        .deprecate= FALSE
      )
  } else {
    HSP_histo(
        PLODs,
        lb= lb, # missing not allowed
        ub= NULL,
        # fullsib_cut only needed if HSP_distro_show, so NULL def might work
        fullsib_cut= fullsib_cut,
        bin= bin,
        HSPmean = 'HSP' %in% mean_show,
        HSPdist = HSP_distro_show,
        POPmean = 'POP' %in% mean_show,
        FSPmean = 'FSP' %in% mean_show,
        UPmean = 'UP' %in% mean_show,
        main= main,
        ...,
        .deprecate= FALSE
      )
  }
}
, doc =  mvbutils::docattr( r"{
histoPLOD      package:kinference

Histogram PLODs


DESCRIPTION

'histoPLOD' shows a histogram of PLODs across pairwise comparisons, where PLOD has been calculated by 'find_HSPs' (qv), 'find_POPs_lglk' (qv), or something similar. The "histogram" can either be on a log-scale or the more familiar unlogged scale; the problem with the latter is that there are often so many UPs (we have hundreds of millions for tuna) that the interesting kin-bumps become squashed to the point of invisibility. The plots should show various bumps for the different kinships, overlapping to some extent; the theoretical locations of the bump-centres are shown by vertical lines for the most common kinships.

The typical workflow would be to first call with 'log=TRUE', to check that the UP bump looks good and that CK bumps (if any) are centred in the right places. For the UP bump only, the width can also be predicted theoretically, and so the entire predicted curve is shown (and it had better be a good match to the empirical distribution, otherwise there is a QC problem). Then you can focus in on interesting bits and actual pairs, by setting 'log=FALSE' and zooming with the 'lb' and 'ub' parameters.

For 'log=FALSE', the lower bound 'lb' should be set to exclude almost all of the UP bump, which will otherwise swamp the signal from the close-kin pairs. Specifically for HSPs or other 2KPs, it is also possible to show the entire "expected" distribution by setting 'HSP_distro_show=TRUE'. However, unlike the UP bump, its variance and vertical scaling have to be calculated empirically. 'histoPLOD' does that from PLODs between that mean and some upper limit 'fullsib_cut', which must be chosen manually. Hopefully the FSP/HSP gap is clear, so it won't matter much what you choose within that.

The 'kinference-vignette' has examples.


USAGE

histoPLOD( PLODs, log = c(FALSE, TRUE)[0],
    mean_show = c( 'POP', 'UP',
        if( !isFALSE( PLODs@trustplot_FSP_HSP_means)) c( 'HSP', 'FSP')),
    UP_distro_show = c(SPA = TRUE, Normal = FALSE),
    HSP_distro_show = !log && !is.null( PLODs$mean_HSP),
    lb,
    ub = NULL,
    fullsib_cut = NULL,
    bin = 5,
    main = deparse1(substitute(PLODs), width.cutoff = 50),
    ...)


ARGUMENTS

  PLODs: dataframe from 'find_HSPs' or conceivably a future 'find_kin3' etc

  log: TRUE or FALSE to call 'PLOD_loghisto' or 'HSP_histo' respectively

  mean_show: for which kinships should expected PLODs values be shown? There's no harm in showing all of them, but some kinships won't make sense in particular applications, so you can turn them off with this argument. (Note that the default won't show FSP or HSP "means" for 'find_POPs_lglk', since they are calculated wrongly at present.)

  main: graph title, defaulting to the name of 'PLODs' argument

  ...: passed to the plotting routine, which is 'plot' for 'log=TRUE' and 'hist' for 'log=FALSE'.

  UP_distro_show: Only for 'log=TRUE', this controls which approximations to show for the theoretical PLODs distribution of UPs. In practice, they look damn similar! The remaining args only apply when 'log=FALSE':

  HSP_distro_show: Only for 'log=FALSE':show the theoretical PLODs distro for HSPs, based on _empirical_ variance and number of "definite" 2KPs.

  lb, ub: only for 'log=FALSE'. 'lb' is a _mandatory_ cutoff for which PLODs to include. 'ub' is optional, but can be used to visually exclude very high PLODs:eg to exclude POPs when the real interest lies in HSPs.

  fullsib_cut: only if 'log=FALSE' and 'HSP_distro_show=TRUE', then use this to determine which PLODs to include when calculating empirical variance of "HSP PLODs".

  bin: bin width for histogram. Only for 'log=FALSE' since most PLODs (which will be for UPs) are already binned during 'find_HSPs' (qv).
}")

)

"HSP_histo" <-
structure( function(
  kinpairs,
  lb,
  ub = NULL,
  fullsib_cut,
  bin = 5,
  HSPmean = TRUE,
  HSPdist = TRUE,
  POPmean = TRUE,
  FSPmean = TRUE,
  UPmean = FALSE, # normally waaaay off to left
  main = "",
  ...,
  .deprecate= TRUE
){
  if( .deprecate  && are_we_deprecating_yet){
    Deprecate( 'histoPLOD') # unified interface
  }

## MVB added UP option, and tidied a bit in 7/2023, but coulda dun more!
  old_palette <- palette()
  on.exit( palette( old_palette)) # rude to just force it!
  kincols <- kinPalette( setPalette=TRUE)

  PLOD <- kinpairs$PLOD
  if( is.null( ub)){ # default
    ub <- max(PLOD)+bin
  }

  hist.plod <- hist(
    PLOD %such.that% (. %in.range% c( lb, ub)),
    breaks=seq(lb, ub, bin),
    col="lightgrey", xlab="PLOD", main = main, ...)

  # Save lotas hist.plod$... later on;
  extract.named( hist.plod[ cq( mids, breaks, counts)])

  # next is a bit scary; can't we just use kincols[] which has full values?
  # it assumes numbers are right
  kincol <- c( POP=1, HSP=8, FSP=9, UP=5)
  want_means <- sapply( names( kincol),
      function( kinship)
          get( kinship %&% 'mean') && (
          ('mean_' %&% kinship) %in% atts( kinpairs))
    )
  # so want_means[ 'POP'] will contain value of POPmean, etc

  for( kin2show in names( want_means)[ want_means]){
    abline( v=attr( kinpairs, 'mean_' %&% kin2show),
      col=kincol[ kin2show], lwd=2)
  }

  if( HSPdist && want_means[ 'HSP']){
    E.hsp <- kinpairs@mean_HSP
    V.hsp <- mean(sqr(PLOD[PLOD>E.hsp & PLOD < fullsib_cut]-E.hsp))
    obs.num <- counts
    # next should use diff()
    exp.num <- 2*sum(PLOD>E.hsp & PLOD<fullsib_cut)*
        (pnorm( breaks[-1], E.hsp, sqrt(V.hsp))-
         pnorm( breaks[-length(breaks)], E.hsp, sqrt(V.hsp)))
    points( mids, exp.num,pch=16, col=kincol['HSP'], type='b')
  }

  if( any( want_means)){ # MVB version
    legendBits <- kincol[ names( want_means)[ want_means]]
    legend("topright", legend = names( legendBits),
        lwd = 2, lty = 1,
        col = legendBits, bg = "white"
      )
  }


  if( FALSE){ # old code...
    if( HSPmean) { abline(v= kinpairs@mean_HSP, col=8, ,lwd=2) }
    if( POPmean) { abline(v = kinpairs@mean_POP, col=1, lwd=2) }
    if( FSPmean) { abline(v = kinpairs@mean_FSP, col=9, lwd=2) }

    legendBits <- data.frame(
        allNames = c("POP","HSP","FSP"), allNumbers = c(1,8,9))
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
}
, secret_doc =  mvbutils::docattr( r"{
HSP_histo      package:kinference

PLOD histogram


DESCRIPTION

Plots an absolute-frequency histogram for the output of 'find_kinpairs', with the lower bound set by the user. Lower bounds should be set to exclude (as much as possible) the UP bump, as this will otherwise swamp the signal from the HSP bump. Users must manually set a lower bound for full-sibling PLODs ('fullsib_cut') on order to exclude full-siblings from the variance estimate for HSP PLODs.


USAGE

HSP_histo(
  kinpairs,
  lb,
  ub = max(kinpairs$PLOD) + bin,
  fullsib_cut,
  bin = 5,
  HSPmean = TRUE,
  HSPdist = TRUE,
  POPmean = TRUE,
  FSPmean = TRUE,
  main = "",
  ...
)


ARGUMENTS

  kinpairs: the output of a call to 'find_kinpairs'

  lb: PLOD lower bound for plot extent. Should exclude the UP bump

  ub: PLOD upper bound for plot extent. Defaults to maximum PLOD score plus a little padding

  fullsib_cut: PLOD score above which there are only full-sibs

  bin: hist bin width. Default 5. lb, ub, and bin together define 'breaks', so you can't pass 'breaks' via '...'

  HSPmean: plot the mean PLOD for kinpairs? Default TRUE

  HSPdist: plot the distribution of PLOD for kinpairs? Default TRUE

  POPmean: plot the mean PLOD for POPs? Default TRUE

  FSPmean: plot the mean PLOD for FSPs? Default TRUE

  main: graph title, passed straight to 'hist()'

  ...: additional pars, passed to 'hist()'


SEE.ALSO

PLOD_loghisto
}")

)

"HSP_paircomps_lots" <-
function(pair_geno, LOD, geno1, geno2, symmo, eta, min_keep_PLOD, keep_n, minbin, binterval, nbins) {
    .Call(`_kinference_HSP_paircomps_lots`, pair_geno, LOD, geno1, geno2, symmo, eta, min_keep_PLOD, keep_n, minbin, binterval, nbins)
}


"ilglk_geno" <-
structure( function(snpg, hist_pars=list(), showPlot = TRUE) {
  define_genotypes()
  useN <- 'missing' # in case it is!
  extract.named( snpg@locinfo[ cq( pbonzer, useN)])

  p0 <- pbonzer[,'O'] + pbonzer[,'C']
  pA <- pbonzer[,'A']
  pB <- pbonzer[,'B']

  n_samps <- nrow( snpg)
  n_loci <- ncol( snpg)
  minfo_fields <- cq( Our_plate, Our_sample) %that.are.in% names( snpg$info)

  if( my.all.equal( snpg@diplos, genotypes6)) {
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
  } else if( my.all.equal( snpg@diplos, genotypes4_ambig)) {
      snpg4 <- snpg
  } else if( my.all.equal( snpg@diplos, genotypes3)) {
    useN[] <- 3
    # Make a dummy 4-way
    snpg4 <- snpgeno(
      NULL,
      diplos = genotypes4_ambig,
      info= snpg4@info,
      locinfo = snpg4@locinfo
    )
    snpg4[ snpg==AAO] <- AAO
    snpg4[ snpg==BBOO] <- BBO
    snpg4[ snpg==AB] <- AB
  } else {
stop( "diplos( snpg) not one of the known ones")
  }

  snpg3 <- snpgeno(
      NULL,
      diplos = genotypes3,
      info= snpg4@info,
      locinfo = snpg4@locinfo
  )
  snpg3[ snpg4 == OO] <- BBOO
  snpg3[ snpg4==AB] <- AB
  snpg3[ snpg4==AAO] <- AAO
  snpg3[ snpg4==BBO] <- BBOO

  pgeno4 <- matrix( 0, n_loci, 4, dimnames=list( NULL, genotypes4_ambig))
  pgeno4[ , OO] <- sqr( p0)
  pgeno4[ , AB] <- 2*pA*pB
  pgeno4[ , AAO] <- 2*pA*p0 + pA*pA
  pgeno4[ , BBO] <- 2*pB*p0 + pB*pB
  lpgeno4 <- log( pgeno4 + 1e-300) # lambda in doco

  pgeno3 <- matrix( 0, n_loci, 3, dimnames=list( NULL, genotypes3))
  pgeno3[ , AB] <- 2*pA*pB
  pgeno3[ , AAO] <- 2*pA*p0 + pA*pA
  pgeno3[ , BBOO] <- 2*pB*p0 + pB*pB + p0*p0
  lpgeno3 <- log( pgeno3 + 1e-300) # lambda in doco

  ## These are slowish in vecless 1.0
  ## ... but less error-prone to write

  make_SPA_bits <- function( LPGENO){
      if( !nrow( LPGENO)){ # ie no loci with this useN
          K <- dK <- ddK <- function( tt) 0*tt
    return( returnList( K, dK, ddK))
      }


      K <- function( tt) {
          ttp1 <- tt+1
          KK[ j] := SUM_ %[l]% log( SUM_ %[g]% exp( ttp1[ j] * LPGENO[ l, g]))
          return( c( KK))
      }

      dK <- function( tt) {
          ttp1 <- tt+1
          etp1l[ j, l, g] := exp( ttp1[ j] * LPGENO[ l, g])
          num[ j, l] := LPGENO[ l, g] %[g]% etp1l[ j, l, g]
          denom[ j, l] := SUM_ %[g]% etp1l[ j, l, g]
          dKK[ j] := SUM_ %[l]% (num[ j, l] / denom[ j, l])
          return( c( dKK))
      }

      ddK <- function( tt) {
          ttp1 <- tt+1
          etp1l[ j, l, g] := exp( ttp1[ j] * LPGENO[ l, g])
          num[ j, l] := LPGENO[ l, g] %[g]% etp1l[ j, l, g]
          denom[ j, l] := SUM_ %[g]% etp1l[ j, l, g]
          num2[ j, l] := sqr( LPGENO[ l, g]) %[g]% etp1l[ j, l, g]
          dKK[ j] := SUM_ %[l]% ( num2[ j, l] / denom[ j, l] -  sqr( num[ j, l] / denom[ j, l]))
          return( c( dKK))
      }
      K <- compile_vecless( K( -1)) # apparently..!
      dK <- compile_vecless( dK( -1)) # apparently..!
      ddK <- compile_vecless( ddK( -1)) # apparently..!

    returnList( K, dK, ddK)
    }

  SPA_bits3 <- make_SPA_bits( lpgeno3[useN==3, ,drop=FALSE])
  SPA_bits4 <- make_SPA_bits( lpgeno4[useN>3, ,drop=FALSE])
  add34 <- function( Fname) function( tt) SPA_bits3[[ Fname]]( tt) + SPA_bits4[[ Fname]]( tt)

  K  <- add34( 'K')
  dK  <- add34( 'dK')
  ddK  <- add34( 'ddK')

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
## amark

  if( any( useN==3)){
    ilglk <- indiv_lglk_geno(
        lpgeno= lpgeno3[ useN == 3,,  drop = FALSE],
        geno= snpg3[ , useN == 3, drop = FALSE]) ## calc for all individuals, at loci where useN == 3
  } else {  ilglk <- rep( 0, nrow( snpg)) }

  if( any( useN > 3)) {
    ilglk <- ilglk + indiv_lglk_geno(
        lpgeno= lpgeno4[ useN > 3,, drop  =FALSE],
        geno= snpg4[ , useN > 3, drop = FALSE])
  }

  if( FALSE) { # "manual" check on calcs
    lp <- lpgeno[ cbind( rep( 1 %upto% n_loci, n_samps), snpg)]
    dim( lp) <- c( n_loci, n_samps)
    ilglk_manual <- colSums( lp)
  }

  dens_SPA <- renorm_SPA( K, dK, ddK, return_what='func', already_vectorized=TRUE)

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
, doc =  mvbutils::docattr( r"{
ilglk_geno      package:kinference

Check individual multilocus genotypes for typicality


DESCRIPTION

'ilglk_geno' computes the per-sample log-likelihood across its entire genotype, i.e. sum log Pr g(i,l); and compares the distribution across individuals to the theoretical distribution given allele frequencies. Some mismatch is normal (and can arise just from noise in allele-frequency estimates), but substantial mismatch is bad. You get to define "substantial". 'ilglk_geno' can also detect outlier individuals, usually with lglks that are much too low rather than too high; it's not obvious what could generate genomes that are "too typical" at the individual level.

You can use 'locator(1)' to click the histogram to figure out where to adjust the 'xlim/ylim' values to change the range of the data to inspect more closely. In other words, you can then re-run the function with its '...hist_par' argument set accordingly.


USAGE

ilglk_geno(snpg, hist_pars = list(), showPlot = TRUE)


ARGUMENTS

  snpg: a 'snpgeno'. Genotype encoding (see 'diplos') must be one of 4-way, 3-way, or 6-way, as determined by 'diplos' and the '$locinfo$useN' field.

  hist_pars: 'list()' passed to 'hist' for controlling histogram, e.g. 'hist_pars=list(xlim=c(-12000, -6000))', or use 'FALSE' to not plot.

  showPlot: show the histogram? Defaults to TRUE, but overrideen by 'hist_pars=FALSE'.


VALUE

Vector of log-likelihood for each individual; also usually (but optionally), a histogram of log-likelihood values across individuals.


EXAMPLES

## get rid of really bad loci
pvals <- check6and4( bluefin, thresh_pchisq_6and4 = c( 0.001, 0.0001))
bluefin_1 <- bluefin[ , pvals$pval4 > 0.01] # drastic QC!
## check for samples that are not like the others
ilglks <- ilglk_geno( bluefin_1)
## looks like anything with a lglk < -1030 is definitely abnormal
bluefin_2 <- bluefin_1[ ilglks > -1030,]
ilglks <- ilglk_geno( bluefin_2) ## much better, but not perfect -
## see other cleaning steps in the vignette
}")

)

"indiv_lglk_geno" <-
function(lpgeno, geno) {
    .Call(`_kinference_indiv_lglk_geno`, lpgeno, geno)
}


"K_indiv" <-
function(tt, geno, vec_LOD, Pg) {
    .Call(`_kinference_K_indiv`, tt, geno, vec_LOD, Pg)
}


"kin_power" <-
structure( function( lociar,
    want_LOD_table= TRUE, # T/F
    k, # 0.5 for HSPs
    hack_LOD= NULL,
    sd_half_range= 10
){
############
  r"--{
  This comment was originally in the docu, but seems better to put it in the code, to avoid confusing the punters! Anyhow, here it is (andbut it won't be a-happenin' for "kinference 1" AFAICS).

  'kin_power' (and downstream) should get a refactor. It's daft to store LODs for only one specific kin; it'd be better to always calculate P1share and P2share as well as P0share (which is PUP), and then compute whatever-is-needed later on-the-fly. As-is, we are re-computing P1 and P2 based on LOD and PUP OTF instead (which is also unsafe, because LOD could have been calculated with k != 0.5).
  }--"

  define_genotypes()
  li <- lociar$locinfo
  if( my.all.equal( genotypes4_ambig, lociar@diplos)) {
      if(is.null( li$useN)) {
          li$useN <- 4L
      }
  }
  if( all(li$useN == 4) & is.null(li$snerr)) {
      snerrmat <- matrix(0, ncol(lociar), 4,
                         dimnames = list(NULL, c("AA2AO", "AO2AA", "BB2BO", "BO2BB")))
      li$snerr <- snerrmat
  }

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

  s6 <- predict_hsp_util( g6p0, g6p1, g6p2, want_LOD_table, k=k,
      hack_LOD=hack_LOD$LOD6) # $ works on NULLs too

  # For the 4-ways, must condense g6p's
  if( exists( 'genotypes4_ambig', inherits=FALSE)) { # TRUE unless overridden sneakily...
    extract.named( map6to4( g6p0, g6p1, g6p2))
    s4 <- predict_hsp_util( g4p0, g4p1, g4p2, want_LOD_table, k=k,
        hack_LOD=hack_LOD$LOD4) %without.name% "matto"

    ### 3way too:
    extract.named( map6to3( g6p0, g6p1, g6p2))
    s3 <- predict_hsp_util( g3p0, g3p1, g3p2, want_LOD_table, k=k,
        hack_LOD=hack_LOD$LOD3) %without.name% "matto"

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
  lociar@kinPower_checksum <- calc_kinPower_checksum( li)

  # Need to run prepare_PLOD_SPA() again
  lociar@Kenv <- lociar@PLODSPA_checksum <- NULL
  lociar <- prepare_PLOD_SPA(lociar, sd_half_range= sd_half_range)
return( lociar)
}
, doc =  mvbutils::docattr( r"{
kin_power      package:kinference

Locus selection for kin-finding


DESCRIPTION

'kin_power' is an essential preliminary step for various QC and PLOD calculations on real datasets; it calculates summary quantities for each locus based on allele frequency estimates. For real data, the results per se aren't of direct interest to most users. However, as its name suggests, 'kin_power' can also be used to predict how well a proposed set of loci would work for finding 2KPs. (It can also do the same for 3KPs or other specified weaker kin, although the answer will always be "nope" for those cases anyway.)

'kin_power' returns its input 'snpgeno' object after adding extra columns to the 'locinfo' attribute, related to the per-locus mean and variance of LOD (presumably an HSP/UP PLOD, though not inevitably) for different true kinships. It respects the per-locus decision about how precisely to genotype ('useN=6/4/3').


USAGE

kin_power( lociar, want_LOD_table=TRUE, k,
  hack_LOD= NULL, sd_half_range= 10)


ARGUMENTS

  lociar: 'snpgeno' objects with the necessary ingredients

  want_LOD_table: can't think why you'd set this to FALSE

  k: target average kinship for LOD; 0.5 for HSPs, 0.25 for HTPs, etc.

  hack_LOD: Don't mess around with this; it's for internal black magic

  sd_half_range: Normally leave this alone, but in case the final 'prepare_PLOD_SPA' (qv) step gives an error (rare but possible), try reducing it below the default.


VALUE

A 'snpgeno' object with augmented columns in '$locinfo'. Those are:

  E_UP, V_UP: mean & variance for UPs

  E_HSP, E_POP, E_FSP: as you would expect

  Ediff: 'E_HSP-E_POP' ie the "absolute" power of that locus

  sdiff: '(E_HSP-E_POP)/sqrt(V_UP)' which is arguably better than 'Ediff' for ranking loci It also attaches 'LOD', 'PUP', and 'ev01' elements (each a matrix) to '$locinfo'. They have been made dull (see 'mvbutils::make_dull') to improve your viewing experience, but they work fine for all normal purposes (and you can always 'unclass' them to remove the S3 class 'dull').


EXAMPLES

## Next will fail cozza missing pre-calculated objects necessary for kin-finding:
try( hsps <- find_HSPs( dropbears, keep_thresh = 0))
dropbears_1 <- kin_power( dropbears, k = 0.5)
## works now
hsps <- find_HSPs( dropbears_1, keep_thresh = 0)
}")

)

"kinPalette" <-
structure( function( kinships= names( kincolours), setPalette= FALSE){
  kincolours <- c(
    POP= "#0D0887FF",
    GGP= "#48039FFF",
    HCP= "#7401A8FF",
    FCP= "#9D189DFF",
    UP= "#BF3984FF",
    HTP= "#DA596AFF",
    FTP= "#EE7B51FF",
    HSP= "#FBA238FF",
    FSP= "#FCCE25FF")
  if( setPalette){
    palette( kincolours)
  }


return( kincolours[ kinships])
}
, doc =  mvbutils::docattr( r"{
kinPalette      package:kinference

Colors for different kinships


DESCRIPTION

'kinPalette' allows consistent colours for different kinships when plotting; it's used by various 'kinference' plotting functions, and you can also use it yourself to add lines, points, etc. The colours are taken from 'viridisLite::viridis' (see *References*); see the code for the hex values.

'kinPalette' returns a named vector of hex values, either for all kinships or just for those you specify. So you can do something like this:

%%#
# plot something...
abline( v=17, col=kinPalette('UP') # use colour for UPs


.GOOD.MANNERS

Optionally, 'kinPalette' will also call 'grdevices::palette' for you, overwriting the existing numeric color definitions. Afterwards, 'plot(...,col=2)' will show a different color. This isn't normally a good idea (you can always refer to the kin-colors by name; numbers are flaky) so the default is not to (from 'kinference' v1.1 onwards). If you really want to do it yourself, then a perhaps-better approach is:

%%#
kp <- kinPalette()
old_palette <- palette( kp)
# on.exit( old_palette) # if inside a funciotn

so that any palette changes are temporary and can be undone.


USAGE

kinPalette( kinships= names( kincolours), setPalette= FALSE)


ARGUMENTS


  kinships: which kinships to return colors for. Default is all of them.

  setPalette: set to TRUE if you really want to set graphics palette.


VALUE

A character vector of hex codes, with names "POP" etc. If the 'kinships' argument is specified, then just the corresponding elements will be returned.


REFERENCES

S. Garnier, N. Ross, A. Camargo, B. Rudis, K. Woo, & M. Sciaini. (2023). sjmgarnier/viridisLite: CRAN release v0.4.2 (v0.4.2CRAN). Zenodo. https://doi.org/10.5281/zenodo.7890875

(NB I had to abbreviate the authors' first names to avoid stupid unnecessary warning from R CMD CHECK about non-ASCII characters... in a raw string... in the days of mandatory UTF-8... sigh)
}")

)

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


"make_6way_into_4way" <-
function( snpg){
  define_genotypes()
  if( my.all.equal( snpg@diplos, genotypes4_ambig)){
return( snpg)
  }

stopifnot( my.all.equal( snpg@diplos, genotypes6))

  temp_snpg <- snpg
  temp_snpg@diplos <- genotypes4_ambig
  temp_snpg[ snpg==AO] <- AAO
  temp_snpg[ snpg==AA] <- AAO
  temp_snpg[ snpg==BO] <- BBO
  temp_snpg[ snpg==BB] <- BBO
  # Need to do OO & AB too, since encoding differs in 4way vs 6way
  temp_snpg[ snpg==OO] <- OO
  temp_snpg[ snpg==AB] <- AB

return( temp_snpg)
}


"make_ij_character" <-
function( pairs, xsnpg){
  rwid <- xsnpg@info[[ rowid_field( xsnpg)]]
  pairs$i <- rwid[ pairs$i]
  pairs$j <- rwid[ pairs$j]
return( pairs)
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


"opti3ready" <-
structure( function( x2, kin){
  xopti <- kin_power( x2, k=0.25)
  xprod <- x3 <- xopti[1,] # save space

  # xopti has all the LODs and PUPs, but not ev01 needed by autopick
  # Distort and manipulate the LODs... we will need Covar[ LOD2, LOD3]
  for( iuse in names( x2$locinfo) %that.match% '^LOD[0-9]+$'){
    xopti$locinfo[[ iuse]] <- x3$locinfo[[ iuse]] - x2$locinfo[[ iuse]]
    xprod$locinfo[[ iuse]] <- x3$locinfo[[ iuse]] * x2$locinfo[[ iuse]]
  }

  # ev01 needs recalculating for new "LOD"
  # predict_hsp_util does the work, via kin_power
  # k doesn't matter here; Pr[g|coin] is calculated directly
  # We will need E[ PLOD_HSP * PLOD_HTP], from xprod
  ev01_prod <- kin_power( xprod, k=0.5, hack_LOD=xprod$locinfo)$locinfo$ev01
  ev01_opti <- 0 * ev01_prod

  for( coin in cq( 0, 1)){
    eco <- 'e' %&% coin
    e2i <- x2$locinfo$ev01[ , eco]
    e3i <- x3$locinfo$ev01[ , eco]
    eprodi <- ev01_prod[ , eco]

    vco <- 'v' %&% coin
    v2i <- x2$locinfo$ev01[ , vco]
    v3i <- x3$locinfo$ev01[ , vco]

    # E[X-Y] = E[X] - E[Y]
    # V[X-Y] = V[X] + V[Y] - 2C[X,Y]
    # C[X,Y] = E[XY] - E[X]*E[Y]
    ev01_opti[ ,eco] <- e2i - e3i
    ev01_opti[ ,vco] <- v2i + v3i - 2*(eprodi - e2i*e3i)
  }

  xopti$locinfo$ev01 <- ev01_opti

  # Recalculate optiPLODs, just for
  kin <- find_HSPs( xopti, kin[,'i'], kin[,'j']) # not really HSPs...

returnList( xopti, kin)
}
, secret_doc =  mvbutils::docattr( r"{
opti3ready      package:kinference

Prepare for autopick_threshold. DOESN'T WORK YET!!


DESCRIPTION

When choosing threshold to eliminate 3rd-order kin, 'opti3ready' can be used first to recompute an optimal PLOD for HSP::HTP, rather than just using the HSP::UP PLOD. NOT WORKING YET so I've removed it from exports.


USAGE

opti3ready(x2, kin)


ARGUMENTS

  x2: a 'snpgeno'

  kin: result of a previous call to 'find_HSPs' (qv) on 'x2'


DETAILS

Some tricky manoeuvering to calculate the needful...


SEE.ALSO

'doc2Rd', 'flatdoc'


EXAMPLES

# Not compulsory to have an EXAMPLES -- you can put examples into other sections.
# Here's how to make a "don't run" example:
## Not run:

reformat.my.hard.drive()
## End(Not run)
}")

)

"paircomps" <-
function(pair_geno, LOD, geno1, geno2, symmo, granulum, granulum_loci) {
    .Call(`_kinference_paircomps`, pair_geno, LOD, geno1, geno2, symmo, granulum, granulum_loci)
}


"PLOD_loghisto" <-
structure( function(
  hsps,
  UP= TRUE,
  HSP= TRUE,
  POP= TRUE,
  FSP= TRUE,
  showUP= c(SPA= TRUE, Normal= FALSE),
  main= '',
  ...,
  .deprecate= TRUE
){
  if( .deprecate && are_we_deprecating_yet){
    Deprecate( 'histoPLOD') # unified interface
  }

  old_palette <- palette()
  on.exit( palette( old_palette)) # rude to just force it!
  kincols <- kinPalette( setPalette=TRUE)

  ## Far from fully tidied; cf HSP_histo
  binmids <- hsps@bins + (hsps@bins[2] - hsps@bins[1])/2

  ## c++ gives bins that are out in a different direction for bins and n_in_bin
  x <- hsps@bins
  y <- hsps@n_PLODs_in_bin

  l <- list( ...)
  ylim <- l$ylim
  ylim <- if( is.null( ylim)) c( 0, log10(max(y))) else c( 0, ylim[ 2])

  plot( x[-1], log10( head( pmax( y, 0.1), -1)),
     ylim= ylim, main= main, ...,
     type= "S", xlab= "PLOD", ylab= "log10(Frequency)")
  # was: hsps@bins[-1], log10( hsps@n_PLODs_in_bin[1:(length(hsps@n_PLODs_in_bin)-1)]),

  # Pass thru _some_ graphical pars
  # linargs <- list( ...)[ cq( lty, lwd, col)]
  # linargs <- linargs %SUCH.THAT% length(.) # eliminate pure NULLs
  # ... which makes the lines() call look confusing...
  # do.call( 'lines', c( list(
  #     tail( x, -1), log10( y[-1]), type='s'),
  #    linargs))

  # MVB: I can't quite bear the original code here which uses hardcode numbers
  # and repetition! eg
  # if( UP) { abline(v= hsps@mean_UP, col= 5, lwd= 2) }
  # so I tidied a bit...
  shown_kinships <- cq( UP, HSP, POP, FSP)
  shown_kinships <- shown_kinships[
      sapply( shown_kinships, get, envir=environment(), inherits=FALSE)]
  for( kin in shown_kinships){
    abline( v= attr( hsps, 'mean_' %&% kin), col= kincols[ kin], lwd= 2)
  }
  if( showUP["SPA"]) {
      lines( binmids,log10( diff( hsps@binprobs)*sum( y[ binmids<0])),
            lwd=2,col=5)
  }
  if( showUP["Normal"]) {
    lines( x, log10( diff( c( 0,
        pnorm( binmids, mean= hsps@mean_UP, sd= sqrt( hsps@var_UP)) *
        sum( y[ binmids<0])))),
        lwd= 2, col= kincols[ 'UP']) ## Normal approx
  }

  # MVB pu-R-ist recode here ;) to avoid heavvvy data.frame
  legendBits <- kincols[ shown_kinships] # c( UP=5, POP=1, HSP=8, FSP=9)
  legend("topright",
      legend= names( legendBits),
      col= legendBits,
      lwd= 2, lty= 1, bg= "white")

  if( FALSE){ # old code
    legendBits <- data.frame(allNames= c("UP","POP","HSP","FSP"),
        allNumbers= c(5,1,8,9))
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
, secret_doc =  mvbutils::docattr( r"{
PLOD_loghisto      package:kinference

PLOD histogram on log-scale


DESCRIPTION

Plots a log-frequency histogram for the output of 'find_HSPs', with the expected mean PLOD for unrelated pairs, the expected distribution of unrelated pairs, and the expected mean PLOD for HSPs. Colours correspond to 'kinPalette' (qv).


USAGE

PLOD_loghisto(
  hsps,
  UP = TRUE,
  HSP = TRUE,
  POP = TRUE,
  FSP = TRUE,
  showUP = c(SPA = TRUE, Normal = FALSE),
  main= deparse1( substitute( hsps), width_cutoff=50),
  ...
)


ARGUMENTS

  hsps: the output of a call to 'find_HSPs'

  UP, HSP, POP, FSP: whether plot the expected (mean) PLOD for pairs of that type? Defaults TRUE

  showUP: plot the expected density curve for unrelated pairs using the SPA approximation (default TRUE), Normal approximation (default FALSE), both, or neither. Either approximation will plot in colour 5, a light magenta.

  main: optional title for plot

  ...: additional pars, passed to 'plot'

}")

)

"PLOD_oddness_oneway" <-
function(hsps, snpg, lb = min(hsps$PLOD)-10, ub = max(hsps$PLOD)+10, bin = 5,
             CLOD_prop = 0.001, ilglk_prop = 0.001, hetz_prop = 0.001, ...) {

    OLD_DOCO <- r"-----{
#' Oddness metrics
#'
#' Plots the percentage of all pairs in each bin with an unusually low CLOD
#' score, ilglk stat, or hetz stat, across the range of PLOD.
#' \code{PLOD_oddness_oneway} shows the percentage of cases where one member
#' has a low score, and \code{PLOD_oddness_twoway} shows the percentage of
#' cases where both members of the pair have a low score.
#'
#'
#' @aliases PLOD_oddness_oneway PLOD_oddness_twoway
#' @param hsps the output of a call to \code{find_HSPs}
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
#' @param ... additional pars, passed to \code{plot}. \code{ylim} and
#' \code{breaks} are set internally, so you cannot pass them via \code{...}.
#' @noRd
#' @seealso HSP_oddness_oneway
#' @keywords internal
}-----"
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


"predict_hsp_util" <-
function( pIBD0, pIBD1, pIBD2, want_LOD_table=FALSE, k=0.5, hack_LOD=NULL) {
  # This version ignores the possibility of errors involving AB or OO...
  # ... which should be pretty rare

  define_genotypes()
  nl <- nrow( pIBD1)
  Phsp <- pIBD1 * k + pIBD0 * (1-k)  # NOT an HSP unless k=0.5 !!
  Pup <- pIBD0
  Ppop <- pIBD1
  Pfsp <- pIBD0*0.25 + pIBD1*0.5 + pIBD2*0.25

  if( is.null( hack_LOD)){
    LOD <- log( Phsp / Pup)
    LOD[ Pup==0] <- 0 # if Pup=0 then p*log(p) = 0; only happens when r=0
  } else {
    # unpack 2D LOD into 3D; inverse of want_LOD_table really!
    # NB the 2D version of LOD is collapsed (not symmetrical) for space
    # Maybe there is some cleverer way, I dunno
    cn <- cbind(
        c( slice.index( hack_LOD,1)),
        colnames( hack_LOD)[ slice.index( hack_LOD, 2)]
      )
    cn <- sprintf( '%s:%s', cn[,1], cn[,2])
    LOD <- 0 * Phsp # right dims
    dns <- cbind(
        c( slice.index( LOD, 1)),
        dimnames( LOD)[[2]][ slice.index( LOD, 2)],
        dimnames( LOD)[[3]][ slice.index( LOD, 3)]
      )
    dn1 <- sprintf( '%s:%s/%s', dns[,1], dns[,2], dns[,3])
    dn2 <- sprintf( '%s:%s/%s', dns[,1], dns[,3], dns[,2])
    mm <- match( dn1, cn, 0)
    mm[ mm==0] <- match( dn2[ mm==0], cn, 0)
    LOD[] <- hack_LOD[ mm]
  }

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
    # wrong for some reason:
    # double_wanted <- wanted %in% mg[ duplicated( c( mg))]
    gpPUP[ ,double_wanted] <- gpPUP[,double_wanted] * 2

    dimnames( gpLOD) <- dimnames( gpPUP) <- list( dimnames( pIBD0)[[1]], mg@what)
    gpLOD@mg <- mg # why not
  }

  # EPLOD is sum( LOD * Pup) but we want to keep it by locus for now
  # expected value of the PLOD for HSP (mean of distn)
  E_HSP[l] := LOD[l,i,j] %[i,j]% Phsp[l,i,j]
  # expected value of the PLOD for UP (mean of distn)
  E_UP[l] := LOD[l,i,j] %[i,j]% Pup[l,i,j]
  E2_UP[l] := (LOD*LOD)[l,i,j] %[i,j]% Pup[l,i,j]
  V_UP <- E2_UP - sqr( E_UP)
  Ediff <- E_HSP - E_UP

  E_POP[l] := LOD[l,i,j] %[i,j]% Ppop[l,i,j]
  E_FSP[l] := LOD[l,i,j] %[i,j]% Pfsp[l,i,j]

  # In R4.3 I was having printing problems because of retained xtensor class.
  # I've tried to fix that, but just to be on the safe side...
  things <- cq( Ediff, V_UP, E_POP, E_FSP, E_UP, E_HSP)
  for( thing in things){
    assign( thing, unclass( get( thing)), environment())
  }

  ## Code for mean and variance of LOD difference given coinheritance ---
  ## results needed for var.PLOD.kin
  P1 <- 2 * Phsp - Pup # of genopair, given exactly 1 coinherited allele
  e1[ l]:= LOD[l,i,j] %[i,j]% P1[l,i,j]
  e2_1[ l]:= (LOD*LOD)[l,i,j] %[i,j]% P1[l,i,j]
  v1 <- e2_1 - sqr( e1)
  e0 <- E_UP
  v0 <- V_UP
  matto <- cbind( e0, v0, e1, v1)

  # Standardized difference ie locus power: not so useful post hoc,
  #  but possibly interesting for 6 vs 4 comps
  sdiff <- (E_HSP - E_UP) / sqrt( V_UP)

#  Ediff <- unclass( Ediff)
#  V_UP <- unclass( V_UP)
#  sdiff <- unclass( sdiff)

  retval <- data.frame( Ediff, V_UP, sdiff, matto, E_POP, E_FSP, E_UP, E_HSP)
    ## matto comes in as 4 named columns, not as matto.
    ## E_UP is not 100% efficient, but bonus points for readability
  if( want_LOD_table) {
    retval@LOD <- gpLOD
    retval@PUP <- gpPUP
  }
return( retval)
}


"prepare_PLOD_SPA" <-
structure( function( geno6, n_pts_SPA_renorm=201, sd_half_range=10) {
    ## To be run after kin_power( ..., want_LOD_table=TRUE)

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
      e$sd_half_range <- sd_half_range

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

          extract.named( renorm_SPA_cumul( K, dK, ddK, n_pts=n_pts_SPA_renorm, sd_half_range=sd_half_range))
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
, doc =  mvbutils::docattr( r"{
prepare_PLOD_SPA      package:kinference

Prepare for kin-finding


DESCRIPTION

'prepare_PLOD_SPA' is something you used to have to run before using some kin-finding/QC tools, to set up your 'snpgeno' object for fancy maths woooo (saddlepoint approximations). There are no meaningful options, you just have to run this. It can be _slightly_ slow which is why it was a separate step. However, nowadays I don't think you need to run it at all, because it's built into 'kin_power' (qv).


USAGE

prepare_PLOD_SPA(geno6, n_pts_SPA_renorm= 201, sd_half_range= 10)


ARGUMENTS

  geno6: a 'snpgeno' object that has been thru 'kin_power'

  n_pts_SPA_renorm: how accurate to make the approximation. Default should be fine.

  sd_half_range: How many SD's into the tails to push the approximation. The default of 10 is massively far. Normally this is fine, but if you get an error (probably from an NA cropping up during extreme calculations), then trying making it smaller.


VALUE

Another 'sngeno' object with an environment 'Kenv', which contains functions (with their own preloaded data) allowing null distributions (eg PLODs for true UPs) to be calculated. Various sanity checks are incorporated to try to stop you from stuffing up with out-of-synch loci etc later.
}")

)

"re_est_ALF" <-
function( snpg) {
## check to be called after load_whopper loads entire dataset
# Pretty weird, since it expects 6-way input, but blurs it to (in effect) 4-way!
# Why not use est_ALF_6way() instead?

  define_genotypes()
  n_samp <- nrow( snpg)
  n_loci <- ncol( snpg)
  gamb <- matrix( character(), n_samp, n_loci)
  gamb@diplos <- genotypes_ambig
  if( my.all.equal( diplos( snpg), genotypes6)){
    gamb[ snpg==AB] <- AB
    gamb[ snpg==OO] <- OO
    gamb[ snpg==AO] <- AAO
    gamb[ snpg==AA] <- AAO
    gamb[ snpg==BO] <- BBO
    gamb[ snpg==BB] <- BBO
  } else if( my.all.equal( diplos( snpg), genotypes4_ambig)){
    # gamb should be ABCO but snpg is ABO so can't rely on "obvious" mapping
    gamb[ snpg==AB] <- AB
    gamb[ snpg==OO] <- OO
    gamb[ snpg==AAO] <- AAO
    gamb[ snpg==BBO] <- AAO
  } else if( !my.all.equal( diplos( snpg), genotypes_C)) {
stop( "Genotypes must either be 6-way, 4-way ABO, or 4-way ABCO ie potentially with 3rd allele")
  }

  ## snpg@geno_amb <- gamb # required by...
  snpg <- est_ALF_ABCO( snpg, geno_amb = gamb)
  snpg$locinfo$pbonzer <- snpg$locinfo$pambig
  snpg$locinfo$pambig <- NULL
return( snpg)
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
}


"split_FSPs_from_HSPs" <-
structure( function( snpg, candipairs) {
  # For pairs already picked as HSPs, ie PLOD(HSP,UP) > eta: they might be FSPs

  # Don't need full pairwise screening for FSPs (do post hoc on a few hundred
  # HSPs), hence all in R.

  define_genotypes()

  # HSPs normally from 'find_HSPs'; or can be M*2 matrix of rows in snpg that are poss HSPs
  # if former, make latter

  if( candipairs %is.a% 'data.frame') {
    candipairs <- as.matrix( candipairs[ cq( i, j)])
  }
  sibg <- just_sibg <- snpg[ c( candipairs),]

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

  nsib <- nrow( candipairs)
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
      p12fsp[l, i] <-
          kappa_fsp[1] * P_k0[l, g1[i, l], g2[i, l]] +
          kappa_fsp[2] * P_k1[l, g1[i, l], g2[i, l]] +
          kappa_fsp[3] * P_k2[l, g1[i, l], g2[i, l]]
      # P[g_1l g_2l | HSP]
      p12hsp[l, i] <-
          kappa_hsp[1] * P_k0[l, g1[i, l], g2[i, l]] +
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
  ret <- data.frame(
      PLOD_FH= PLOD_FH,
      i= candipairs[,1],
      j= candipairs[,2]
    )

  # Next 2 are COMPLETELY WRONG !!!
##  ret@E_FSP <- EPLOD_FH_F
##  ret@E_HSP <- EPLOD_FH_H

  ret@E_FSP <- EPLOD_FH_F
  ret@E_HSP <- EPLOD_FH_H

  ret@call <- sys.call()

  return(ret)
}
, doc =  mvbutils::docattr( r"{
split_FSPs_from_HSPs      package:kinference
split_FSPs_from_POPs
split_HSPs_from_HTPs

Discriminate between kinships of known close-kin


DESCRIPTION

These are for pairs already picked as likely close-kin via one of the 'find_XXX' functions, but whose exact kinship is uncertain; e.g., they might clearly be either POPs or FSPs, but it's not obvious which. The 'split_XXX_from_YYY' functions apply a more powerful likelihood-based test statistic to each pair, to help decide what it is. All these functions use 4-way genotypes, i.e. not relying on 6-way genotyping. They probably should be adapted to cope with 3-way (ie not trusting double-nulls) but currently they aren't (so they do trust double-nulls).


USAGE

split_FSPs_from_HSPs(snpg, candipairs)
split_FSPs_from_POPs( snpg, candipairs, gerr, use_obsolete_version=FALSE)
split_HSPs_from_HTPs(snpg, candipairs)


ARGUMENTS

  snpg: a 'snpgeno' object

  candipairs: normally, a dataframe with rows being pairs and columns _i_ and _j_ (and possibly others) e.g. from 'find_POPs' or 'find_HSPs'. Can also be a 2-column matrix (each row again one pair).

  gerr: genotyping error rate, where 0.01 would mean '1%'. It's there to make the POP case robust; you have to choose it, but the precise value should not matter. See *Details*.

  use_obsolete_version: the original code for POPs-vs-FSPs was based on a different non-likelihood-based statistic. It turned out to have low statistical power, but did not require specifying a 'gerr'. For replicability purposes, you can still run it by setting this parameter to 'TRUE', but otherwise don't bother.


DETAILS

The idea of 'split_FSPs_from_POPs'- though this is not the only possible workflow- is that pairs which are _either_ POPs _or_ FSPs should stand out very clearly from everything else, via 'find_POPs'. Then the job is to pick between those possibilities. The workflow is supposed to be:

  nail POPs/FSPs first with 'find_POPs';

  pick between them with 'split_FSPs_from_POPs', making use of eg age data too;

  look for HSPs (and potentially some HTPs) and filter out already-known POPs and FSPs;

  filter out HTPs from the remaining set of HSPs with 'split_HSPs_from_HTPs' and/or 'autopick_threshold'.

However, an equally reasonable workflow might be:

  nail HSPs and everything stronger (and potentially some HTPs) with 'find_HSPs';

  split HSPs/HTPs from POPs/FSPs with 'split_FSPs_from_HSPs';

  filter out HTPs from the remaining set of HSPs with 'split_HSPs_from_HTPs' and/or 'autopick_threshold';

  use 'split_FSPs_from_POPs' to do just what it says. (Again, age data may help in marginal cases.)

All 'split_<blah>' functions return expected values under different possible kin-types (but not variances, since these cannot be predicted for all kin-types).

The 'gerr' parameter in 'split_FSPs_from_POPs' is there to alleviate the problem that a single locus displaying apparent Mendelian exclusion is in theory reason enough to prove that a pair is _not_ a POP (if a likelihood-based criterion is used). But, of course, we can have genotyping errors (and, in a reasonably big dataset, mutations). Allowing for a small amount of error gives the method much mor flexibility, without paying a high price in statistical efficiency (provided 'gerr' is small). The technical interpretation is that, if a genotyping error occurs at a locus, then the true genotype is replaced by a randomly-drawn genotype from the marginal distro of genotypes at that locus. Real genotyping errors don't work like that, but it is mathematically convenient and achieves the desired effect of robustifying the FSP-vs-POP statistic. The value to use is up to you; you can experiment; if you really want to estimate it, then look at replicate genotypes.


VALUE

A dataframe with columns 'PLOD_FH' (or similar depending on what is being split), 'i', and 'j', containing the indices of each pair. Theoretical means for the two kinships in the function's name are returned as attributes (variances cannot be predicted), as is the 'call'.


EXAMPLES

library( atease) # so x@att--- easier than attr( x, 'att')
## bluefin data
## stripped-down data-cleaning for example - see the
## vignette for approach with real data!
pvals <- check6and4( bluefin, thresh_pchisq_6and4 = c( 0.001, 0.0001))
bluefin_1 <- bluefin[ , pvals$pval4 > 0.01]
ilglks <- ilglk_geno( bluefin_1)
bluefin_2 <- bluefin_1[ ilglks > -1050,]
bluefin_3 <- est_ALF_ABO_quick( bluefin_2)
bluefin_4 <- bluefin_3[ , bluefin_3@locinfo$pbonzer[,"B"] > 0.02]
bluefin_5 <- kin_power( bluefin_4, k = 0.5)
dups <- find_duplicates( bluefin_5, max_diff_loci = 20)
bluefin_6 <- bluefin_5[ -c(drop_dups_pairwise_equiv( dups)) ,]
## via find_HSPs for isolated split of 'FSPsOrPOPs'
hsps <- find_HSPs( bluefin_6, keep_thresh = 0)
histoPLOD( hsps, log= FALSE, lb= 0, fullsib_cut= 70, HSP_distro_show= TRUE)
## note gap between PLOD  60 -- 70, separating likely HSPs
## from likely POPs and FSPs
FSPsOrPOPs <- hsps[ hsps$PLOD > 70,]
HSPsOrHTPs <- hsps[ hsps$PLOD < 70,]
## ... implicitly >0 from keep_thresh in find_HSPs() call
## via find_POPs_lglk for 'laser-focus on POPs'
maybePOPs <- find_POPs_lglk( bluefin_6, gerr = 0.01, keep_thresh = -30)
hist( maybePOPs$PLOD, breaks = 20)
abline( v = maybePOPs@mean_POP, col = kinPalette("POP"), lwd = 2)
abline( v = maybePOPs@mean_UP, col = kinPalette("UP"), lwd = 2)
## note the 6 pairs clustered around the expected mean PLOD for
## POPs (all > 60). They are probably POPs, but may include FSPs.
maybePOPs <- maybePOPs[ maybePOPs$PLOD > 60,]
## split_FSPs_from_POPs
### using FSPsOrPOPs
splitFSPsOrPOPs <- split_FSPs_from_POPs( bluefin_6,
    candipairs= FSPsOrPOPs, gerr = 0.01)
hist(splitFSPsOrPOPs$PLOD)
## add line of expected PLOD if FSP
abline( v = splitFSPsOrPOPs@E_FSP, lwd = 2, col = kinPalette("FSP"))
## add line of expected PLOD if POP
abline( v = splitFSPsOrPOPs@E_POP, lwd = 2, col = kinPalette("POP"))
## Expected PLOD is -ve for true POPs, +ve for true FSPs
### using maybePOPs
splitmaybePOPs <- split_FSPs_from_POPs( bluefin_6,
    candipairs= maybePOPs, gerr = 0.01)
hist(splitmaybePOPs$PLOD)
## add line of expected PLOD if FSP
abline( v = splitmaybePOPs@E_FSP, lwd = 2, col = kinPalette("FSP"))
## add line of expected PLOD if POP
abline( v = splitmaybePOPs@E_POP, lwd = 2, col = kinPalette("POP"))
## Expected PLOD is -ve for true POPs, +ve for true FSPs
## split_HSPs_from_HTPs
splitHSPsOrHTPs <- split_HSPs_from_HTPs( bluefin_6, HSPsOrHTPs)
hist(splitHSPsOrHTPs$PLOD, breaks = 15)
## add line of expected PLOD if HSP
abline( v = splitHSPsOrHTPs@E_HSP, lwd = 2, col = kinPalette("HSP"))
## add line of expected PLOD if HTP
abline( v = splitHSPsOrHTPs@E_HTP, lwd = 2, col = kinPalette("HTP"))
## Expected PLOD is +ve for true HSPs, -ve for true HTPs
}")

)

"split_FSPs_from_POPs" <-
function(
    snpg,
    candipairs,
    gerr,
    use_obsolete_version=FALSE
){
## For pairs already picked as first-order kin (fsp or pops), eg via plod(hsp,up) > eta
## don't need full pairwise screening (do post hoc on a few hundred candidates), hence all in r.
## Based almost entirely on code from split_fsps_from_hsps(), just changed at the end...

  if( use_obsolete_version){
    # Non-lglk, based on ppn IBS *genotypes*
return( OLD_split_FSPs_from_POPs( snpg, candipairs))
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

  if( candipairs %is.a% 'data.frame') {
    candipairs <- as.matrix( candipairs[ cq( i, j)])
  }
  sibg <- just_sibg <- snpg[ c( candipairs),]

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


  nsib <- nrow( candipairs)
  nloci <- ncol(sibg)

  # if only R were zero-indexed... (SARCASM, PRESUMABLY!)
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
  ret <- data.frame(
      PLOD_FP= PLOD_FP,
      i= candipairs[,1],
      j= candipairs[,2])

## Comment from FSP/HSP version was: Next 2 are COMPLETELY WRONG !!!
## ??? Is that true..??? They look OK to me!

  ret@E_FSP <- EPLOD_FP_F
  ret@E_POP <- EPLOD_FP_P

  ret@call <- mvb.sys.call(1) # don't ask, just works

return( ret)
}


"split_HSPs_from_HTPs" <-
function( snpg, candipairs) {
  # For pairs already picked as possible HSPs, they might be HTPs

  # Don't need full pairwise screening for FSPs (do post hoc on a few hundred
  # HSPs), hence all in R.

  define_genotypes()

  # HSPs normally from 'find_HSPs'; or can be M*2 matrix of rows in snpg that are poss HSPs
  # if former, make latter

  if( candipairs %is.a% 'data.frame') {
    candipairs <- as.matrix( candipairs[ cq( i, j)])
  }
  sibg <- just_sibg <- snpg[ c( candipairs),]

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

  nsib <- nrow( candipairs)
  nloci <- ncol(sibg)

  # if only R were zero-indexed... HA! NO WAY! I *did* spot this comment
  # 0-based is madness; if there has to be a single base, it must be 1 !
  # But of course R should let indices start
  # at _any_ value, at user's discretion... hence
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
                    i       = candipairs[,1],
                    j       = candipairs[,2])

  # Next 2 are COMPLETELY WRONG !!!
##  ret@E_HSP <- EPLOD_ST_F
##  ret@E_HTP <- EPLOD_ST_H

  ret@E_HSP <- EPLOD_ST_HS
  ret@E_HTP <- EPLOD_ST_HT

  ret@call <- sys.call()

  return(ret)
}


"upairid" <-
structure( function( i, j, sep='/'){
  if( missing( j) && (i %is.a% 'list')){ # called direct
    j <- i$j
    i <- i$i
  }
  paste( pmin( i, j), pmax( i, j), sep=sep)
}
, doc =  mvbutils::docattr( r"{
upairid      package:kinference

Label kin-pairs


DESCRIPTION

'find_HSPs' etc return _pairs_ of kin, labelled by 'i' and 'j'. If you want to see whether one or more pairs also occur in another set of kin, eg from calling 'find_POPs_lglk', it's tempting to 'paste' together each set's 'i' and 'j', and then use 'match'. But that will go fail if the order has switched, so that 'i' in the first set corresponds to 'j' in the second, and so on- which can definitely happen with 'find_XXX' (for good reasons). To avoid that, 'upairid' sorts each (i,j) pair into order before pasting them, so that you can use 'match( upairid( first_set), upairid( second_set))' with no need to worry about order-switching.


USAGE

upairid(i, j, sep='/')


ARGUMENTS

  i, j: If 'i' is a list (presumably from a call to 'find_XXX') and 'j' is missing, then i. Otherwise, 'i' and 'j' can be vectors of numbers or characters.

  sep: character (or string) to put between each 'i' and 'j'.


VALUE

A character vector of pairwise-sorted, pasted, IDs.


EXAMPLES

upairid( 1:5, 5:1)
d1 <- list( i=1:3, j=4:6)
d2 <- list( i=4:6, j=1:3) # swapped
match( paste( d1$i, d1$j), paste( d2$i, d2$j)) # nope
match( upairid( d1), upairid( d2)) # yep
}")

)

"var_PLOD_kin" <-
structure( function(
    linfo,
    emp_V_HSP= V_noX( C_equiv, 2),
    n_meio,
    debug=FALSE,
    C_equiv=NULL ) {
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
  # as produced by kin_power
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

  count_l <- linfo$count # 's only present if linfo is simulated
  # eg 100 loci like this one, 100 like the next...
  if( is.null( count_l)){
    count_l <- rep( 1, nrow( linfo))
  }

  L <- sum( count_l)
  pi <- count_l / L # sum(pi)==1

  # Extract columns of ev01. Slightly odd code--- maybe don't need to bother
  ev01 <- unclass( linfo$ev01) # get rid of "dull" attr
  for( thing in THINGS ) {
    assign(thing %&% "_l", linfo$ev01[,thing])
    # ... the "_l" is a pseudo subscript...
  }

  # What are the "typical" properties of a locus?
  e0 <- pi %**% e0_l
  e1 <- pi %**% e1_l
  v0 <- pi %**% v0_l
  v1 <- pi %**% v1_l

  # noX = no crossover, i_e. all in separate equal chromos
  V_noX <- function( C, meioses) {
      p <- 2 ^ (1-meioses) # Marginal prob of coin = 1/2 @ HSPs, 1/8 @ HCPs
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
        find.root( V_noX, target=emp_V_HSP, start=1, step=1,
            fdirection='decreasing', min.x=1, max.x=L, meioses=2)
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
    # Assumes we are HSP or descendent...
    # ... ie stepladder; not extension ladder a la GGP

    pHSP_11 <- (1/4) * (1+exp( -4*rho*dee))

    # Extrapolate to number of meioses, m rather than HSP
    pm_11 <- ((0.25 * (1 + exp( -2*rho*dee))) ^ (meioses-2L)) * pHSP_11
    pm <- 2^(1L-meioses) # AKA pm_1: marginal prob of coin at a single locus

    pm_10 <- pm - pm_11
    # pm_01 <- pm_10 is just symmetry
    pm_00 <- 1 - pm_11 - 2*pm_10

    EL <- L * ( e1_1*pm + e1_0*(1-pm))
    EL2 <- L * ( e2_1*pm + e2_0*(1-pm)) +
        2 * (L-dee) %**%
        ( sqr( e1_1)*pm_11 + 2*e1_1*e1_0*pm_10 +  sqr( e1_0)*pm_00)
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
, doc =  mvbutils::docattr( r"{
var_PLOD_kin      package:kinference

Predict variance of PLOD for HTPs and HCPs


DESCRIPTION

'var_PLOD_kin' infers the extent of genetic linkage based on the empirical (i.e. observed) variance of PLODs among (presumed) true HSPs, by comparing it to the (smaller) variance that would theoretically apply if all loci were inherited independently. It then predicts the empirical variance for HTPs and/or HC1Ps, based on one of two extreme assumptions about the detailed nature of recombination. The truth will likely be between those extremes. 'var_PLOD_kin' is called internally by 'autopick_threshold' (qv), and few users will want to call it directly.

Although the biggest linkage effect comes from the main "order" of kinship (2KP, 3KP, etc, as controlled by the 'n_meio' argument), this is one situation that is also affected somewhat by the "subtype" of kinship (HSP vs GGP vs FTP, etc). For CKMR purposes, the commonest type of kin of given order are typically those born closest in time, so for a given order of kinship, the 'var_PLOD_kin' always assumes that the relevant kin-pairs will consist of the subtype with _single_ shared ancestor and minimax number of generations since shared ancestor. This means HSPs for 'n_meio=2' (FTPs have 2 shared ancestors; GGPs entail 2 generations of gap, whereas HSPs have only 1), HTPs for 'n_meio=3', HC1Ps for 'n_meio=4', etc.

'var_PLOD_kin' can also be used for "simulation", e.g. to explore the effects of different types of linkage, or of different numbers of markers; see the 'linfo' and 'C_equiv' arguments, plus the code itself, for hints.


.NATURE.OF.LINKAGE

The two extreme alternatives considered by 'var_PLOD_kin' each use one parameter to capture "linkage", but the parameter has quite different meanings. The two versions are:

  separate chromosomes, each containing the same number of loci, with _no_ recombination. The linkage parameter is the number of such chromosomes.

  one single chromosome containing all the loci, equally spaced, with random recombination at multiple places. The linkaage parameter is the probability of recombination between adjacent loci.


USAGE

var_PLOD_kin(
  linfo,
  emp_V_HSP = V_noX(C_equiv, 2),
  n_meio,
  debug = FALSE,
  C_equiv = NULL
)


ARGUMENTS

  linfo: either a 'snpgeno' object, or its '$locinfo' attribute (or a fake one). Respectively, either 'linfo$locinfo' or 'linfo' itself should be a dataframe containing columns 'e0', 'e1', 'v0', 'v1', 'count'. Each row is one "type" of locus, i.e., with roughly the same values of e/v 0/1, and 'count' says how many such loci there are (for simulation purposes; with real data, 'count==1'). ' e/v 0/1' are means and variances of the per-locus LOD (note no P) when the locus is or isn't co-inherited. See _Details_

  emp_V_HSP: empirical variance of PLOD among definite HSPs, normally from 'autopick_threshold' (qv). You're supposed to be running this on real data, so that 'emp_V_HSP' is an actual number; however, for testing purposes, you can set up an artificial version via 'C_equiv' below.

  n_meio: Target number of meioses:2 for HSPs, 3 for HTPs, 4 for HCPs, etc. This is by far the main driver of variance, but technically the only one; see _Subtypes of kin_.

  debug: Logical flag, not for regular folks to use! Defaults to FALSE. The 'debug' package needs to be loaded beforehand.

  C_equiv: for artificial test, with 'emp_V_HSP' set to the no-crossover variance from 'C_equiv' chromos (need not be integer). Ignored if 'emp_V_HSP' is set.


VALUE

Matrix with two rows 'V0' and 'Vx', and one column for each element of 'n_meio' (which is always augmented to include 2), named "M2" etc. The two rows pertain respectively to the no-crossover multiple-chromosome scenario, and the single-chromosome multiple-crossover scenario. The matrix also has an attribute 'info', which is a numeric vector of elements named 'V_UP', 'V_HSP', 'C_hat', and 'rho_hat' (note that 'V_HSP' should duplicate the first column of the matrix) 'C_hat' is estimated equivalent number of chromosomes for the no-crossover scenario, and 'rho_hat' is per-locus crossover rate for the all-crossover scenario.


GRUESOME.DETAILS

The "per-locus LOD" (whose properties are stored in the columns 'e0', 'e1', 'v0', 'v1' in 'linfo') is created by calling 'kin_power' (qv). The normal use-case would be from 'kin_power(...,k=0.5)', so that the (P)LOD pertains to 2KP/UP comparisons. However, if you called 'kin_power(...,k=0.25)' then the (P)LOD would be designed for 3KP/UP comparisons, and so on. In fact, you could even hand-tweak the calculations to contain LODs for 3KP/2KP comparisons, which _might_ in principle improve the resolution (but you'd have to fiddle manually; you could actually do it based on two calls to 'kin_power', one with 'k=0.5' and one with 'k=0.25', and manipulating the results).


EXAMPLES

# COMPLETELY MADE-UP e/v values! Nothing to do with genetics :)
var_PLOD_kin( data.frame( count=45, ev01= I( cbind( e0=-1, e1=2, v0=0.03, v1=0.02))), C_equiv=22, n_meio=3:4)
#        M2    M3    M4
#  V0 208.2 156.6  91.9
#  Vx 208.2 186.2 101.8
#  attr(,"info")
#      V_UP    V_HSP    C_hat  rho_hat
#    1.3500 208.2273  22.0000   0.2616
}")

)
