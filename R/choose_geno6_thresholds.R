## TODO
# - document return items
#' Genotyping pre-processed DartCap
#'
#' @aliases geno6way
#'
#' \code{choose_geno6_thresholds} is to be run on a generous subset of the data, to figure out the thresholds for single nulls etc etc. Then, actually \bold{doing} the genotyping is handled by \code{geno6way}. You can call the latter separately, so that you can apply thresholds etc derived from one subset of the data to a new dataset, so that your genotyping is consistent.
#'
#' Classifies each fish/locus (a biallelic SNP perhaps with nulls) into one of six categories: OO, AO, BO, AA, BB, AB. Input should come from previous steps of \code{check_baits}, either "stand-alone" to compute thresholds and error-rates, or applying preset thresholds-and-error-rates to fresh count data that's had some processing. \code{geno6} \emph{discards} the raw counts and returns instead a \code{snpgeno} object with associated locus information.
#'
#' If in stand-alone mode and assuming \code{plots=TRUE}, then fits to each locus are shown graphically, and you get the chance to make a comment; any non-null comments are recorded to a file. The idea is that you can look at the file afterwards to see which loci you really didn't like.
#'
#' Can also be run to "shut up and genotype" a dataset using previously-determined thresholds (eg if you have new plates of fish).
#'
#' @usage geno6(lociar, OOthresh_tc = lociar@args$geno_deambig$OOthresh_tc,  nquants_bump, max_dat_quantile, nhistbreaks = 101, distro,  minpO_rethresh, m_for_rethresh, max_refits, plots = NULL, minbin, show_progress = interactive())
#'
#' @param lociar a \code{loc.ar} from previous calls to \code{check_baits}, in p'tic with attributes
#' @param OOthresh_tc ???
#' @param nquants_bump ???
#' @param max_dat_quantile ???
#' @param nhistbreaks ???
#' @param distro ???
#' @param minpO_rethresh ???
#' @param m_for_rethresh ???
#' @param max_refits ???
#' @param plots ???
#' @param minbin ???
#' @param show_progress ???
#'
#' @return A \code{snpgeno} object with the original the \code{locinfo} attribute augmented by the :
#' \item{ snerr:}{}
#' \item{ best_cut:}{}
#' \item{ pbonzer:}{}
#' \item{ perr:}{}
#'
#' Not compulsory. Other section headings, e.g. AUTHOR, should also go here. Use \bold{single} quotes around object names and code fragments, e.g. \code{bit.of.code()}. Use \bold{double} quotes for "text" or "file.name". See \code{\link{doc2Rd}} for full details of format.
#'
#' @importFrom atease @
#' @importFrom mvbutils cq my.all.equal %in.range% scatn %&% FOR
#' @importFrom graphics par
#' @importFrom vecless :=
#' @importFrom stats dnorm pnorm pt dt quantile qnorm nlminb
#' @importFrom utils flush.console tail
#' @importFrom handy2 sqr find.root rel.delta
#' @importFrom gbasics logit
#' @export
"choose_geno6_thresholds" <- function( lociar,
    li= NULL, # could supply presets here
    OOthresh_tc= lociar@args$geno_deambig$OOthresh_tc,
    nquants_bump, # 10,
    max_dat_quantile, # 0.95,
    nhistbreaks=101,
    distro, # 'normal' or 't6' etc
    minpO_rethresh, # 0.05 since 0.05 * 0.05 = bugger all
    m_for_rethresh, # eg 5 means: split widest gap containing 5 fish
    max_refits, # at least 1, maybe 2 to shift the upper-tail outlier-trap "hix"
    minbin, # only for plots
    show_progress= interactive(),
    plots=TRUE,
    comment_file= NULL){

  # Updated 2016 version for check_baits() pipeline
  # Assumes ambig genotypes (AAO = AA|AO) with poss 3rd allele (C), have
  # been set earlier on plotting controlled by "..._args"

  n_loci <- ncol( lociar)
  n_fish <- nrow( lociar)

  define_genotypes()

  dimnames( lociar)[[3]] <- cq( A, B, C)

  if( is.null( li)) { # use own data
    li <- lociar@locinfo
  } else { # check preset definitions match this OK...
    mm <- match( lociar@locinfo$consensus, li$consensus, 0)
    stopifnot( all( mm>0))
    li <- li[ mm,]
    stopifnot( all( lociar@locinfo$FullAltSeq==li$FullAltSeq))
    stopifnot( all( lociar@locinfo$FullRefSeq==li$FullRefSeq))
    # NA-friendly
    stopifnot( my.all.equal( lociar@locinfo$FullThirdSeq, li$FullThirdSeq))
  }

  c3 <- unclass( lociar)
  locus_name <- li$Locus
  mht <- li$med_het_tot
  pambig <- li$pambig
  if( is.null( mht)) {
    stop( "Need med_het_tot column in locinfo/li")
  }

  # counts normed by locus
  c3[f,l,all] := c3[f,l,all] / mht[ l]
  cA <- c3[,,A]
  cB <- c3[,,B]
  cC <- c3[,,C]
  # guard against length-1 dims cRazies
  dim( cA) <- dim( cB) <- dim( cC) <- dim( c3)[1:2]

  geno <- lociar@geno_amb
  gobs <- matrix( 0, n_loci, length( genotypes_ambig),
                 dimnames=list( NULL, genotypes_ambig))
  for( ig in genotypes_ambig) {
    gobs[,ig] <- colSums( geno==ig)
  }

  dname <- paste(deparse(substitute(lociar), width.cutoff=30, nlines=1),
                 collapse='')

  opar <- par(mar=c( 2, 3, 2, 0)+0.1, no.readonly=TRUE)

  allres <- NULL
  comments <- structure( character( length( mht)), names=locus_name) # for graphics "annotation"

  if( tolower( distro)=='normal') {
    CDF <- pnorm
    PDF <- dnorm
  } else {
    tdof <- as.numeric( sub( 't', '', tolower( distro)))
    CDF <- function( x, mean, sd, lower) pt( (x-mean)/sd, df=tdof, lower.tail=lower)
    PDF <- function( x, mean, sd, lower) (1/sd) * dt( (x-mean)/sd, df=tdof)
  }

  # Function needed later: optimal split is when tail-probs are equal,
  # so their diff should be zero
  diffprob <- function( x, pX, pO) {
    # num & denom should be mult by pX, but this simpler form avoids rounding error
    cond_pXX <- pX / (pX + 2*pO)
    cond_pXX * CDF( x, mean=mu[2], sd=sigma[2], lower=TRUE) - 
      (1-cond_pXX) * CDF( x, mean=mu[1], sd=sigma[1], lower=FALSE)
  }

  alid <- li$AlleleID
  loclabs <- NULL

  if( !show_progress) {
    tt <- tempfile()
    sink( tt)
    on.exit( {sink(); unlink( tt)})
  }

  # More careful about totals, to zap splatter from 3rd allele
  tot2 <- make_tot2( cA, cB, cC, geno)

  pbonzer <- structure( rep( 0, 4), names=cq( A, B, C, O))

  new_OOthresh <- rep( OOthresh_tc, n_loci)
  for( i in 1:n_loci) {
    if( show_progress) {
      cat( '\r', i); flush.console()
    }

    dat <- unname( tot2[,i])

    # as.character() needed so that match() works
    genoi <- as.character( geno[,i])
    is_OO <- genoi==OO
    is_hetz <- genoi %in% c( AB, BC, AC)

    bumpiq <- seq( 0, 1, length=nquants_bump+1)
    qdat <- quantile( dat[ !is_OO], bumpiq)
    dqd <- diff( qdat)
    big_bumpi <- which.min( dqd)
    qdat <- c( qdat, tail( qdat, 1)) # to avoid OOR on next step
    big_bump <- mean( qdat[ big_bumpi+0:1])


    # Outlier protection: counts above hix get treated as just "above hix"
    # hix SHOULD be substantially more than mean of higher bump,
    # otherwise silly... so, make sure it is
    # Rarely, can end up too low after fitting bump properly, in which case
    # it will be upped and then a refit
    hix <- quantile( dat[ is_hetz], max_dat_quantile)
    # arbitrary; for mean counts of 10000000 it could be too large...
    #  will refit if required
    hix <- max( hix, big_bump * 1.8)
    setup_counts_and_nhi()

    # Bumps are where quantiles are closest; max 2 bumps in theory
    # Biggest bump is-- presumably-- either single-nulls or no-nulls;
    # try both poss
    # Don't use too many quantiles... 10 might be too many with just one plate


    # Sigma based on ppn of (single) Normal within interval at bump
    # Ignore centering...
    # Big bump is 50%-100% of prob mass, centred around median of its own Normal
    # So Pnorm( next_quant, mean= bb, sd=sigma) - 0.5 = dquant * (1--2)
    # => qnorm( dquant*(1--2) + 0.5) = (next_quant-bb) / sd
    # => sd = (next_quant-bb) / qnorm

    sigma_hat <- dqd[ big_bumpi] /  qnorm( (bumpiq[2]-bumpiq[1])*(1:2) + 0.5)
    sigma_hat <- max( sigma_hat) # min seems to go too small...

    # Work with trimmed data to avoid crazy ****
    # Would be better to censor values above say 0.95
    # (ie not too high a quantile), but normalmixEM can't

    # Use bespoke constrained/"ML" fit using OO as well

    # These will be overwritten by lglk_shebang
    # FUCKING names
    mulo <- unname( if( big_bump > hix/2) big_bump/2 else big_bump)
    siglo <- sigma_hat
    muhi <- sighi <- NULL

    nobs <- gobs[i,]
    environment( lglk_shebang) <- environment()

    phat <- pA <- pB <- pC <- pO <- NULL # set by lglk_shebang
    penscale <- 0
    repeat{
      pstart <- unname( c(
          logit( 2* mulo / hix),
          log( siglo),
          logit( pambig[i,1]),
          logit( pambig[i,2] / (1-pambig[i,1]) ),
          # else no 5th param
          if( pambig[i,3]>0) logit( pambig[i,3] / (1-pambig[i,1]-pambig[i,2]) )
        ))

      if( !is.finite( lglk_shebang( pstart))) {
        siglo <- siglo * 2
      } else {
    break
      }
    } # until OK startval

    # Set penscale, following eg of est_ALF_ABCO
    testo <- numeric( 3)
    for( tempi in 1:3) {
      testo[ tempi] <- lglk_shebang( pstart+0.01*c( -1, 0, 1)[tempi])
    }
    openscale <- max( abs( diff( testo))) / 1e2

    # in case changing OOthresh & having to regenotype, or changed hix
    for( n_refits in 0:(max_refits-1)) {
      penscale <- openscale
      res <- nlminb( pstart, NEG( lglk_shebang),
          control=list( trace=6))

      # Should be able to reduce penscale now
      repeat{
        ophat <- phat
        penscale <- penscale / 10
        res <- nlminb( res$par, NEG( lglk_shebang),
            control=list( trace=6))
        if( rel.delta( phat, ophat) < 1e-2)
      break
      }

      # For now, organize the results to match 'normalmixEM'
      # ... so they can be extracted
      lglk_shebang( res$par) # ensure up-to-date
      res$lambda <- c( 2*pO*(1-pO), sqr( 1-pO))
      res$mu <- c( mulo, muhi)
      res$sigma <- c( siglo, sighi)
      # not comparable with normalmixEM since different data
      res$loglik <- (-res$objective)

      pbonzer[] <- c( pA, pB, pC, pO)

      if( n_refits == max_refits) { # ... then no need to check
        break
      }

      # May revise OOthresh and/or hix
      refit_needed <- FALSE # unless changed by checks below
      check_OOthresh <- pO > minpO_rethresh # otherwise not worth bothering with
      phix <- CDF( hix, muhi, sighi, lower=TRUE)
      # 0.6 since muhi is constrained <= hix; 0.98 for outliers
      shift_hix <- ! (phix %in.range% c( 0.6, 0.95))
      # Can still result in constraint being hit, but only
      # for crap locus I think

      if( check_OOthresh) {
        # Try to split "widest gap". Need a rough upper limit to start:
        # 2X as many in bump as at OO
        new_OOthresh[ i] <- find_new_OOthresh(pO, mulo, siglo, dat,
                                              m_for_rethresh)
        new_is_OO <- dat < new_OOthresh[i]
        refit_needed <- !all( new_is_OO == is_OO)
        if( refit_needed) { # because genotypes have changed
          scatn( 'OO changed from %i to %i: regenotyping %s', sum( is_OO),
                sum( new_is_OO), locus_name[ i])
          geno[,i] <- genoi <- geno_deambig_ABC( lociar[,i,,drop=FALSE],
              mht=mht[i],
              OOthresh_tc=new_OOthresh[i],
              het_cut=lociar@het_cut,
              tc_hist_pars= NULL,
              ppnA_hist_pars= NULL,
              return_what= 'just_geno')

          # Prepare for refit
          for( ig in genotypes_ambig) {
            nobs[ig] <- sum( genoi==ig)
          }

          # Total counts shouldn't change much, except tot[OO] still has
          # A+B counts, and may get regenoed to AAO which includes only A counts
          tot2[,i] <- make_tot2( cA[,i], cB[,i], cC[,i], genoi)
          dat <- tot2[,i]

          setup_counts_and_nhi()
        } # if genotypes change
      } # if  check_OOthresh

      if( shift_hix) {
        scatn( 'Changing the outlier limit; refitting')
        refit_needed <- TRUE
        hix <- qnorm( pmax( 0.6, pmin( 0.98, phix)), mean=muhi, sd=sighi)
        setup_counts_and_nhi()
      }

      if( !refit_needed) {
        break
      }

      pstart <- res$par  # start from same place, except...
      pstart[1] <- logit( mulo / (hix/2))
    } # for n_refits

    # Now choose cutoff for distinguishing ZO from ZZ where Z is Ref or
    # Alt (A or B) and O is "true" null ie not allele C
    # If allele C is present, will always be scored as ZO (or OO)

    # Has been fitted to TOTAL counts, so including allele C etc
    # Choose allele-specific cutoff to min errors based on already-calculated
    # pA etc, using *widths* of distros just calced here...
    # ... but only A & B alleles

    environment( diffprob) <- list2env( res) # and also variables defined here

    # *Conditional* error rate for XXO as XX or XO (does not apply if C present)
    best_cut <- c( A=0, B=0)
    snerr <- structure( rep( 0, 4), names=cq( AA2AO, AO2AA, BB2BO, BO2BB))
    for( al in cq( A, B)) {
      pal <- get( 'p' %&% al)
      best_cut[ al] <- if( diffprob( new_OOthresh[i], pX=pal, pO=pO) > 0)
        # did have OOthresh_tc instead here, though shouldn't matter much
          new_OOthresh[i]
        else # Paige's cutoff: starting value for root-finder
          find.root(diffprob, start=mulo*sqrt(2), step=mulo/6, pX=pal, pO=pO)

      # snerr: technically wrong if on OOthresh bdy, but small...
      snerr[gsub('X', al, 'XX2XO')] <- diffprob(best_cut[al], pX=1, pO=0)
      # NB sign !!!
      snerr[gsub('X', al, 'XO2XX')] <- (-diffprob(best_cut[al], pX=0, pO=1))
    }

    if( plots) {
      geno6_plot_counts()
    }

    # Obscure syntax, but it works
    this_res <- with( res, returnList(
        mu, sigma, lambda, loglik, best_cut, snerr, pbonzer))
    if( is.null( allres)) {
      allres <- data.frame(FOR(this_res, I( matrix( ., nrow=n_loci,
                           ncol=length(.), byrow=TRUE,
                           dimnames=list(NULL, names( .))))))
    }
    for( thing in names( this_res)) {
      allres[[thing]][i,] <- this_res[[thing]]
    }
  } # for loci

  if( !show_progress) {
    sink()
    unlink( tt)
    on.exit()
  }

  for( scalar_thing in which( sapply( allres, ncol)==1)) {
    allres[[ scalar_thing]] <- c( allres[[ scalar_thing]])
  }

  li$OOthresh_tc <- lociar@locinfo$OOthresh_tc <- new_OOthresh

  p6 <- with(allres, calc_g6probs(pbonzer[,'A'], pbonzer[,'B'],
                                  pbonzer[,'C'], snerr))

  lociar@locinfo <- cbind( li, allres[ cq( snerr, best_cut, pbonzer)])
  lociar@locinfo$perr <- p6@perr # gets messed up if in previous c() call
  return( lociar)
}
